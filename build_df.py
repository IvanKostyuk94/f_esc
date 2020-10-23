import matplotlib.pyplot as plt
import importlib.util
import pandas as pd
import numpy as np
import h5py
import tables
import os
from crashpy.dataclasses.simulation import LoadedSimulation as Sim
from crashpy.utilities import crashMemMap
from get_spectrum import get_spectra

"""
Load data from h5py file into a dict.

This assumes all datasets contained in the h5py file can be converted to
numpy arrays.

Parameters
----------
f : h5py.File
    h5py file object to load data from.

Returns
-------
dd : dict
    Dictionary of all group and datasets contained in 'f'. Hierarchical
    structure is kept identical.
"""
def h5toDict(f):

    d = {}

    def visitor(name, node):
        if isinstance(node, h5py.Dataset):
            d[name] = np.array(node)
        return 
    
    f.visititems(visitor)

    # restore hierarchical structure as the visitor flattens it
    dd = {}
    for key, val in d.items():
        keys = key.split('/')
        ddd = dd
        for k in keys[:-1]:
            if k in ddd:
                ddd = ddd[k]
            else:
                ddd[k] = {}
                ddd = ddd[k]
        ddd[keys[-1]] = val

    return dd

# dict of dicts, first level is TNG simulation, second level is simulation
# snapshot. For each snapshot, load its header (key 'header') and a halo
# dataframe (key 'df').
def construct_halo_dict(simname, configs, halo_keys = ['no_dust'], with_dust=False):
    halos = {}
    if not with_dust:
        halo_keys = halo_keys
        configs = configs
        
    else:
        halo_keys = ['no_dust', 'dust']
        configs = configs

    if not (len(halo_keys)==len(configs)):
        raise ValueError(f'Length of halo_keys and configs has to be equal, received the shapes {len(halo_keys)} and {len(configs)} instead.')

    for key, config in enumerate(configs):
        
        sim = config.get_sim(simname)
        snaps = config.to_process[simname]
        halos[halo_keys[key]] = {}
        
        for snap in snaps:
            df = pd.read_pickle(config.selected_halo_df_file(simname, config.snap_name(snap)))
            df['csim_path', 0] = np.nan
            for i in range(config.n_rcs):
                df['f_esc', i] = np.nan

            # Remove halos with defective is_binary test
            defect_IDs = []
            for ID in config.completed_haloIDs(simname, config.snap_name(snap)):
                rundir = config.rundir(simname, config.snap_name(snap), config.halo_name(ID))
                df.loc[ID, 'csim_path'] = rundir
                try:
                    csim = Sim(rundir)
                except:
                    print(f'An error occured in {rundir}')
                    defect_IDs.append(ID)
                ls = []
                for i, pf in enumerate(csim.getAllphysfiles()):
                    try:
                        with h5py.File(config.f_esc_file(pf),'r') as f:
                            fesc = h5toDict(f)
                        ls.append(fesc)
                    except:
                        pass
                try:
                    df.loc[ID, 'f_esc'] = np.array(ls)
                except:
                    print(f'Could not load the output of {ID}')
                    defect_IDs.append(ID)

            dic = {}

            dic['header'] = sim.snap_cat[snap].header
            dic['df'] = df
            finished_IDs = config.completed_haloIDs(simname, config.snap_name(snap))
            for ID in defect_IDs:
                finished_IDs.remove(ID)
            dic['IDs'] = np.array(finished_IDs)

            halos[halo_keys[key]][config.snap_name(snap)] = dic 

    return halos

# Construct a dataframe 
def construct_freq_dataframe(dictionary, configs, settings=['dust', 'no_dust'], name='freq_f_esc'):
    store = pd.HDFStore(name,'w')

    for setting, config in zip(settings, configs):
        f_esc = []
        halo_masses = []
        metal = []
        rel_star_mass = []
        rel_gas_mass = []
        rel_dust_mass = []
        all_IDs = []
        Q0 = []
        a_star = []
        redshifts = []
        halo_radii = []
        Temperature = []
        xHII = []
        xHeII = []
        xHeIII = []

        per_freq = []
        per_source = []
        emitted_photons = []
        escaped_photons = []
        frequencies = []
        n_iterations = []


        for i,snapshot in enumerate(dictionary[setting].keys()):
            IDs = dictionary[setting][snapshot]['IDs']
            
            f_esc_elements = []
            per_freq_elements = []
            per_source_elements = []
            emitted_photons_elements = []
            escaped_photons_elements = []
            frequencies_elements = []
            n_iterations_elements = []
            Temperature_elements = []
            xHII_elements = []
            xHeII_elements = []
            xHeIII_elements = []

            input_df = dictionary[setting][snapshot]['df']
            for ID in IDs:
                f_esc_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['cum'])
                per_freq_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['per_freq'])
                per_source_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['per_source'])
                emitted_photons_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['emitted_photons'])
                escaped_photons_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['escaped_photons'])
                frequencies_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['freqs'])
                n_iterations_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['n_iterations'])

                average_quantities = get_average_quantities(ID, config, z[i])
                Temperature_elements.append(average_quantities[0])
                xHII_elements.append(average_quantities[1])
                xHeII_elements.append(average_quantities[2])
                xHeIII_elements.append(average_quantities[3])
            
            group_mass_elements = input_df.loc[IDs, ('GroupMass',0)]
            metal_elements = input_df.loc[IDs, ('GroupStarMetallicity', 0)]#/1e-3
            star_mass_elements = input_df.loc[IDs, ('GroupMassType', 4)]#/1e-3
            gas_elements = input_df.loc[IDs, ('gas_mass',0)]#/1e-1
            dust_mass_elements = input_df.loc[IDs, ('dust_mass',0)]#/1e-5
            Q0_elements = input_df.loc[IDs, ('Q_0',0)]#/1e53
            a_star_elements = input_df.loc[IDs, ('a_star',0)]#/1e8
            halo_radii_elements = input_df.loc[IDs, ('Group_R_Crit200',0)]

            f_esc.extend(f_esc_elements)
            per_freq.extend(per_freq_elements)
            per_source.extend(per_freq_elements)
            emitted_photons.extend(emitted_photons_elements)
            escaped_photons.extend(escaped_photons_elements)
            halo_masses.extend(group_mass_elements)
            metal.extend(metal_elements)
            rel_star_mass.extend(star_mass_elements/group_mass_elements)
            rel_gas_mass.extend(gas_elements/group_mass_elements)
            rel_dust_mass.extend(dust_mass_elements/group_mass_elements)
            redshifts.extend(np.full(len(IDs),z[i]))
            all_IDs.extend(IDs)
            Q0.extend(Q0_elements)
            a_star.extend(a_star_elements)
            frequencies.extend(frequencies_elements)
            n_iterations.extend(n_iterations_elements)
            halo_radii.extend(halo_radii_elements)
            Temperature.extend(Temperature_elements)
            xHII.extend(xHII_elements)
            xHeII.extend(xHeII_elements)
            xHeIII.extend(xHeIII_elements)

            
        dataset = pd.DataFrame({'ID':all_IDs, 'z':redshifts, 
                            'HaloMass':halo_masses, 'Metallicity':metal, 
                            'FractionStars':rel_star_mass, 'FractionGas':rel_gas_mass,
                            'FractionDust':rel_dust_mass, 'Q0':Q0, 'aStar':a_star, 
                            'HaloRadii':halo_radii, 'f_esc':f_esc,
                            'Temperature': Temperature, 'xHII':xHII, 'xHeII':xHeII, 'xHeIII':xHeIII,
                            'per_freq':per_freq, 'per_source':per_source, 'emitted_photons':emitted_photons, 
                            'escaped_photons':escaped_photons, 'frequencies':frequencies, 'n_iterations':n_iterations})
        # Set f_esc to float64 instead of df object (not done automatically for some reason)
        dataset = dataset.astype(dtype= {"f_esc":"float64"})
        store[setting] = dataset
    store.close()

    return

def redshift_to_snap(redshift):
    correspondense = {6:'sn013', 8:'sn008', 10:'sn004'}
    return correspondense[redshift]

def get_simulation_path(halo_id, conf, redshift):
    snap = redshift_to_snap(redshift)
    conf_dir = os.path.join('/ptmp/mpa/mglatzle/TNG_f_esc', conf)
    simulation_path =  os.path.join(conf_dir, f'run/L35n2160TNG/{snap}/g{halo_id}/Output/phys_ic00_rt05.out')
    return simulation_path

def get_average_quantities(halo_id, conf, redshift):
    path = get_simulation_path(halo_id, conf, redshift)
    average_quantities = []
    try:
        halo = crashMemMap(path, 'all')
        for i in range(4):
            quantity = np.average(halo[i])
            average_quantities.append(quantity)
    except:
        average_quantities = [np.nan, np.nan, np.nan, np.nan] 
    return average_quantities


def build_fid_df(simname):
    # load the config module for the fiducial 2 dust configuration
    spec_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2d/config.py")
    config_dust = importlib.util.module_from_spec(spec_dust)
    spec_dust.loader.exec_module(config_dust)

    # load the config module for the fiducial 2 no dust configuration
    spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2/config.py")
    config_no_dust = importlib.util.module_from_spec(spec_no_dust)
    spec_no_dust.loader.exec_module(config_no_dust)

    halos = construct_halo_dict(simname, [config_no_dust, config_dust])
    construct_freq_dataframe(dictionary=halos, name = 'df_f_esc_freq.h5', configs = ['fid2', 'fid2d'])
    return

def build_div_esc(simname):
    spec_0_3 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_3_e-1/config.py")
    spec_0_5 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_5_e-1/config.py")
    spec_0_7 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_7_e-1/config.py")
    spec_1_0 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/full_esc/config.py")
    
    config_0_3 = importlib.util.module_from_spec(spec_0_3)
    config_0_5 = importlib.util.module_from_spec(spec_0_5)
    config_0_7 = importlib.util.module_from_spec(spec_0_7)
    config_1_0 = importlib.util.module_from_spec(spec_1_0)
    
    spec_0_3.loader.exec_module(config_0_3)
    spec_0_5.loader.exec_module(config_0_5)
    spec_0_7.loader.exec_module(config_0_7)
    spec_1_0.loader.exec_module(config_1_0)

    halo_keys = ['0.3','0.5','0.7','1.0']
    configs = [config_0_3, config_0_5, config_0_7, config_1_0]

    halos = construct_halo_dict(simname, configs=configs, halo_keys = halo_keys)
    configs = ['esc_3_e-1','esc_5_e-1','esc_7_e-1','full_esc']
    
    construct_freq_dataframe(dictionary=halos, configs=configs, settings=halo_keys)
    return


def build_full_esc_df(simname):
    spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/full_esc/config.py")
    config_no_dust = importlib.util.module_from_spec(spec_no_dust)
    spec_no_dust.loader.exec_module(config_no_dust)

    halos = construct_halo_dict(simname, config_no_dust) #config_dust
    construct_freq_dataframe(dictionary=halos, name = 'df_full_esc_freq.h5', settings=['no_dust'],  configs = ['full_esc'])
    return

if __name__ == "__main__":
    # Set some global variables such as the redshifts of the snapshots and the simulation name
    z = [6,8,10]
    simname = 'L35n2160TNG'
    h = 6.62607004e-34
    e = 1.60217662e-19
    j_to_erg = 1e7
    norm = 1e52

    #build_fid_df(simname)
    build_div_esc(simname)
    #build_full_esc_df(simname)
