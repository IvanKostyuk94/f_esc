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


def completed_haloIDs(sim_name, snap_name):
    """
    Return list of halos that have completed runs as indicated by CRASH check
    file.
    """

    def test_check(rdir):
        checkfile = _os.path.join(
            rdir, os.path.basename(exe_path())+".check"
        )
        if os.path.isfile(checkfile):
            with open(checkfile) as f:
                checkval = int(f.read())
            if checkval == 0:
                return True
        return False

    run_dirs = filter(
        os.path.isdir,
        _glob.glob(_os.path.join(
            rundir(sim_name, snap_name), "*"
        ))
    )
    completed = filter(
        lambda x: test_check(x),
        run_dirs
    )
    return list(map(halo_number, map(_os.path.basename, completed)))
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
        error_message = f'Length of halo_keys and configs has to be equal, received the shapes {len(halo_keys)} and {len(configs)} instead.' 
        raise ValueError(error_message)

    for key, config in enumerate(configs):
        print(f'Working on dictionary for {key}')
        
        sim = config.get_sim(simname)
        snaps = config.to_process[simname]
        halos[halo_keys[key]] = {}
        
        for snap in snaps:
            print(f'Working on snapshot {snap}')

            df = pd.read_pickle(config.selected_halo_df_file(simname, config.snap_name(snap)))
            df['csim_path', 0] = np.nan
            for i in range(config.n_rcs):
                df['f_esc', i] = np.nan

            # Begin insert for testing
            # counter = 0 
            # End insert for testing
            for ID in config.completed_haloIDs(simname, config.snap_name(snap)):
                # Begin insert for testing
                # counter += 1
                # if counter > 10:
                #    continue
                # End insert for testing

                rundir = config.rundir(simname, config.snap_name(snap), config.halo_name(ID))
                df.loc[ID, 'csim_path'] = rundir
                try:
                    csim = Sim(rundir)
                except:
                    print(f'An error occured in {rundir}')
                    continue
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
                    print(f'Located at: {rundir}')

            dic = {}

            dic['header'] = sim.snap_cat[snap].header
            dic['df'] = df
            finished_IDs = config.completed_haloIDs(simname, config.snap_name(snap))
            dic['IDs'] = np.array(finished_IDs)

            halos[halo_keys[key]][config.snap_name(snap)] = dic 

    return halos

# Construct a dataframe 
def construct_freq_dataframe(dictionary, configs, settings, name='freq_f_esc', fesc_galaxy=False):
    # Begin comment out to try different approach
    #store = pd.HDFStore(name,'w')
    # End comment out

    for setting, config in zip(settings, configs):
        print(f'Working on dataframe for setting {setting}')
        f_esc = []
        f_esc_0_2 = []
        halo_masses = []
        metal = []
        rel_star_mass = []
        rel_gas_mass = []
        rel_dust_mass = []
        all_IDs = []
        Q0 = []
        redshifts = []
        halo_radii = []
        Temperature = []
        xHII = []
        xHeII = []
        xHeIII = []
        grid_size = []
        bh_mass = []
        bh_growth = []
        sfr = []
        density = []
        clumping = []

        per_freq = []
        per_source = []
        emitted_photons = []
        escaped_photons = []
        frequencies = []
        n_iterations = []

        per_freq_0_2 = []
        per_source_0_2 = []
        emitted_photons_0_2 = []
        escaped_photons_0_2 = []
        frequencies_0_2 = []
        n_iterations_0_2 = []


        for i,snapshot in enumerate(dictionary[setting].keys()):
            print(f'Working an snapshot {snapshot}')

            IDs = dictionary[setting][snapshot]['IDs']
            
            f_esc_elements = []
            per_freq_elements = []
            per_source_elements = []
            emitted_photons_elements = []
            escaped_photons_elements = []
            frequencies_elements = []
            n_iterations_elements = []
            
            f_esc_elements_0_2 = []
            per_freq_elements_0_2 = []
            per_source_elements_0_2 = []
            emitted_photons_elements_0_2 = []
            escaped_photons_elements_0_2 = []
            frequencies_elements_0_2 = []
            n_iterations_elements_0_2 = []
            
            Temperature_elements = []
            xHII_elements = []
            xHeII_elements = []
            xHeIII_elements = []
            grid_size_elements = []
            density_elements = []
            clumping_elements = []
            
            input_df = dictionary[setting][snapshot]['df']

            # begin test
            # IDs = IDs[:9]
            # end test
            for ID in IDs:
                f_esc_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['cum'])
                per_freq_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['per_freq'])
                per_source_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['per_source'])
                emitted_photons_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['emitted_photons'])
                escaped_photons_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['escaped_photons'])
                frequencies_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['freqs'])
                n_iterations_elements.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['1.0e0']['n_iterations'])
                
                if fesc_galaxy:
                    f_esc_elements_0_2.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['2.0e-1']['cum'])
                    per_freq_elements_0_2.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['2.0e-1']['per_freq'])
                    per_source_elements_0_2.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['2.0e-1']['per_source'])
                    emitted_photons_elements_0_2.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['2.0e-1']['emitted_photons'])
                    escaped_photons_elements_0_2.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['2.0e-1']['escaped_photons'])
                    frequencies_elements_0_2.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['2.0e-1']['freqs'])
                    n_iterations_elements_0_2.append(input_df.loc[ID, ('f_esc',4)]['5.0e-2']['2.0e-1']['n_iterations'])

                average_quantities = get_average_quantities(ID, config, z[i])
                Temperature_elements.append(average_quantities[0])
                xHII_elements.append(average_quantities[1])
                xHeII_elements.append(average_quantities[2])
                xHeIII_elements.append(average_quantities[3])
                grid_size_elements.append(average_quantities[4])
                density_elements.append(average_quantities[5])
                clumping_elements.append(average_quantities[6])
            
            group_mass_elements = input_df.loc[IDs, ('GroupMass',0)]
            metal_elements = input_df.loc[IDs, ('GroupStarMetallicity', 0)]#/1e-3
            star_mass_elements = input_df.loc[IDs, ('GroupMassType', 4)]#/1e-3
            gas_elements = input_df.loc[IDs, ('gas_mass',0)]#/1e-1
            dust_mass_elements = input_df.loc[IDs, ('dust_mass')]#/1e-5
            Q0_elements = input_df.loc[IDs, ('Q_0',0)]#/1e53
            halo_radii_elements = input_df.loc[IDs, ('Group_R_Crit200',0)]
            bh_mass_elements = input_df.loc[IDs, ('GroupBHMass', 0)]
            bh_growth_elements = input_df.loc[IDs, ('GroupBHMdot', 0)]
            sfr_elements = input_df.loc[IDs, ('GroupSFR', 0)]

            f_esc.extend(f_esc_elements)
            per_freq.extend(per_freq_elements)
            per_source.extend(per_source_elements)
            emitted_photons.extend(emitted_photons_elements)
            escaped_photons.extend(escaped_photons_elements)
            frequencies.extend(frequencies_elements)
            n_iterations.extend(n_iterations_elements)

            f_esc_0_2.extend(f_esc_elements_0_2)
            per_freq_0_2.extend(per_freq_elements_0_2)
            per_source_0_2.extend(per_source_elements_0_2)
            emitted_photons_0_2.extend(emitted_photons_elements_0_2)
            escaped_photons_0_2.extend(escaped_photons_elements_0_2)
            frequencies_0_2.extend(frequencies_elements_0_2)
            n_iterations_0_2.extend(n_iterations_elements_0_2)

            halo_masses.extend(group_mass_elements)
            metal.extend(metal_elements)
            rel_star_mass.extend(star_mass_elements/group_mass_elements)
            rel_gas_mass.extend(gas_elements/group_mass_elements)
            rel_dust_mass.extend(dust_mass_elements/group_mass_elements)
            redshifts.extend(np.full(len(IDs),z[i]))
            all_IDs.extend(IDs)
            Q0.extend(Q0_elements)
            halo_radii.extend(halo_radii_elements)
            Temperature.extend(Temperature_elements)
            xHII.extend(xHII_elements)
            xHeII.extend(xHeII_elements)
            xHeIII.extend(xHeIII_elements)
            grid_size.extend(grid_size_elements)
            bh_mass.extend(bh_mass_elements)
            bh_growth.extend(bh_growth_elements)
            sfr.extend(sfr_elements)
            density.extend(density_elements)
            clumping.extend(clumping_elements)

        if fesc_galaxy:
            dataset = pd.DataFrame({'ID':all_IDs, 'z':redshifts, 
                                'HaloMass':halo_masses, 'Metallicity':metal, 
                                'FractionStars':rel_star_mass, 'FractionGas':rel_gas_mass,
                                'FractionDust':rel_dust_mass, 'Q0':Q0, 
                                'HaloRadii':halo_radii, 'f_esc':f_esc, 'f_esc_0_2':f_esc_0_2,
                                'Temperature': Temperature, 'xHII':xHII, 'xHeII':xHeII, 'xHeIII':xHeIII, 'GridSize':grid_size,
                                'BHMass':bh_mass, 'BHGrowth':bh_growth, 'SFR':sfr,
                                'density':density, 'clumping':clumping,
                                'per_freq':per_freq, 'per_source':per_source, 'emitted_photons':emitted_photons, 
                                'escaped_photons':escaped_photons, 'frequencies':frequencies, 'n_iterations':n_iterations,
                                'per_freq_0_2':per_freq_0_2, 'per_source_0_2':per_source_0_2, 'emitted_photons_0_2':emitted_photons_0_2, 
                                'escaped_photons_0_2':escaped_photons_0_2, 'frequencies_0_2':frequencies_0_2, 'n_iterations_0_2':n_iterations_0_2})
        else:
            dataset = pd.DataFrame({'ID':all_IDs, 'z':redshifts, 
                                'HaloMass':halo_masses, 'Metallicity':metal, 
                                'FractionStars':rel_star_mass, 'FractionGas':rel_gas_mass,
                                'FractionDust':rel_dust_mass, 'Q0':Q0, 
                                'HaloRadii':halo_radii, 'f_esc':f_esc,
                                'Temperature': Temperature, 'xHII':xHII, 'xHeII':xHeII, 'xHeIII':xHeIII, 'GridSize':grid_size,
                                'BHMass':bh_mass, 'BHGrowth':bh_growth, 'SFR':sfr,
                                'density':density, 'clumping':clumping,
                                'per_freq':per_freq, 'per_source':per_source, 'emitted_photons':emitted_photons, 
                                'escaped_photons':escaped_photons, 'frequencies':frequencies, 'n_iterations':n_iterations})
        # Set f_esc to float64 instead of df object (not done automatically for some reason)
        dataset = dataset.astype(dtype = {"f_esc":"float64"})
        if fesc_galaxy:
            dataset = dataset.astype(dtype = {'f_esc_0_2':'float64'})

        dataset.to_pickle(name)
    # Begin comment out to try different approach
        #store[setting] = dataset
    #store.close()
    # End comment out

def redshift_to_snap(redshift):
    correspondense = {6:'sn013', 8:'sn008', 10:'sn004'}
    return correspondense[redshift]

def get_simulation_path(halo_id, conf, redshift):
    snap = redshift_to_snap(redshift)
    conf_dir = os.path.join('/ptmp/mpa/mglatzle/TNG_f_esc', conf)
    simulation_path =  os.path.join(conf_dir, f'run/L35n2160TNG/{snap}/g{halo_id}/Output/phys_ic00_rt05.out')
    density_path = os.path.join(conf_dir, f'run/L35n2160TNG/{snap}/g{halo_id}/Input/dens_ic00.in')
    return simulation_path, density_path


def clumping_factor(density_map):
    volume = density_map.shape[0]**3
    C = np.sum(np.square(density_map))*volume/(np.sum(density_map)**2)
    return C

def get_average_quantities(halo_id, conf, redshift):
    path_sim, path_dens = get_simulation_path(halo_id, conf, redshift)
    dens = crashMemMap(path_dens, 'all')[0]
    average_quantities = []
    try:
        halo = crashMemMap(path_sim, 'all')
        for i in range(len(halo)):
            dens_weighted = halo[i]*dens 
            quantity = np.sum(dens_weighted)/np.sum(dens)
            average_quantities.append(quantity)
        average_quantities.append(int(halo[0].shape[0]))
        average_quantities.append(np.average(dens))
        average_quantities.append(clumping_factor(dens))
    except:
        print('Could not obtain average quantities for halo')
        print(halo_id)
        print(conf)
        print(redshift)
        print('-'*80)
        average_quantities = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan] 
    return average_quantities


def build_fid_df(simname, name='df_f_esc_freq.h5'):
    # load the config module for the fiducial 2 dust configuration
    spec_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2d/config.py")
    config_dust = importlib.util.module_from_spec(spec_dust)
    spec_dust.loader.exec_module(config_dust)

    # load the config module for the fiducial 2 no dust configuration
    spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2/config.py")
    config_no_dust = importlib.util.module_from_spec(spec_no_dust)
    spec_no_dust.loader.exec_module(config_no_dust)

    halos = construct_halo_dict(simname, [config_no_dust, config_dust], with_dust=True)
    construct_freq_dataframe(dictionary=halos, name = name, configs = ['fid2', 'fid2d'], settings=['no_dust', 'dust'])
    return

def build_no_dust_df(simname, name='df_no_dust.h5'):
    # load the config module for the fiducial 2 no dust configuration
    spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2/config.py")
    config_no_dust = importlib.util.module_from_spec(spec_no_dust)
    spec_no_dust.loader.exec_module(config_no_dust)

    print('Building halo dictionary')
    halos = construct_halo_dict(simname, configs = [config_no_dust], with_dust=False)
    print('Finished buiding halo dictionary')
    print('Constructing halo dataframe')
    # Note here 'configs' corresponds to the names of the runs and not to the config object as above
    construct_freq_dataframe(dictionary=halos, name = name, configs = ['fid2'], settings=['no_dust'], fesc_galaxy=True)
    return

def build_div_esc(simname, name='var_esc.h5'):
    spec_0_3 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_3e-1/config.py")
    spec_0_5 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_5e-1/config.py")
    spec_0_7 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_7e-1/config.py")
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
    configs = ['esc_3e-1','esc_5e-1','esc_7e-1','full_esc']
    
    construct_freq_dataframe(dictionary=halos, name = name, configs=configs, settings=halo_keys)
    return


def build_full_esc_df(simname, name='df_full_esc_freq.h5'):
    spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/full_esc/config.py")
    config_no_dust = importlib.util.module_from_spec(spec_no_dust)
    spec_no_dust.loader.exec_module(config_no_dust)

    halos = construct_halo_dict(simname, config_no_dust) #config_dust
    construct_freq_dataframe(dictionary=halos, name = name, settings=['no_dust'],  configs = ['full_esc'])
    return

if __name__ == "__main__":
    # Set some global variables such as the redshifts of the snapshots and the simulation name
    z = [6,8,10]
    simname = 'L35n2160TNG'
    h = 6.62607004e-34
    e = 1.60217662e-19
    j_to_erg = 1e7
    norm = 1e52

    build_no_dust_df(simname)#, name='TestyMcTestface')
    #build_div_esc(simname)
    #build_full_esc_df(simname)
