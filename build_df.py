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
def construct_halo_dict(simname, config):
    halos = {}
    print(f'Working on dictionary')
    
    sim = config.get_sim(simname)
    snaps = config.to_process[simname]
    
    for snap in snaps:
        print(f'Working on snapshot {snap}')

        df = pd.read_pickle(config.selected_halo_df_file(simname, config.snap_name(snap)))
        df['csim_path', 0] = np.nan
        for i in range(config.n_rcs):
            df['f_esc', i] = np.nan

        # Begin insert for testing
        # counter = 0 
        # End insert for testing
        
        # Begin patch only consider IDs actually included
        IDs_used = []
        # End patch
        
        completed_IDs = config.completed_haloIDs(simname, config.snap_name(snap))
        
        for ID in completed_IDs:
            # Begin insert for testing
            # counter += 1
            # if counter > 10:
            #    continue
            # End insert for testing
            # Begin patch only consider IDs actually included
            IDs_used.append(ID)
            # End patch

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
                # This is a temporary patch since not all halo IDs have the f_esc calculated yet
                IDs_used.remove(ID)
                #completed_IDs.remove(ID)
                print(f'Could not load the output of {ID}')
                print(f'Located at: {rundir}')

        dic = {}

        dic['header'] = sim.snap_cat[snap].header
        dic['df'] = df
        #finished_IDs = completed_IDs
        # This is only used until f_esc is finished
        finished_IDs = IDs_used
        dic['IDs'] = np.array(finished_IDs)

        halos[config.snap_name(snap)] = dic 

    return halos

# Construct a dataframe 
def construct_freq_dataframe(dictionary, config, name, simname, fesc_galaxy=False):
    # Be careful with the hard coded redshifts 
    z = [6,8,10]
    
    print(f'Working on dataframe for setting {config}')
    
    f_esc = []
    f_esc_0_2 = []
    halo_masses = []
    metal = []
    gas_metal = []
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


    for i,snapshot in enumerate(dictionary.keys()):
        print(f'Working an snapshot {snapshot}')

        IDs = dictionary[snapshot]['IDs']
        
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
        
        input_df = dictionary[snapshot]['df']

        # begin test
        # IDs = IDs[:9]
        # end test
        for ID in IDs:
            try:
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

            except:
                f_esc_elements.append(np.NaN)
                per_freq_elements.append(np.NaN)
                per_source_elements.append(np.NaN)
                emitted_photons_elements.append(np.NaN)
                escaped_photons_elements.append(np.NaN)
                frequencies_elements.append(np.NaN)
                n_iterations_elements.append(np.NaN)

                if fesc_galaxy:
                    f_esc_elements_0_2.append(np.NaN)
                    per_freq_elements_0_2.append(np.NaN)
                    per_source_elements_0_2.append(np.NaN)
                    emitted_photons_elements_0_2.append(np.NaN)
                    escaped_photons_elements_0_2.append(np.NaN)
                    frequencies_elements_0_2.append(np.NaN)
                    n_iterations_elements_0_2.append(np.NaN)


            average_quantities = get_average_quantities(ID, config, z[i], simname)
            Temperature_elements.append(average_quantities[0])
            xHII_elements.append(average_quantities[1])
            xHeII_elements.append(average_quantities[2])
            xHeIII_elements.append(average_quantities[3])
            grid_size_elements.append(average_quantities[4])
            density_elements.append(average_quantities[5])
            clumping_elements.append(average_quantities[6])
        
        group_mass_elements = input_df.loc[IDs, ('GroupMass',0)]
        metal_elements = input_df.loc[IDs, ('GroupStarMetallicity', 0)]#/1e-3
        metal_gas_elements = input_df.loc[IDs, ('GroupGasMetallicity', 0)]
        star_mass_elements = input_df.loc[IDs, ('GroupMassType', 4)]#/1e-3
        gas_elements = input_df.loc[IDs, ('gas_mass',0)]#/1e-1
        dust_mass_elements = input_df.loc[IDs, ('dust_mass', 0)]#/1e-5
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
        gas_metal.extend(metal_gas_elements)
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
            'HaloMass':halo_masses, 'Metallicity':metal, 'GasMetallicity':gas_metal,
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
        try:
            dataset = pd.DataFrame({'ID':all_IDs, 'z':redshifts, 
                'HaloMass':halo_masses, 'Metallicity':metal, 'GasMetallicity':gas_metal, 
                                'FractionStars':rel_star_mass, 'FractionGas':rel_gas_mass,
                                'FractionDust':rel_dust_mass, 'Q0':Q0, 
                                'HaloRadii':halo_radii, 'f_esc':f_esc,
                                'Temperature': Temperature, 'xHII':xHII, 'xHeII':xHeII, 'xHeIII':xHeIII, 'GridSize':grid_size,
                                'BHMass':bh_mass, 'BHGrowth':bh_growth, 'SFR':sfr,
                                'density':density, 'clumping':clumping,
                                'per_freq':per_freq, 'per_source':per_source, 'emitted_photons':emitted_photons, 
                                'escaped_photons':escaped_photons, 'frequencies':frequencies, 'n_iterations':n_iterations})
        except:
            print(f'ID: {len(all_IDs)}')
            print(f'z: {len(redshifts)}')

            print(f'HaloMass: {len(halo_masses)}')
            print(f'Metalicity: {len(metal)}')
            
            print(f'FractionStars: {len(rel_star_mass)}')
            print(f'FractionGas: {len(rel_gas_mass)}')
            
            print(f'FractionDust: {len(rel_dust_mass)}')
            print(f'Q0: {len(Q0)}')

            print(f'HaloRadii: {len(halo_radii)}')
            print(f'f_esc: {len(f_esc)}')
            
            print(f'Temperature: {len(Temperature)}')
            print(f'xHII: {len(xHII)}')
            print(f'xHeII: {len(xHeII)}')
            print(f'xHeIII: {len(xHeIII)}')
            print(f'GridSize:{len(grid_size)}')

            print(f'BHMass: {len(bh_mass)}')
            print(f'BHGrowth: {len(bh_growth)}')
            print(f'SFR: {len(sfr)}')
            
            print(f'density: {len(density)}')
            print(f'clumping: {len(clumping)}')

            print(f'per_freq: {len(per_freq)}')
            print(f'per_source: {len(per_source)}')
            print(f'emitted_photons: {len(emitted_photons)}')
            
            print(f'escaped_photons: {len(escaped_photons)}')
            print(f'frequencies: {len(frequencies)}')
            print(f'n_iterations: {len(n_iterations)}')
            exit()


    # Set f_esc to float64 instead of df object (not done automatically for some reason)
    dataset = dataset.astype(dtype = {"f_esc":"float64"})
    if fesc_galaxy:
        dataset = dataset.astype(dtype = {'f_esc_0_2':'float64'})

    dataset.to_pickle(name)
    return


def redshift_to_snap(redshift):
    correspondense = {6:'sn013', 8:'sn008', 10:'sn004'}
    return correspondense[redshift]

def get_simulation_path(halo_id, conf, redshift, simname):
    snap = redshift_to_snap(redshift)
    conf_dir = os.path.join('/ptmp/mpa/mglatzle/TNG_f_esc', conf)
    simulation_path =  os.path.join(conf_dir, f'run/{simname}/{snap}/g{halo_id}/Output/phys_ic00_rt05.out')
    density_path = os.path.join(conf_dir, f'run/{simname}/{snap}/g{halo_id}/Input/dens_ic00.in')
    return simulation_path, density_path


def clumping_factor(density_map):
    volume = density_map.shape[0]**3
    C = np.sum(np.square(density_map))*volume/(np.sum(density_map)**2)
    return C

def get_average_quantities(halo_id, conf, redshift, simname):
    path_sim, path_dens = get_simulation_path(halo_id, conf, redshift, simname)
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


def build_df(run_names=['fid2'], filename=None, simname='L35n2160TNG', fesc_galaxy=False):
    for run_name in run_names:
        df_dir = '/u/ivkos/analysis/dfs'
        if filename == None:
            df_dir = '/u/ivkos/analysis/dfs'
            outputpath = os.path.join(df_dir, f'{run_name}.pickle')
        else:
            outputpath = os.path.join(df_dir, filename)
        # load the config module for the fiducial 2 no dust configuration
        spec = importlib.util.spec_from_file_location("module.name",f"/freya/ptmp/mpa/mglatzle/TNG_f_esc/{run_name}/config.py")
        config = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config)

        print('Building halo dictionary')
        halos = construct_halo_dict(simname, config=config)
        print('Finished buiding halo dictionary')
        print('Constructing halo dataframe')
        # Note here 'configs' corresponds to the names of the runs and not to the config object as above
        construct_freq_dataframe(dictionary=halos, name=outputpath, config=run_name, simname=simname, fesc_galaxy=fesc_galaxy)
    return

if __name__ == "__main__":

    simname_tng = 'L35n2160TNG'
    all_runs_tng = ['fid2', 'fid2d', 'esc_3e-1', 'esc_5e-1', 'esc_7e-1', 'full_esc', 'large_radii', 'numerical_1e4', 'numerical_1e6', 'new_numerical_1e6', 'new_numerical_5e6', 'TNG50_2', 'TNG50_3'] 

    simname_tng2 = 'L35n1080TNG'
    all_runs_tng = ['TNG50_2'] 

    simname_tng3 = 'L35n540TNG'
    all_runs_tng = ['TNG50_3'] 

    build_df(run_names=['full_esc'], filename=None, simname=simname_tng, fesc_galaxy=True)


# Currently not working
# def build_fid_df(simname, name='df_f_esc_freq.h5'):
#     # load the config module for the fiducial 2 dust configuration
#     spec_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2d/config.py")
#     config_dust = importlib.util.module_from_spec(spec_dust)
#     spec_dust.loader.exec_module(config_dust)

#     # load the config module for the fiducial 2 no dust configuration
#     spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2/config.py")
#     config_no_dust = importlib.util.module_from_spec(spec_no_dust)
#     spec_no_dust.loader.exec_module(config_no_dust)

#     halos = construct_halo_dict(simname, [config_no_dust, config_dust], with_dust=True)
#     construct_freq_dataframe(dictionary=halos, name = name, configs = ['fid2', 'fid2d'], settings=['no_dust', 'dust'])
#     return


# def build_div_esc(simname, name='var_esc.h5'):
#     spec_0_3 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_3e-1/config.py")
#     spec_0_5 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_5e-1/config.py")
#     spec_0_7 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/esc_7e-1/config.py")
#     spec_1_0 = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/full_esc/config.py")
    
#     config_0_3 = importlib.util.module_from_spec(spec_0_3)
#     config_0_5 = importlib.util.module_from_spec(spec_0_5)
#     config_0_7 = importlib.util.module_from_spec(spec_0_7)
#     config_1_0 = importlib.util.module_from_spec(spec_1_0)
    
#     spec_0_3.loader.exec_module(config_0_3)
#     spec_0_5.loader.exec_module(config_0_5)
#     spec_0_7.loader.exec_module(config_0_7)
#     spec_1_0.loader.exec_module(config_1_0)

#     halo_keys = ['0.3','0.5','0.7','1.0']
#     configs = [config_0_3, config_0_5, config_0_7, config_1_0]

#     halos = construct_halo_dict(simname, configs=configs, halo_keys = halo_keys)
#     configs = ['esc_3e-1','esc_5e-1','esc_7e-1','full_esc']
    
#     construct_freq_dataframe(dictionary=halos, name = name, configs=configs, settings=halo_keys)
#     return


# def build_full_esc_df(simname, name='df_full_esc_freq.h5'):
#     spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/full_esc/config.py")
#     config_no_dust = importlib.util.module_from_spec(spec_no_dust)
#     spec_no_dust.loader.exec_module(config_no_dust)

#     halos = construct_halo_dict(simname, config_no_dust) #config_dust
#     construct_freq_dataframe(dictionary=halos, name = name, settings=['no_dust'],  configs = ['full_esc'])
#     return

