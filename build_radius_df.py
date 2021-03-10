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
    f_esc_radii = {0.2:{}, 0.3:{}, 0.4:{}, 0.5:{}, 0.6:{}, 0.7:{}, 0.8:{}, 0.9:{}, 1.0:{}, 1.1:{}, 1.2:{}, 1.3:{}, 1.4:{}, 1.5:{}}
    f_esc_notation = {0.2:'2.0e-1', 0.3:'3.0e-1', 0.4:'4.0e-1', 0.5:'5.0e-1', 0.6:'6.0e-1', 0.7:'7.0e-1', 
    0.8:'8.0e-1', 0.9:'9.0e-1', 1.0:'1.0e0', 1.1:'1.1e0', 1.2:'1.2e0', 1.3:'1.3e0', 1.4:'1.4e0', 1.5:'1.5e0'}
    
    print(f'Working on dataframe for setting {config}')
    for key in f_esc_radii:
        f_esc_radii[key]['f_esc'] = []
        f_esc_radii[key]['per_freq'] = []
        f_esc_radii[key]['per_source'] = []
        f_esc_radii[key]['emitted_photons'] = []
        f_esc_radii[key]['escaped_photons'] = []
        f_esc_radii[key]['frequencies'] = []
        f_esc_radii[key]['n_iterations'] = []

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


    for i,snapshot in enumerate(dictionary.keys()):
        print(f'Working an snapshot {snapshot}')

        IDs = dictionary[snapshot]['IDs']

        f_esc_radii_elements = {0.2:{}, 0.3:{}, 0.4:{}, 0.5:{}, 0.6:{}, 0.7:{}, 0.8:{}, 0.9:{}, 1.0:{}, 1.1:{}, 1.2:{}, 1.3:{}, 1.4:{}, 1.5:{}}
        for key in f_esc_radii_elements:
            f_esc_radii_elements[key]['f_esc_elements'] = []
            f_esc_radii_elements[key]['per_freq_elements'] = []
            f_esc_radii_elements[key]['per_source_elements'] = []
            f_esc_radii_elements[key]['emitted_photons_elements'] = []
            f_esc_radii_elements[key]['escaped_photons_elements'] = []
            f_esc_radii_elements[key]['frequencies_elements'] = []
            f_esc_radii_elements[key]['n_iterations_elements'] = []
        
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
            for key in f_esc_radii_elements:
                f_esc_radii_elements[key]['f_esc_elements'].append(input_df.loc[ID, ('f_esc',4)]['5.0e-2'][f_esc_notation[key]]['cum'])
                f_esc_radii_elements[key]['per_freq_elements'].append(input_df.loc[ID, ('f_esc',4)]['5.0e-2'][f_esc_notation[key]]['per_freq'])
                f_esc_radii_elements[key]['per_source_elements'].append(input_df.loc[ID, ('f_esc',4)]['5.0e-2'][f_esc_notation[key]]['per_source'])
                f_esc_radii_elements[key]['emitted_photons_elements'].append(input_df.loc[ID, ('f_esc',4)]['5.0e-2'][f_esc_notation[key]]['emitted_photons'])
                f_esc_radii_elements[key]['escaped_photons_elements'].append(input_df.loc[ID, ('f_esc',4)]['5.0e-2'][f_esc_notation[key]]['escaped_photons'])
                f_esc_radii_elements[key]['frequencies_elements'].append(input_df.loc[ID, ('f_esc',4)]['5.0e-2'][f_esc_notation[key]]['freqs'])
                f_esc_radii_elements[key]['n_iterations_elements'].append(input_df.loc[ID, ('f_esc',4)]['5.0e-2'][f_esc_notation[key]]['n_iterations'])

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

        for key in f_esc_radii:
            f_esc_radii[key]['f_esc'].extend(f_esc_radii_elements[key]['f_esc_elements'])
            f_esc_radii[key]['per_freq'].extend(f_esc_radii_elements[key]['per_freq_elements'])
            f_esc_radii[key]['per_source'].extend(f_esc_radii_elements[key]['per_source_elements'])
            f_esc_radii[key]['emitted_photons'].extend(f_esc_radii_elements[key]['emitted_photons_elements'])
            f_esc_radii[key]['escaped_photons'].extend(f_esc_radii_elements[key]['escaped_photons_elements'])
            f_esc_radii[key]['frequencies'].extend(f_esc_radii_elements[key]['frequencies_elements'])
            f_esc_radii[key]['n_iterations'].extend(f_esc_radii_elements[key]['n_iterations_elements'])

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

    try:
        dataset_dict = {('ID', 0):all_IDs, ('z', 0):redshifts, 
                            ('HaloMass', 0):halo_masses, ('Metallicity', 0):metal, ('GasMetallicity', 0):gas_metal, 
                            ('FractionStars', 0):rel_star_mass, ('FractionGas', 0):rel_gas_mass,
                            ('FractionDust', 0):rel_dust_mass, ('Q0', 0):Q0, 
                            ('HaloRadii', 0):halo_radii, 
                            ('Temperature', 0): Temperature, ('xHII', 0):xHII, ('xHeII', 0):xHeII, ('xHeIII', 0):xHeIII, ('GridSize', 0):grid_size,
                            ('BHMass', 0):bh_mass, ('BHGrowth', 0):bh_growth, ('SFR', 0):sfr,
                            ('density', 0):density, ('clumping', 0):clumping}
        for key in f_esc_radii:
            dataset_dict[('f_esc', key)] = f_esc_radii[key]['f_esc']
            dataset_dict[('per_freq', key)] = f_esc_radii[key]['per_freq']
            dataset_dict[('per_source', key)] = f_esc_radii[key]['per_source']
            dataset_dict[('emitted_photons', key)] = f_esc_radii[key]['emitted_photons']
            dataset_dict[('escaped_photons', key)] = f_esc_radii[key]['escaped_photons']
            dataset_dict[('frequencies', key)] = f_esc_radii[key]['frequencies']
            dataset_dict[('n_iterations', key)] = f_esc_radii[key]['n_iterations']
        
        dataset = pd.DataFrame(dataset_dict)
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

        exit()


    # Set f_esc to float64 instead of df object (not done automatically for some reason)
    for key in f_esc_radii:
        dataset = dataset.astype(dtype = {('f_esc',key):"float64"})

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
    all_runs_tng = ['fid2', 'fid2d', 'esc_3e-1', 'esc_5e-1', 'esc_7e-1', 'full_esc', 'large_radii', 'numerical_1e4', 'numerical_1e6', 'TNG50_2', 'TNG50_3'] 

    simname_tng2 = 'L35n1080TNG'
    all_runs_tng = ['TNG50_2'] 

    simname_tng3 = 'L35n540TNG'
    all_runs_tng = ['TNG50_3'] 

    build_df(run_names=['large_radii'], filename=None, simname=simname_tng)