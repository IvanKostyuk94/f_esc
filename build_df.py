import matplotlib.pyplot as plt
import importlib.util
import pandas as pd
import numpy as np
import h5py
import tables
from crashpy.dataclasses.simulation import LoadedSimulation as Sim
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
def construct_halo_dict(simname, config_dust, config_no_dust):
    halos = {}
    halo_keys = ['dust', 'no_dust']
    for key, config in enumerate([config_dust, config_no_dust]):
        
        sim = config.get_sim(simname)
        snaps = config.to_process[simname]
        halos[halo_keys[key]] = {}
        for snap in snaps:

            df = pd.read_pickle(config.selected_halo_df_file(simname, config.snap_name(snap)))
            df['csim_path', 0] = np.nan
            for i in range(config.n_rcs):
                df['f_esc', i] = np.nan
            for ID in config.completed_haloIDs(simname, config.snap_name(snap)):
                rundir = config.rundir(simname, config.snap_name(snap), config.halo_name(ID))
                df.loc[ID, 'csim_path'] = rundir
                csim = Sim(rundir)
                ls = []
                for i, pf in enumerate(csim.getAllphysfiles()):
                    with h5py.File(config.f_esc_file(pf),'r') as f:
                        fesc = h5toDict(f)
                    ls.append(fesc)
                df.loc[ID, 'f_esc'] = np.array(ls)

            dic = {}

            dic['header'] = sim.snap_cat[snap].header
            dic['df'] = df
            dic['IDs'] = np.array(config.completed_haloIDs(simname, config.snap_name(snap)))

            halos[halo_keys[key]][config.snap_name(snap)] = dic 

    return halos

# Construct a dataframe 
def construct_dataframe(dictionary, settings=['dust', 'no_dust'], name='f_esc.h5'):
    store = pd.HDFStore(name,'w')

    for setting in settings:
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
        for i,snapshot in enumerate(dictionary[setting].keys()):
            IDs = dictionary[setting][snapshot]['IDs']
            
            f_esc_elements = []
            input_df = dictionary[setting][snapshot]['df']
            for ID in IDs:
                f_esc_elements.append(input_df.loc[ID, ('f_esc',0)]['5.0e-2']['1.0e0']['cum'])
            
            group_mass_elements = input_df.loc[IDs, ('GroupMass',0)]
            metal_elements = input_df.loc[IDs, ('GroupStarMetallicity', 0)]#/1e-3
            star_mass_elements = input_df.loc[IDs, ('M_star', 0)]#/1e-3
            gas_elements = input_df.loc[IDs, ('gas_mass',0)]#/1e-1
            dust_mass_elements = input_df.loc[IDs, ('dust_mass',0)]#/1e-5
            Q0_elements = input_df.loc[IDs, ('Q_0',0)]#/1e53
            a_star_elements = input_df.loc[IDs, ('a_star',0)]#/1e8

            f_esc.extend(f_esc_elements)
            halo_masses.extend(group_mass_elements)
            metal.extend(metal_elements)
            rel_star_mass.extend(star_mass_elements/group_mass_elements)
            rel_gas_mass.extend(gas_elements/group_mass_elements)
            rel_dust_mass.extend(dust_mass_elements/group_mass_elements)
            redshifts.extend(np.full(len(IDs),z[i]))
            all_IDs.extend(IDs)
            Q0.extend(Q0_elements)
            a_star.extend(a_star_elements)
            
        dataset = pd.DataFrame({'ID':all_IDs, 'z':redshifts, 
                            'HaloMass':halo_masses, 'Metalicity':metal, 
                            'FractionStars':rel_star_mass, 'FractionGas':rel_gas_mass,
                            'FractionDust':rel_dust_mass, 'Q0':Q0, 'aStar':a_star, 'f_esc':f_esc})
        store[setting] = dataset
    store.close()
    
    return 

if __name__ == "__main__":
    # Set some global variables such as the redshifts of the snapshots and the simulation name
    z = [6,8,10]
    simname = 'L35n2160TNG'
    h = 6.62607004e-34
    e = 1.60217662e-19
    j_to_erg = 1e7
    norm = 1e52

    # load the config module for the fiducial 2 dust configuration
    spec_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2d/config.py")
    config_dust = importlib.util.module_from_spec(spec_dust)
    spec_dust.loader.exec_module(config_dust)

    # load the config module for the fiducial 2 no dust configuration
    spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2/config.py")
    config_no_dust = importlib.util.module_from_spec(spec_no_dust)
    spec_no_dust.loader.exec_module(config_no_dust)

    halos = construct_halo_dict(simname, config_dust, config_no_dust)
    construct_dataframe(dictionary=halos, name = 'df_f_esc.h5')
