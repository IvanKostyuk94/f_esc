import pandas as pd
import numpy as np
import os
from crashpy.dataclasses.simulation import LoadedSimulation as Sim
from crashpy.utilities import crashMemMap
from get_spectrum import get_spectra
import pickle
from matplotlib.colors import Normalize
import illustris_python as il

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089, Ob0=0.0486, Tcmb0=2.725)

def lum_str_to_float(string):
    string = string.replace('d', 'e')
    return float(string)

def redshift_to_snap(redshift):
    snapnames = {6:'sn013', 8:'sn008', 10:'sn004'}
    return snapnames[redshift]

def scale_to_red(a):
    return 1/a-1

def star_ages(ID, redshift, conf, star_df):
    snap = redshift_to_snap(redshift)
    path_sources = f'/ptmp/mpa/mglatzle/TNG_f_esc/{conf}/run/L35n2160TNG/{snap}/g{ID}/Input/sources_ic00.in'
    
    star_ID_names = pd.read_csv(path_sources, delim_whitespace=True, header=None, usecols=[4])
    star_IDs = star_ID_names[4].apply(lambda element: int(element[1:-5]))

    a_formation = star_df.loc[star_IDs]['formation']
    z_formation = scale_to_red(a_formation)
    ages = (cosmo.age(redshift).value-cosmo.age(z_formation).value+0.005)*1000
    return ages

def add_age(df, conf, star_dic):
    ages = []
    empty_ages = []
    # for index, row in df.iterrows():
    #     empty_ages.append([])
    #     df['StellarAges'] = pd.Series(empty_ages)
    for index, row in df.iterrows():
        age = star_ages(ID=row.ID, redshift=row.z, conf=conf, star_df=star_dic[row.z])
        ages.append(age)
        #df.StellarAges[index] = age
        try:
            if len(age) != len(row.per_source):
                raise ValueError(f'The number of age elements and sources should be the same, got {len(age)} and {len(row.per_source)} instead')
        except:
            print(age)
            print(row.per_source)
    df['StellarAges'] = pd.Series(ages)
    return

def build_star_mass_df(basePath):
    stars_6 = il.snapshot.loadSubset(basePath, 13, 'stars', fields=['ParticleIDs', 'Masses', 'GFM_StellarFormationTime'])
    stars_8 = il.snapshot.loadSubset(basePath, 8, 'stars', fields=['ParticleIDs', 'Masses', 'GFM_StellarFormationTime'])
    stars_10 = il.snapshot.loadSubset(basePath, 4, 'stars', fields=['ParticleIDs', 'Masses', 'GFM_StellarFormationTime'])
    
    stars_6_df = pd.DataFrame({'mass':stars_6['Masses'],'formation':stars_6['GFM_StellarFormationTime']}, index=stars_6['ParticleIDs'])
    stars_8_df = pd.DataFrame({'mass':stars_8['Masses'], 'formation':stars_8['GFM_StellarFormationTime']}, index=stars_8['ParticleIDs'])
    stars_10_df = pd.DataFrame({'mass':stars_10['Masses'], 'formation':stars_10['GFM_StellarFormationTime']}, index=stars_10['ParticleIDs'])
    
    stars_6_df.sort_index(inplace=True)
    stars_8_df.sort_index(inplace=True)
    stars_10_df.sort_index(inplace=True)
    
    star_dic = {6:stars_6_df, 8:stars_8_df, 10:stars_10_df}
    return star_dic

# Calculate the star mass in the galaxy (0.2R_vir)
def star_mass_gal(ID, redshift, conf, star_df, side_length):
    snap = redshift_to_snap(redshift)
    path_sources = f'/ptmp/mpa/mglatzle/TNG_f_esc/{conf}/run/L35n2160TNG/{snap}/g{ID}/Input/sources_ic00.in'
    coord = pd.read_csv(path_sources, delim_whitespace=True, header=None, usecols=[0,1,2])+0.5-side_length/2
    rel_dist = np.sqrt(np.sum(coord**2, axis=1))*2/side_length

    star_ID_names = pd.read_csv(path_sources, delim_whitespace=True, header=None, usecols=[4])
    star_IDs = star_ID_names[4].apply(lambda element: int(element[1:-5]))

    star_masses = np.array(star_df.loc[star_IDs]['mass'])
    
    tot_star_mass = 0
    for i in range(len(rel_dist)):
        if rel_dist[i]<0.2:
            tot_star_mass += star_masses[i]/h
    return tot_star_mass*1e10/h

def redshift_to_snap(redshift):
    correspondense = {6:'sn013', 8:'sn008', 10:'sn004'}
    return correspondense[redshift]

def get_simulation_path(halo_id, conf, redshift):
    snap = redshift_to_snap(redshift)
    conf_dir = os.path.join('/ptmp/mpa/mglatzle/TNG_f_esc', conf)
    simulation_path =  os.path.join(conf_dir, f'run/L35n2160TNG/{snap}/g{halo_id}/Output/phys_ic00_rt05.out')
    density_path = os.path.join(conf_dir, f'run/L35n2160TNG/{snap}/g{halo_id}/Input/dens_ic00.in')
    return simulation_path, density_path

def gal_clumping_factor(density_map):
    volume = len(density_map)
    C = np.sum(np.square(density_map))*volume/(np.sum(density_map)**2)
    return C

def get_average_quantities(halo_id, conf, redshift, df_star_masses):
    path_sim, path_dens = get_simulation_path(halo_id, conf, redshift)
    dens = crashMemMap(path_dens, 'all')[0]
    r_gal = dens.shape[0]/2*0.2
    center = np.array(dens.shape)/2
    cut_low = int(np.floor(center[0]-r_gal))
    cut_up = int(np.ceil(center[0]+r_gal))
    galaxy_box = dens[cut_low:cut_up, cut_low:cut_up, cut_low:cut_up] 
    gal_size = galaxy_box.shape[0]
    gal_center = galaxy_box.shape[0]/2
    inside_gal = []
    for i in range(gal_size):
       for j in range(gal_size):
           for k in range(gal_size):
               if np.sqrt((i-gal_center)**2+(j-gal_center)**2+(k-gal_center)**2)<=r_gal:
                   inside_gal.append(galaxy_box[i,j,k])
    inside_gal = np.array(inside_gal)
    c_gal = gal_clumping_factor(inside_gal)
    sigma_gas_gal = inside_gal.sum()/(np.pi*r_gal**2)
    star_mass_in_gal = star_mass_gal(halo_id, redshift, conf, star_df=df_star_masses, side_length=dens.shape[0])
    sigma_star_gal = star_mass_in_gal/(np.pi*r_gal**2)
    return (c_gal, sigma_gas_gal, sigma_star_gal)

def update_df(df_name, name_star_mass_dic=None, new_df_name=None, conf='fid2'):
    path_to_data = '/u/ivkos/analysis/dfs/data'
    path_to_dfs = '/u/ivkos/analysis/dfs'
    
    if new_df_name == None:
        new_df_name = df_name.split('.')[0]+'_updated.pickle'
    path_to_df = os.path.join(path_to_dfs, df_name)
    path_to_new_df = os.path.join(path_to_dfs, new_df_name)
    
    if name_star_mass_dic == None:
        name_star_mass_dic = 'star_masses.pickle'
    path_star_masses = os.path.join(path_to_data, name_star_mass_dic)
    with open(path_star_masses, 'rb') as handle:
        star_masses = pickle.load(handle)

    df = pd.read_pickle(path_to_df)
    gal_clump = []
    surface_gas_gal = []
    surface_star_gal = []

    for i, halo in df.iterrows():
        halo_id = halo['ID']
        z = halo['z']
        print(f'Working on halo {halo_id} at z={z}')
        df_star_masses = star_masses[z]
        c_gal, sig_gas, sig_star = get_average_quantities(halo_id, conf, z, df_star_masses)
        gal_clump.append(c_gal)
        surface_gas_gal.append(sig_gas)
        surface_star_gal.append(sig_star)
    
    df['clump_gas'] = gal_clump
    df['sigma_gas_gal'] = surface_gas_gal
    df['sigma_star_gal'] = surface_star_gal
    add_age(df, conf=conf, star_dic=star_masses)
    
    df.to_pickle(path_to_new_df)
    return

if __name__ == "__main__":
    path_to_data = '/u/ivkos/analysis/dfs/data'
    basePath = '/virgo/simulations/IllustrisTNG/L35n2160TNG/output'
    basePath_2 = '/virgo/simulations/IllustrisTNG/L35n1080TNG/output'
    basePath_3 = '/virgo/simulations/IllustrisTNG/L35n540TNG/output'

    # mass_tng = build_star_mass_df(basePath)
    # path_to_dump = os.path.join(path_to_data, 'star_masses.pickle')
    # with open(path_to_dump, 'wb') as handle:
    #     pickle.dump(mass_tng, handle)

    # mass_tng2 = build_star_mass_df(basePath_2)
    # path_to_dump2 = os.path.join(path_to_data, 'star_masses_tng2.pickle')
    # with open(path_to_dump2, 'wb') as handle:
    #     pickle.dump(mass_tng2, handle)

    # mass_tng3 = build_star_mass_df(basePath_3)
    # path_to_dump3 = os.path.join(path_to_data, 'star_masses_tng3.pickle')
    # with open(path_to_dump3, 'wb') as handle:
    #     pickle.dump(mass_tng3, handle)


    h = 0.6774

    df_name = 'esc_analysis.pickle'
    name_star_mass_dic = 'star_masses.pickle'
    new_df_name = 'esc_analysis_updated.pickle'
    conf = 'esc_analysis'
    update_df(df_name, name_star_mass_dic=None, new_df_name=None, conf=conf)




