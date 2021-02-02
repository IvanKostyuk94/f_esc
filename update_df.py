import pandas as pd
import numpy as np
import os
from crashpy.dataclasses.simulation import LoadedSimulation as Sim
from crashpy.utilities import crashMemMap
from get_spectrum import get_spectra
import pickle
from matplotlib.colors import Normalize
import illustris_python as il

basePath = '/virgo/simulations/IllustrisTNG/L35n2160TNG/output'
h = 0.6774
def build_star_mass_df():
    stars_6 = il.snapshot.loadSubset(basePath, 13, 'stars', fields=['ParticleIDs', 'Masses'])
    stars_8 = il.snapshot.loadSubset(basePath, 8, 'stars', fields=['ParticleIDs', 'Masses'])
    stars_10 = il.snapshot.loadSubset(basePath, 4, 'stars', fields=['ParticleIDs', 'Masses'])
    
    stars_6_df = pd.DataFrame({'mass':stars_6['Masses']}, index=stars_6['ParticleIDs'])
    stars_8_df = pd.DataFrame({'mass':stars_8['Masses']}, index=stars_8['ParticleIDs'])
    stars_10_df = pd.DataFrame({'mass':stars_10['Masses']}, index=stars_10['ParticleIDs'])
    
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

with open('star_masses.pickle', 'rb') as handle:
    star_masses = pickle.load(handle)

def update_df(df_name, name_star_mass_dic, new_df_name, conf='fid2'):
    with open('star_masses.pickle', 'rb') as handle:
        star_masses = pickle.load(handle)

    df = pd.read_pickle(df_name)
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
    
    df.to_pickle(new_df_name)
    return

df_name = 'df_no_dust_ages.pickle'
name_star_mass_dic = 'star_masses.pickle'
new_df_name = 'df_no_dust_updated.pickle'
update_df(df_name, name_star_mass_dic, new_df_name)




