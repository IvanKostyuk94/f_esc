import pandas as pd
import numpy as np
import os
import pickle


def redshift_to_snap(redshift):
    snapnames = {6:'sn013', 8:'sn008', 10:'sn004'}
    return snapnames[redshift]


def star_metal(ID, redshift, conf, star_df):
    snap = redshift_to_snap(redshift)
    path_sources = f'/ptmp/mpa/mglatzle/TNG_f_esc/{conf}/run/L35n2160TNG/{snap}/g{ID}/Input/sources_ic00.in'
    
    star_ID_names = pd.read_csv(path_sources, delim_whitespace=True, header=None, usecols=[4])
    star_IDs = star_ID_names[4].apply(lambda element: int(element[1:-5]))

    metal = star_df.loc[star_IDs]['metallicity']
    return metal


def add_metal(df, conf, star_dic):
    metals = []
    for index, row in df.iterrows():
        print(f'Working on halo {row.ID} at redshift {row.z}')
        metal = star_metal(ID=row.ID, redshift=row.z, conf=conf, star_df=star_dic[row.z])
        metals.append(metal)
        try:
            if len(metal) != len(row.per_source):
                raise ValueError(f'The number of metal elements and sources should be the same, got {len(metal)} and {len(row.per_source)} instead')
        except:
            print(metal)
            print(row.per_source)
    df['StellarMetallicities'] = pd.Series(metals)
    return


df = pd.read_pickle('dfs/esc_analysis_updated.pickle')

with open('extended_star_dict.pickle', 'rb') as f:
    star_dict = pickle.load(f)


add_metal(df, 'esc_analysis', star_dict)

df.to_pickle('dfs/esc_analysis_updated_metals.pickle')