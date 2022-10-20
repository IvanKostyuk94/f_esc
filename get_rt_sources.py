import pandas as pd
import os

def z_to_snap(z):
    z_snap = {6:'sn013', 8:'sn008', 10:'sn004'}
    return z_snap[z]

def get_sources_path(id, z):
    base = '/ptmp/mpa/mglatzle/TNG_f_esc/new_main/input/L35n2160TNG'
    snap = z_to_snap(z)
    halo = 'g'+str(id)
    source_file = 'sources_ic00.in'
    source_path = os.path.join(base, snap, halo, source_file)
    return source_path

def get_num_sources_from_path(source_path):
    with open(source_path, 'r') as f:
        num = len(f.readlines())
    return num

def get_num_sources(id,z):
    source_path = get_sources_path(id, z)
    num = get_num_sources_from_path(source_path)
    return num

def df_num_sources(df):
    num_sources = []
    for _, row in df.iterrows():
        num = get_num_sources(row['ID'], row['z'])
        num_sources.append(num)
    df['Num_sources_rt'] = num_sources