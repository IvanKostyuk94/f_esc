df = pd.read_pickle('dfs/full_esc_updated.pickle')
snaps = ['sn013', 'sn008', 'sn004']
red_dict = {'sn013':6, 'sn008':8, 'sn004':10}
testing2 = {}
for snap in snaps:
    testing2[snap] = []
    snap_path = f'/ptmp/mpa/mglatzle/TNG_f_esc/full_esc/input/L35n2160TNG/{snap}'
    full_esc_path = f'/ptmp/mpa/mglatzle/TNG_f_esc/full_esc/input/L35n2160TNG/{snap}/sel_halos_df.pickle'
    fid2_path = f'/ptmp/mpa/mglatzle/TNG_f_esc/fid2/input/L35n2160TNG/{snap}/sel_halos_df.pickle'
    halo_df_full_esc = pd.read_pickle(full_esc_path)
    halo_df_fid2 = pd.read_pickle(fid2_path)
    counter=0
    for _, row in df[df.z==red_dict[snap]].iterrows():
#         if counter>100:
#             break
#         counter += 1
        try:
            ID = row['ID']
            cell_volume = np.prod(halo_df_full_esc.loc[ID]['size']/halo_df_full_esc.loc[ID]['shape'])
            halo_gas_mass = halo_df_full_esc.loc[ID]['GroupMassType',0]
            dens_path = snap_path+f'/g{ID}/dens_ic00.in'
            dens = crashMemMap(dens_path, 'all')[0]
            dens_sum = np.sum(dens)
            dens_conv_sum = tng_gas_dens(dens_sum, z=red_dict[snap])
            gas_mass_grid = dens_conv_sum*cell_volume
            halo_df_full_esc.loc[ID]['gas_mass',0] = gas_mass_grid
            halo_df_fid2.loc[ID]['gas_mass',0] = gas_mass_grid
            testing2[snap].append(gas_mass_grid/halo_gas_mass)
        except:
            print(ID)
            continue
    halo_df_full_esc.to_pickle(full_esc_path)
    halo_df_fid2.to_pickle(fid2_path)
        #conv_gas_mass = utils.phys_gas_number_dens(halo_gas_mass/cell_volume, red_dict[snap], TNGcosmo.h).value
    testing2[snap] = pd.Series(testing2[snap])
