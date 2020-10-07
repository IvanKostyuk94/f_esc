import os
from crashpy.crash2VTK.crash2VTK import CRASH_gridfile2RectilinearGridfile

def redshift_to_snap(redshift):
    correspondense = {6:'sn013', 8:'sn008', 10:'sn004'}
    return correspondense[redshift]

def get_simulation_path(halo_id, conf, redshift):
    snap = redshift_to_snap(redshift)
    conf_dir = os.path.join('/ptmp/mpa/mglatzle/TNG_f_esc', conf)
    simulation_path =  os.path.join(conf_dir, f'run//L35n2160TNG/{snap}/g{halo_id}/Output/phys_ic00_rt05.out')
    return simulation_path

def create_vtr(halo_id, conf, redshift):
    simulation_path = get_simulation_path(halo_id, conf, redshift)
    output_path = f'/u/ivkos/analysis/halos_vtr/{conf}_z{redshift}_{halo_id}.vtr'
    CRASH_gridfile2RectilinearGridfile(simulation_path, output_path)
    return

create_vtr(7085,'fid2',6)