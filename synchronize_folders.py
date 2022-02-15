import pandas as pd
import numpy as np
import os
import fileinput
import sys


# Function which copies finished halos from the run folder of 'origin' 
# to the folder for f_esc computing 'destination'
def transfer_to_all_sources(origin, destination, data_dir='/ptmp/mpa/mglatzle/TNG_f_esc'):
    run_dir_in = os.path.join(data_dir, origin, 'run', 'L35n2160TNG')
    run_dir_out = os.path.join(data_dir, destination, 'run', 'L35n2160TNG')
    snap_dirs = ['sn013', 'sn008', 'sn004']
    for snap in snap_dirs:
        snap_dir_in = os.path.join(run_dir_in, snap)
        snap_dir_out = os.path.join(run_dir_out, snap)

        for halo in os.listdir(snap_dir_in):
            halo_dir_in = os.path.join(snap_dir_in, halo)
            halo_dir_out = os.path.join(snap_dir_out, halo)
            
            output_dir_in = os.path.join(halo_dir_in, 'Output')

            last_rt_in = os.path.join(output_dir_in, 'phys_ic00_rt05.out')
            if os.path.isfile(last_rt_in):
                print(f'working on halo {halo} in snap {snap}')
                if not os.path.isdir(halo_dir_out):
                    os.system(f'cp -r {halo_dir_in} {snap_dir_out}')
                else:
                    print(f'{halo} already exists')

# Corrects the path in the 00_fesc... file from the folder of the CRASH run to the folder of the f_esc calculation
def correct_fesc_file(origin, destination, data_dir='/ptmp/mpa/mglatzle/TNG_f_esc', from_freya=True):
    fesc_file = os.path.join(data_dir, destination, '00_f_esc.txt')
    for line in fileinput.input(fesc_file, inplace=1):
        line = line.replace(origin, destination)
        sys.stdout.write(line)
    fileinput.close()
    return

# Replaces all the links and files from the CRASH runs to the paths for the f_esc 
# calculations containing all the sources
def correct_paths(origin, destination, data_dir='/ptmp/mpa/mglatzle/TNG_f_esc', from_freya=True, with_dust=False):
    raven_path = '/raven/ptmp/ivkos'
    freya_path = '/ptmp/mpa/mglatzle/TNG_f_esc'
    run_dir_out = os.path.join(data_dir, destination, 'run', 'L35n2160TNG')
    run_dir_out_input = os.path.join(data_dir, destination, 'input', 'L35n2160TNG')
    
    init_path = os.path.join(data_dir, destination, 'bin', 'CRASH_5.4.3v09.in')
    modules_path = os.path.join(data_dir, destination, 'config', 'Modules')
    random_path = os.path.join(data_dir, destination, 'config', 'SYS_RANDOM.in')
    snap_dirs = ['sn013', 'sn008', 'sn004']
    for snap in snap_dirs:
        snap_dir_out = os.path.join(run_dir_out, snap)
        
        snap_dir_out_input = os.path.join(run_dir_out_input, snap)
        for halo in os.listdir(snap_dir_out):
            print(f'Working on {halo} in snap {snap}')
            halo_dir_out = os.path.join(snap_dir_out, halo)
            
            
            halo_dir_out_input = os.path.join(snap_dir_out_input, halo)

            config_dir_out = os.path.join(halo_dir_out, 'config')
            
            simul_file_main = os.path.join(config_dir_out, 'SIMULATION.in')

            input_dir = os.path.join(halo_dir_out, 'Input')
            old_init = os.path.join(halo_dir_out, 'SYS_INIT.in')
            old_modules = os.path.join(halo_dir_out, 'config', 'Modules')
            old_random = os.path.join(halo_dir_out, 'config', 'SYS_RANDOM.in')
            
            source_file = os.path.join(input_dir, 'sources_ic00.in')
            simul_file = os.path.join(input_dir, 'simul_ic00.in')
            dens_file = os.path.join(input_dir, 'dens_ic00.in')
            ion_file = os.path.join(input_dir, 'ion.in')
            temp_file = os.path.join(input_dir, 'temp_ic00.in')
            if with_dust:
                dust_dir = os.path.join(input_dir, 'dust')
            
            source_file_sources = os.path.join(halo_dir_out_input, 'sources_ic00.in')
            dens_file_sources = os.path.join(halo_dir_out_input, 'dens_ic00.in')
            simul_file_sources = os.path.join(halo_dir_out_input, 'simul_ic00.in')
            ion_file_sources = os.path.join(halo_dir_out_input, 'ion.in')
            temp_file_sources = os.path.join(halo_dir_out_input, 'temp_ic00.in')
            if with_dust:
                dust_dir_sources = os.path.join(halo_dir_out_input, 'dust')
            

            if (destination not in os.path.realpath(source_file)) or (not from_freya):
                for line in fileinput.input(simul_file_main, inplace=1):
                    if from_freya:
                        if (line[:8] == "'IN_DIR=") or (line[:9] == "'OUT_DIR="):
                            line = line.replace(origin, destination)
                    else:
                        if (line[:8] == "'IN_DIR=") or (line[:9] == "'OUT_DIR="):
                            line = line.replace(raven_path, freya_path)
                    sys.stdout.write(line)
            
                os.system(f'rm {source_file}')
                os.system(f'rm {simul_file}')
                os.system(f'rm {dens_file}')
                os.system(f'rm {ion_file}')
                os.system(f'rm {temp_file}')
                os.system(f'rm {old_init}')
                os.system(f'rm {old_modules}')
                os.system(f'rm {old_random}')
                if with_dust:
                    os.system(f'rm {dust_dir}')
                
                os.system(f'ln -s {source_file_sources} {input_dir}')
                os.system(f'ln -s {dens_file_sources} {input_dir}')
                os.system(f'ln -s {simul_file_sources} {input_dir}')
                os.system(f'ln -s {ion_file_sources} {input_dir}')
                os.system(f'ln -s {temp_file_sources} {input_dir}')
                os.system(f'ln -s {init_path} {old_init}')
                os.system(f'ln -s {modules_path} {old_modules}')
                os.system(f'ln -s {random_path} {old_random}')
                if with_dust:
                    os.system(f'ln -s {dust_dir_sources} {input_dir}')

            else:
                print((destination not in os.path.realpath(source_file)) and (not from_freya))
                print(f'{halo_dir_out} is already finished')
    fileinput.close()

    
origin = 'potato'
destination = 'new_1e-1'


#correct_fesc_file(origin=origin, destination=destination)
#transfer_to_all_sources(origin=origin, destination=destination)
#correct_paths(origin=origin, destination='new_3e-1', from_freya=False)
#correct_paths(origin=origin, destination='new_5e-1', from_freya=False)
#correct_paths(origin=origin, destination='new_7e-1', from_freya=False)
correct_paths(origin=origin, destination='new_dust', from_freya=False, with_dust=True)
