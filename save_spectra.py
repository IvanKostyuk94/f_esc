# Function which builds a dictionary containing the information about the frequency dependent escape fractions
# and the escape fractions of individual sources

from build_df import *

def build_freq_dataframe(dictionary, settings=['dust', 'no_dust'], name='freq_f_esc'):
    redshifts = {'sn013':6,'sn008':8,'sn004':10}
    return


# load the config module for the fiducial 2 dust configuration
spec_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2d/config.py")
config_dust = importlib.util.module_from_spec(spec_dust)
spec_dust.loader.exec_module(config_dust)

# load the config module for the fiducial 2 no dust configuration
spec_no_dust = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2/config.py")
config_no_dust = importlib.util.module_from_spec(spec_no_dust)
spec_no_dust.loader.exec_module(config_no_dust)

halos = construct_halo_dict(simname, config_dust, config_no_dust)