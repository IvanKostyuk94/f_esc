import importlib.util
import h5py
from crashpy.dataclasses.simulation import LoadedSimulation as Sim
import numpy as np
# load the config module for the fiducial 2 dust configuration
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


simname = 'L35n2160TNG'
spec = importlib.util.spec_from_file_location("module.name","/freya/ptmp/mpa/mglatzle/TNG_f_esc/fid2/config.py")
config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(config)
snaps = config.to_process[simname]


rundir = config.rundir(simname, config.snap_name(13), config.halo_name(2081))
try:
    csim = Sim(rundir)
    print(rundir)
except:
    print(f'An error occured in {rundir}')
    exit()
ls = []
for i, pf in enumerate(csim.getAllphysfiles()):
    with h5py.File(config.f_esc_file(pf),'r') as f:
        fesc = h5toDict(f)
        print(fesc)
        if not (fesc==None): 
            exit()
    ls.append(fesc)
    # except:
    #     exit()
    #     pass