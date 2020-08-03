import numpy as np 


def get_spectra(ID, snap, setting):
    if setting == 'dust':
        fid = 'fid2d'
    elif setting == 'no_dust':
        fid = 'fid2'
    else:
        raise ValueError(f'The setting must be "dust" or "no_dust" got {setting} instead')

    path = f'/freya/ptmp/mpa/mglatzle/TNG_f_esc/{fid}/run/L35n2160TNG/{snap}/g{ID}/Input/sources_ic00.in'
    path_spectra = f'/freya/ptmp/mpa/mglatzle/TNG_f_esc/{fid}/db/SPECTRA/L35n2160TNG/{snap}/g{ID}/'
    
    with open(path) as sources:
        matrix = [line.split() for line in sources]

    no_of_photons = np.zeros((len(matrix)))
    for i in range(len(matrix)):
        name = matrix[i][4].replace('\'','')+'.ps'
        spec = np.loadtxt(path_spectra+name, skiprows=1)

        no_of_photons[i] = 1e52*float(matrix[i][3].replace('d','e'))

    return no_of_photons