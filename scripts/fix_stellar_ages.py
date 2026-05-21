import pandas as pd
import os
import pickle
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089, Ob0=0.0486, Tcmb0=2.725)


def redshift_to_snap(redshift):
    snapnames = {6: "sn013", 8: "sn008", 10: "sn004"}
    return snapnames[redshift]


def scale_to_red(a):
    return 1 / a - 1


def star_ages(ID, redshift, conf, star_df):
    snap = redshift_to_snap(redshift)
    path_sources = f"/ptmp/mpa/mglatzle/TNG_f_esc/{conf}/run/L35n2160TNG/{snap}/g{ID}/Input/sources_ic00.in"

    star_ID_names = pd.read_csv(
        path_sources, delim_whitespace=True, header=None, usecols=[4]
    )
    star_IDs = star_ID_names[4].apply(lambda element: int(element[1:-5]))

    a_formation = star_df.loc[star_IDs]["formation"]
    z_formation = scale_to_red(a_formation)
    ages = (cosmo.age(redshift).value - cosmo.age(z_formation).value + 0.005) * 1000
    return ages


def add_age(df, conf, star_dic):
    for index, row in df.iterrows():
        print(f"Working on Halo {row.ID} at redhshift {row.z}")

        try:
            if len(row.StellarAges) != len(row.per_source):
                age = star_ages(
                    ID=row.ID, redshift=row.z, conf=conf, star_df=star_dic[row.z]
                )
                df.StellarAges[index] = age
                print(
                    f"Correction at {index} per_source: {len(row.per_source)}, stellar ages:{len(row.StellarAges)}"
                )
        except:
            age = star_ages(
                ID=row.ID, redshift=row.z, conf=conf, star_df=star_dic[row.z]
            )
            df.StellarAges[index] = age
            print(
                f"Correction at {index} per_source: {len(row.per_source)}, stellar ages:{len(row.StellarAges)}"
            )

    # df['StellarAges'] = pd.Series(ages)
    return


df_path = "dfs/full_esc_updated.pickle"
df = pd.read_pickle(df_path)
conf = "full_esc"
path_to_data = "/u/ivkos/analysis/dfs/data"
name_star_mass_dic = "star_masses.pickle"
path_star_masses = os.path.join(path_to_data, name_star_mass_dic)
with open(path_star_masses, "rb") as handle:
    star_masses = pickle.load(handle)

add_age(df=df, conf=conf, star_dic=star_masses)
df.to_pickle(df_path)
