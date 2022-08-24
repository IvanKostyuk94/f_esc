import numpy as np
import os
import pandas as pd

from crashpy.dataclasses import spectrum
import scipy.constants as constants


def get_nearest_lum_source(faint_source, lum_sources):
    distances = []
    for _, source in lum_sources.iterrows():
        distance = np.sqrt(
            (source.x - faint_source.x) ** 2
            + (source.y - faint_source.y) ** 2
            + (source.z - faint_source.z) ** 2
        )
        distances.append(distance)
    distances = np.array(distances)
    min_index = np.argmin(distances)
    min_value = distances[min_index]
    index = lum_sources.iloc[min_index].name
    return min_value, index


def source_merger(source_df, lum_source_idx, faint_source_idx, path_to_halo_spec):
    faint_source = source_df.loc[faint_source_idx]
    lum_source = source_df.loc[lum_source_idx]

    faint_spec_path = os.path.join(
        path_to_halo_spec, faint_source["name"][1:-1] + ".ps"
    )
    lum_spec_path = os.path.join(path_to_halo_spec, lum_source["name"][1:-1] + ".ps")

    spec_faint = spectrum.CrashSpectrum.fromFile(faint_spec_path)
    spec_lum = spectrum.CrashSpectrum.fromFile(lum_spec_path)

    lum_faint = faint_source["Q"]
    lum_lum = lum_source["Q"]

    spec_faint.setEmissivity(lum_faint)
    spec_lum.setEmissivity(lum_lum)

    summed_spectrum = spectrum.CrashSpectrum(
        constants.e / constants.h * spec_lum.spectrum["hnu"],
        spec_lum.spectrum["Lnu"] + spec_faint.spectrum["Lnu"],
    )

    out_path = lum_spec_path
    summed_spectrum.toFile(path=out_path, f=True)

    source_df.loc[lum_source_idx, "Q"] = lum_faint + lum_lum
    source_df.drop(faint_source_idx, inplace=True)
    return


def merge_halo_sources(
    path_to_input, path_to_spectra, snap, halo_id, threshold=0.01, max_dist=1.9
):
    column_names = [
        "x",
        "y",
        "z",
        "a",
        "Q",
        "b",
        "name",
        "c",
        "halo",
        "d",
        "packets",
        "e",
        "loc_uv",
        "f",
        "loc_x",
    ]
    sources_path = os.path.join(path_to_input, snap, f"g{halo_id}", "sources_ic00.in")
    spectra_path = os.path.join(path_to_spectra, snap, f"g{halo_id}")

    source_df = pd.read_table(sources_path, delimiter="\t", names=column_names)
    source_df["Q"] = source_df["Q"].str.replace("d", "e").astype("float")

    lum_sources = source_df[source_df["Q"] > threshold * source_df["Q"].max()]
    faint_sources = source_df[source_df["Q"] < threshold * source_df["Q"].max()]
    for faint_idx, source in faint_sources.iterrows():
        dist, lum_source_idx = get_nearest_lum_source(
            faint_source=source, lum_sources=lum_sources
        )
        if dist > max_dist:
            source_df.drop(faint_idx, inplace=True)
        else:
            source_merger(
                source_df=source_df,
                lum_source_idx=lum_source_idx,
                faint_source_idx=faint_idx,
                path_to_halo_spec=spectra_path,
            )

    source_df["Q"] = source_df["Q"].apply("{:.3e}".format)
    source_df["Q"] = source_df["Q"].str.replace("e", "d")
    source_df.to_csv(sources_path, sep="\t", header=False, index=False)
    return


def get_all_halos(path, snap):
    path_to_snap = os.path.join(path, snap)
    files = os.listdir(path_to_snap)
    halos = [int(halo[1:]) for halo in files if halo.startswith("g")]
    return halos


def all_halo_merge(path_to_input, path_to_spectra, threshold=0.01, max_dist=5):
    snaps = ["sn004", "sn008", "sn013"]
    for snap in snaps:
        halos = get_all_halos(path_to_input, snap)
        for halo in halos:
            print(f"Merging sources in halo {halo} in snap {snap}")
            column_names = [
                "x",
                "y",
                "z",
                "a",
                "Q",
                "b",
                "name",
                "c",
                "halo",
                "d",
                "packets",
                "e",
                "loc_uv",
                "f",
                "loc_x",
            ]
            sources_path = os.path.join(
                path_to_input, snap, f"g{halo}", "sources_ic00.in"
            )
            spectra_path = os.path.join(path_to_spectra, snap, f"g{halo}")
            try:
                source_df = pd.read_table(
                    sources_path, delimiter="\t", names=column_names
                )
            except:
                print(
                    f"Halo {halo} in snap {snap} does not seem to have a sources file!"
                )
            source_df["Q"] = source_df["Q"].str.replace("d", "e").astype("float")

            lum_sources = source_df[source_df["Q"] > threshold * source_df["Q"].max()]
            faint_sources = source_df[source_df["Q"] < threshold * source_df["Q"].max()]
            for faint_idx, source in faint_sources.iterrows():
                dist, lum_source_idx = get_nearest_lum_source(
                    faint_source=source, lum_sources=lum_sources
                )
                if dist > max_dist:
                    source_df.drop(faint_idx, inplace=True)
                else:
                    source_merger(
                        source_df=source_df,
                        lum_source_idx=lum_source_idx,
                        faint_source_idx=faint_idx,
                        path_to_halo_spec=spectra_path,
                    )

            source_df["Q"] = source_df["Q"].apply("{:.3e}".format)
            source_df["Q"] = source_df["Q"].str.replace("e", "d")
            source_df.to_csv(sources_path, sep="\t", header=False, index=False)
    return


def set_packet_num(input_dir, packet_num):
    column_names = [
        "x",
        "y",
        "z",
        "a",
        "Q",
        "b",
        "name",
        "c",
        "halo",
        "d",
        "packets",
        "e",
        "loc_uv",
        "f",
        "loc_x",
    ]
    for snap in os.listdir(input_dir):
        packets = packet_num
        snap_dir = os.path.join(input_dir, snap)
        for halo in os.listdir(snap_dir):
            try:
                print(
                    f"Setting the packet number to {packets} in halo {halo} of snap {snap}"
                )
                halo_dir = os.path.join(snap_dir, halo)
                sources_path = os.path.join(halo_dir, "sources_ic00.in")

                sources = pd.read_table(
                    sources_path, delimiter="\t", names=column_names
                )
                sources.packets = packets
                sources.to_csv(sources_path, sep="\t", header=False, index=False)
            except:
                print(halo)
                continue


def set_loc_esc(input_dir, esc_frac):
    column_names = [
        "x",
        "y",
        "z",
        "a",
        "Q",
        "b",
        "name",
        "c",
        "halo",
        "d",
        "packets",
        "e",
        "loc_uv",
        "f",
        "loc_x",
    ]
    for snap in os.listdir(input_dir):
        snap_dir = os.path.join(input_dir, snap)
        for halo in os.listdir(snap_dir):
            try:
                print(
                    f"Setting the local escape fraction to {esc_frac} in halo {halo} of snap {snap}"
                )
                halo_dir = os.path.join(snap_dir, halo)
                sources_path = os.path.join(halo_dir, "sources_ic00.in")

                sources = pd.read_table(
                    sources_path, delimiter="\t", names=column_names
                )
                sources.loc_uv = esc_frac
                sources.loc_x = esc_frac
                sources.to_csv(sources_path, sep="\t", header=False, index=False)
            except:
                print(halo)
                continue


def remove_heavy_halos(path, star_mass_lim=1e8):
    snaps = ["sn004", "sn008", "sn013"]
    df_name = "sel_halos_df.pickle"
    for snap in snaps:
        halos = get_all_halos(path, snap)
        df_path = os.path.join(path, snap, df_name)
        all_halos = pd.read_pickle(df_path)
        halos_used = all_halos.loc[halos]

        star_masses = halos_used[("GroupMassType", 4)] * 1e10 / 0.6774
        to_remove = np.array(star_masses[star_masses > star_mass_lim].index)
        halos_to_remove = [f"g{halo}" for halo in to_remove]

        for halo in halos_to_remove:
            path_to_halo = os.path.join(path, snap, halo)
            print(f"Removing {path_to_halo}")
            os.system(f"rm -r {path_to_halo}")
    return


def prepare_halos(conf, threshold=0.01, max_dist=5, packet_number="1.000d+08"):
    basepath = "/ptmp/mpa/mglatzle/TNG_f_esc"
    path_to_conf = os.path.join(basepath, conf)
    path_to_input = os.path.join(path_to_conf, "input", "L35n2160TNG")
    path_to_spectra = os.path.join(path_to_conf, "db", "SPECTRA", "L35n2160TNG")

    remove_heavy_halos(path_to_input, star_mass_lim=1e9)
    all_halo_merge(path_to_input, path_to_spectra, threshold=0.01, max_dist=5)
    set_packet_num(path_to_input, packet_number)
    return


prepare_halos("new_1e-1")
prepare_halos("new_5e-1")
prepare_halos("new_7e-1")
prepare_halos("new_dust")
