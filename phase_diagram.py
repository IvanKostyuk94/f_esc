import pandas as pd
import numpy as np
from pyTNG.cosmology import TNGcosmo
from pyTNG import data_interface as _data_interface
import os
import illustris_python as il
import astropy.units as u
from astropy.constants import m_p, k_B
from matplotlib import pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import binned_statistic_2d
from crashpy.utilities import crashMemMap

h = TNGcosmo.h


def get_sim():
    basepath = "/virgotng/universe/IllustrisTNG/"
    sim_name = "L35n2160TNG"
    sim = _data_interface.TNG50Simulation(os.path.join(basepath, sim_name))
    sim_path = os.path.join(basepath, sim_name, "output")
    return sim, sim_path


def get_snap_number(z, numerical=True):
    if numerical:
        z_to_snap = {6: 13, 8: 8, 10: 4}
    else:
        z_to_snap = {6: "sn013", 8: "sn008", 10: "sn004"}
    return z_to_snap[z]


def get_paths(id, z, conf="new_main"):
    snap = get_snap_number(z, numerical=False)
    conf_dir = os.path.join("/ptmp/mpa/mglatzle/TNG_f_esc", conf)
    simulation_path = os.path.join(
        conf_dir, f"run/L35n2160TNG/{snap}/g{id}/Output/phys_ic00_rt05.out"
    )
    input_path = os.path.join(
        conf_dir, f"run/L35n2160TNG/{snap}/g{id}/Input/temp_ic00.in"
    )
    dens_path = os.path.join(
        conf_dir, f"run/L35n2160TNG/{snap}/g{id}/Input/dens_ic00.in"
    )
    return input_path, simulation_path, dens_path


def get_maps(id, z, get_output, conf="esc_analysis"):
    input_path, sim_path, dens_path = get_paths(id, z, conf=conf)
    if get_output:
        temp = crashMemMap(sim_path, 1)
    else:
        temp = crashMemMap(input_path, 1)
    dens = crashMemMap(dens_path, "all")

    maps = {}
    maps["temp"] = temp
    maps["dens"] = dens[0]

    return maps


def get_cell_volume(id, z, maps, conf="esc_analysis"):
    df_dir = "/u/ivkos/analysis/dfs"
    df_path = os.path.join(df_dir, conf + ".pickle")
    df = pd.read_pickle(df_path)
    to_cm = (1 * u.kpc).to(u.cm).value / h
    halo_rad = df[(df.ID == id) & (df.z == z)].HaloRadii.values[0] * to_cm
    cells = maps[next(iter(maps))].shape[0]
    cell_size = 2 * halo_rad / cells
    cell_volume = cell_size**3
    return cell_volume


def get_temperature(x_e, internal_en):
    X_H = 0.76
    mu = 4 * m_p / (1 + 3 * X_H + 4 * X_H * x_e)
    U = internal_en * u.km**2 * u.s ** (-2)
    gamma = 5 / 3
    temperature = (gamma - 1) * U * mu / k_B
    return temperature.to(u.K).value


def get_number_dens(density, z):

    mass_to_kg = (1 * u.Msun).to(u.kg).value * 1e10 / h
    volume_to_cm3 = ((1 * u.kpc).to(u.cm).value / h / (1 + z)) ** 3

    X_H = 0.76
    mu = m_p.value * (4 - 3 * X_H)

    num_dens = density * mass_to_kg / volume_to_cm3 / mu
    return num_dens


def get_gas(z, id):
    _, sim_path = get_sim()
    gas = il.snapshot.loadHalo(sim_path, get_snap_number(z), id, "gas")

    mass_to_sun = 1e10 / h

    gas_dict = {}
    gas_dict["density"] = get_number_dens(gas["Density"], z)
    gas_dict["mass"] = gas["Masses"] * mass_to_sun
    gas_dict["temperature"] = get_temperature(
        gas["ElectronAbundance"], gas["InternalEnergy"]
    )
    return gas_dict


def get_mapped_gas(id, z, get_output, conf="esc_analysis"):
    maps = get_maps(id, z, get_output, conf=conf)
    cell_vol = get_cell_volume(id, z, maps, conf)

    gas = {}
    gas["density"] = np.log10(maps["dens"].flatten())
    gas["temperature"] = np.log10(maps["temp"].flatten())
    gas["mass"] = maps["dens"].flatten() * cell_vol
    return gas


def gas_to_log_scale(gas):
    gas["density"] = np.log10(gas["density"])
    gas["temperature"] = np.log10(gas["temperature"])
    return gas


def get_bin_edges(gas, nx, ny):
    x_bins = np.linspace(gas["density"].min(), gas["density"].max(), num=nx)
    y_bins = np.linspace(gas["temperature"].min(), gas["temperature"].max(), num=ny)
    return x_bins, y_bins


def plot_parameters(params, multiple=False):
    parameters = {}
    parameters["x_labelsize"] = 50
    parameters["y_labelsize"] = 50

    parameters["length_major_ticks"] = 16
    parameters["length_minor_ticks"] = 8
    parameters["width_minor_ticks"] = 3
    parameters["width_major_ticks"] = 4
    parameters["labelsize_x_ticks"] = 35
    parameters["labelsize_y_ticks"] = 35

    parameters["colorbar_labelsize"] = 50
    parameters["colorbar_ticklabelsize"] = 35

    parameters["axes_width"] = 3

    parameters["figure_width"] = 40
    parameters["figure_height"] = 35

    parameters["x_label"] = r"$\log(n) [\mathrm{cm}^{-3}]$"
    parameters["y_label"] = r"$\log(T) [\mathrm{K}]$"
    parameters["bar_label"] = r"$\log(\frac{M}{M_\mathrm{max}})$"

    parameters["nx"] = 45
    parameters["ny"] = 30

    parameters["v_min"] = -4
    parameters["v_max"] = 0

    parameters["x_lim_min"] = -4.8
    parameters["x_lim_max"] = 0

    parameters["y_lim_min"] = 2.8
    parameters["y_lim_max"] = 6.4

    if params != None:
        for element in params:
            parameters[element] = params[element]
    return parameters


def get_hist(gas, parameters):
    x_bins, y_bins = get_bin_edges(gas, parameters["nx"], parameters["ny"])
    hist, *_ = binned_statistic_2d(
        x=gas["density"],
        y=gas["temperature"],
        values=gas["mass"],
        bins=[x_bins, y_bins],
        statistic="sum",
    )
    hist = hist / hist.max()
    return hist, x_bins, y_bins


def get_col_norm(parameters):
    # v_min = np.log10(hist[hist > 0].min())
    # v_max = np.log10(hist.max())
    v_min = parameters["v_min"]
    v_max = parameters["v_max"]
    v_center = (v_max + v_min) / 2

    col_norm = colors.DivergingNorm(vmin=v_min, vcenter=v_center, vmax=v_max)
    return col_norm


def set_ax_params(ax, parameters):
    ax.set_xlabel(parameters["x_label"], size=parameters["x_labelsize"])
    ax.set_ylabel(parameters["y_label"], size=parameters["y_labelsize"])

    ax.set_xlim(parameters["x_lim_min"], parameters["x_lim_max"])
    ax.set_ylim(parameters["y_lim_min"], parameters["y_lim_max"])

    ax.tick_params(
        length=parameters["length_major_ticks"], width=parameters["width_major_ticks"]
    )
    ax.tick_params(
        length=parameters["length_minor_ticks"],
        width=parameters["width_minor_ticks"],
        which="minor",
    )
    return


def create_color_bar(ax, parameters, subfig, f, multiple=False):
    if multiple:
        f.subplots_adjust(right=0.9)
        cax = f.add_axes([0.95, 0.12, 0.05, 0.76])
    else:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = f.colorbar(subfig, cax=cax)
    cbar.set_label(parameters["bar_label"], size=parameters["colorbar_labelsize"])
    cbar.ax.tick_params(labelsize=parameters["colorbar_ticklabelsize"])
    return


def set_plt_params(parameters):
    plt.rc("axes", linewidth=parameters["axes_width"])
    plt.rc("xtick", labelsize=parameters["labelsize_x_ticks"])
    plt.rc("ytick", labelsize=parameters["labelsize_y_ticks"])

    plt.rcParams["figure.figsize"] = (
        parameters["figure_width"],
        parameters["figure_height"],
    )
    plt.tight_layout(rect=(0, 0, 1, 0.7))
    return


def plot_histogram(gas, params=None):

    parameters = plot_parameters(params)
    set_plt_params(parameters)

    hist, x_bins, y_bins = get_hist(gas, parameters)

    f, ax = plt.subplots()
    col_norm = get_col_norm(parameters)
    x_grid, y_grid = np.meshgrid(x_bins, y_bins)

    subfig = ax.pcolormesh(
        x_grid, y_grid, np.log10(hist.T), norm=col_norm, cmap=plt.get_cmap("inferno")
    )

    set_ax_params(ax, parameters)
    create_color_bar(ax, parameters, subfig, f)
    return


def get_halo_gases(zs, ids, conf="esc_analysis"):
    halo_gases = []
    for id, z in zip(ids, zs):
        gases = []
        gases.append(gas_to_log_scale(get_gas(z, id)))
        gases.append(get_mapped_gas(id, z, get_output=False, conf=conf))
        gases.append(get_mapped_gas(id, z, get_output=True, conf=conf))
        halo_gases.append(gases)
    return halo_gases


def plot_multiple_histograms(halo_gases, params=None):

    parameters = plot_parameters(params, multiple=True)
    set_plt_params(parameters)

    col_norm = get_col_norm(parameters)

    f, axs = plt.subplots(
        3, 2, sharex="col", sharey="row", gridspec_kw={"hspace": 0, "wspace": 0}
    )
    for i, gases in enumerate(halo_gases):
        for j, gas in enumerate(gases):
            hist, x_bins, y_bins = get_hist(gas, parameters)
            x_grid, y_grid = np.meshgrid(x_bins, y_bins)
            subfig = axs[j, i].pcolormesh(
                x_grid,
                y_grid,
                np.log10(hist.T),
                norm=col_norm,
                cmap=plt.get_cmap("inferno"),
            )
            if i == 1:
                parameters["y_label"] = ""
            set_ax_params(axs[j, i], parameters)
    create_color_bar(axs, parameters, subfig, f, multiple=True)
    return


def generate_histogram_plot(
    z, id, from_tng=True, get_output=False, plot_params=None, conf="esc_analysis"
):
    if from_tng:
        gas = get_gas(z, id)
        gas = gas_to_log_scale(gas)
    else:
        gas = get_mapped_gas(id, z, get_output, conf)

    plot_histogram(gas, params=plot_params)
    return


def generate_multiple_histograms(zs, ids, plot_params=None, conf="esc_analysis"):
    halo_gases = get_halo_gases(zs, ids, conf)
    plot_multiple_histograms(halo_gases, params=plot_params)
    return


if __name__ == "__main__":
    example_halos_1e7 = [(4779, 6, 0.84), (2606, 8, 0.556570), (519, 10, 0.643313)]
    example_halos_1e8 = [(221, 6, 0.489819), (8, 10, 0.565650)]
