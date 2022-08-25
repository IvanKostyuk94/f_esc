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

h = TNGcosmo.h


def add_StarMass(df):
    df["StarMass"] = df["HaloMass"] * df["FractionStars"] * 1e10 / h
    return


def get_sim():
    basepath = "/virgotng/universe/IllustrisTNG/"
    sim_name = "L35n2160TNG"
    sim = _data_interface.TNG50Simulation(os.path.join(basepath, sim_name))
    sim_path = os.path.join(basepath, sim_name, "output")
    return sim, sim_path


def get_snap_number(z):
    z_to_snap = {6: 13, 8: 8, 10: 4}
    return z_to_snap[z]


def get_temperature(x_e, internal_en):
    X_H = 0.76
    mu = 4 * m_p / (1 + 3 * X_H + 4 * X_H * x_e)
    U = internal_en * u.km**2 * u.s ** (-2)
    gamma = 5 / 3
    temperature = (gamma - 1) * U * mu / k_B
    return temperature.to(u.K).value


def get_gas(z, id):
    _, sim_path = get_sim()
    gas = il.snapshot.loadHalo(sim_path, get_snap_number(z), id, "gas")

    mass_to_g = (1 * u.Msun).to(u.g).value * 1e10 / h
    volume_to_cm3 = ((1 * u.kpc).to(u.cm).value / h / (1 + z)) ** 3

    gas_dict = {}
    gas_dict["density"] = gas["Density"] * mass_to_g / volume_to_cm3
    gas_dict["mass"] = gas["Masses"] * mass_to_g
    gas_dict["temperature"] = get_temperature(
        gas["ElectronAbundance"], gas["InternalEnergy"]
    )
    return gas_dict


def get_bin_edges(gas, nx=30, ny=30):
    x_bins = np.logspace(
        np.log10(gas["density"].min()), np.log10(gas["density"].max()), num=nx
    )
    y_bins = np.logspace(
        np.log10(gas["temperature"].min()), np.log10(gas["temperature"].max()), num=ny
    )
    return x_bins, y_bins


def plot_parameters(params):
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

    parameters["figure_width"] = 25
    parameters["figure_height"] = 15

    parameters["x_label"] = r"$\rho [\frac{\mathrm{g}}{\mathrm{cm}^3}]$"
    parameters["y_label"] = r"$T [\mathrm{K}]$"
    parameters["bar_label"] = r"$M[\mathrm{g}]$"

    if params != None:
        for element in params:
            parameters[element.key] = element.value
    return parameters


def get_hist(gas):
    x_bins, y_bins = get_bin_edges(gas)
    hist, *_ = binned_statistic_2d(
        x=gas["density"],
        y=gas["temperature"],
        values=gas["mass"],
        bins=[x_bins, y_bins],
        statistic="sum",
    )
    return hist, x_bins, y_bins


def get_col_norm(hist):
    v_min = np.log10(hist[hist > 0].min())
    v_max = np.log10(hist.max())
    v_center = (v_max + v_min) / 2

    col_norm = colors.DivergingNorm(vmin=v_min, vcenter=v_center, vmax=v_max)
    return col_norm


def set_ax_params(ax, parameters):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(parameters["x_label"], size=parameters["x_labelsize"])
    ax.set_ylabel(parameters["y_label"], size=parameters["y_labelsize"])

    ax.tick_params(
        length=parameters["length_major_ticks"], width=parameters["width_major_ticks"]
    )
    ax.tick_params(
        length=parameters["length_minor_ticks"],
        width=parameters["width_minor_ticks"],
        which="minor",
    )
    return


def create_color_bar(ax, parameters, subfig, f):
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

    hist, x_bins, y_bins = get_hist(gas)

    f, ax = plt.subplots()
    col_norm = get_col_norm(hist)
    x_grid, y_grid = np.meshgrid(x_bins, y_bins)

    subfig = ax.pcolormesh(
        x_grid, y_grid, np.log10(hist.T), norm=col_norm, cmap=plt.get_cmap("coolwarm")
    )

    set_ax_params(ax, parameters)
    create_color_bar(ax, parameters, subfig, f)
    return


example_halos_1e7 = [(4779, 6, 0.84), (2606, 8, 0.556570), (519, 10, 0.643313)]
example_halos_1e8 = [(221, 6, 0.489819), (8, 10, 0.565650)]
