import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyTNG.cosmology as _cosmo


h = _cosmo.TNGcosmo.h

def set_ax_params(ax):
    ax.tick_params(
        length=16,
        width=4,
    )
    ax.tick_params(
        length=8,
        width=3,
        which="minor",
    )
    ax.tick_params(axis="both", which="both",labelsize=30,)

    # change all spines
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(3)
    return


def plot_temp_solver_comp(df_temp, df_no_temp):

	figsize = [15, 10]
	size = 150
	label_solver = 'with T solver'
	label_no_solver = 'without T solver'
	legendsize = 30
	labelsize = 40

	x_label = r'$\log \left(\frac{M_\star}{M_\odot} \right)$'
	y_label = r'$f_\mathrm{esc}$'


	df_merged = df_no_temp.merge(df_temp[['ID','z','f_esc']], on=['ID', 'z'], how='inner', suffixes=('_no_temp_solver', '_temp_solver'))
	df_merged['StarMass'] = np.log10(df_merged['HaloMass']*df_merged['FractionStars']*1e10/h)

	_, ax = plt.subplots(figsize=figsize)
	ax.scatter(df_merged['StarMass'], df_merged['f_esc_no_temp_solver'], 
			s=size, marker='x', label=label_no_solver, color='grey')
	ax.scatter(df_merged['StarMass'], df_merged['f_esc_temp_solver'], 
			s=size, marker='+', label=label_solver, color='black')
	
	ax.set_xlabel(x_label, size=labelsize)
	ax.set_ylabel(y_label, size=labelsize)
	set_ax_params(ax)
	ax.legend(fontsize=legendsize)
