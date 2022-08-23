# Ionizing photon production and escape fractions during cosmic reionization in the TNG50 simulation 

Toolkit used to analyze data and create plots of Illustris TNG-50 halos post-processed with the radiation transfer code CRASH to obtain the state of radiational equilibrium.

The findings are summarized in this <a href="https://arxiv.org/abs/2207.11278">paper</a>.

To use the notebooks the dataframes containing the halo information are needed. To get the dataframes please contact me under <a href = "mailto: ivkos@mpa-garching.mpg.de">ivkos@mpa-garching.mpg.de</a>.

## Python files

<i>These files interact directly with either the CRASH post-processing runs or TNG50 snapshot data</i>
<ul>
    <li><b>build_df.py:</b> Goes through the folder of simulation runs and creates a <a href="https://pandas.pydata.org/">pandas</a> dataframe conataining the properties of the examined halos</li>
    <li><b>update_df.py:</b></li> Adds several additional columns to the halo dataframe namely the ages of stellar particles in the halo as well as the gas clumping, surface stellar density and surface gas density of the central galaxy, defined to be enclosed by 0.2Rvir.
    <li><b>add_metallicities.py:</b></li> Adds an additional column to the dataframe containing stellar metallicities. Should be moved into <b>update_df.py</b> in the future
    <li><b>build_radius_df.py:</b></li> <b>Deprecated</b> code originally used to analyze the escape fraction at different distances from the halo center. Should this become necessary in the future I would recommend to change <b>build_df.py</b> with elements from this file.
    <li><b>fix_stellar_ages.py:</b></li> <b>Deprecated</b> code which was only used when <b>update_df.py</b> used the wrong stellar particles in a data frame update.
</ul>  

## Notebooks
<i>Most of these files mostly interact with the dataframe summarizing the results of the simulation and can be used idependently from the simulation runs</i>

<ul>
    <li><b>counts_histogram.ipynb:</b> Used to create plots summarizing the halo population as seen in fig. 1 of the paper.
    <li><b>esc_fraction.ipynb:</b> A python implementation of the code used to obtain escape fractions along lines of sight and create maps of densities and escape fractions of the halo as seen in fig. 5. Note that this code requires the full density and ionization maps and therefore needs the full simulation data of a halo. To be more readable this code needs to be refactored and most of it moved to a python file. 
    <li><b>fesc_vs_quant.ipynb:</b> Used for analyzing the average escape fraction as a function of different properties.
    <li><b>clean_up.ipynb:</b> Collection of short scripts to delete unnecessary files in the simulation runs in order to reduce memory usage.
    <li><b>halo_image.ipynb:</b> Used to create projected images of halos as seen in fig. 4. This notebook needs the full halo maps and therefore also the full simulation results.
    <li><b>large_radii.ipynb:</b> Used to analyze the escape fraction as a function of distance from the halo center. Not used in the final project.


    
