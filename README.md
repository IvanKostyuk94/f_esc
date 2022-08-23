# Ionizing photon production and escape fractions during cosmic reionization in the TNG50 simulation 

Toolkit used to analyze data and create plots of Illustris TNG-50 halos post-processed with the radiation transfer code CRASH to obtain the state of radiational equilibrium.

The findings are summarized in this <a href="https://arxiv.org/abs/2207.11278">paper</a>.

## Files

<ul>
    <li><b>build_df.py:</b> Goes through the folder of simulation runs and creates a <a href="https://pandas.pydata.org/">pandas</a> dataframe conataining the properties of the examined halos</li>
    <li><b>update_df.py:</b></li> Adds several additional columns to the halo dataframe namely the ages of stellar particles in the halo as well as the gas clumping, surface stellar density and surface gas density of the central galaxy, defined to be enclosed by 0.2Rvir.
    <li><b>build_radius_df.py:</b></li> Deprecated function originally used to analyze the escape fraction at different distances from the halo center. Should this become necessary in the future I would recommend to change <b>build_df.py</b> with elements from this file.
    <li><b>add_metallicities.py:</b></li> Adds an additional column to the dataframe containing stellar metallicities. Should be moved into <b>update_df.py</b> in the future
    
</ul>  

