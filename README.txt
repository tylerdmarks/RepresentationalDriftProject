The 'Goard Lab 2P pre-processing pipeline' folder contains software for processing 2-photon calcium imaging data in the form of multi-tif files. 
This pipeline (A > B > C) produces .mat files containing neural response DFF data in [neurons x frames] format. The versions used for this study are included here. 
An updated version of this code is available at: https://github.com/ucsb-goard-lab/Two-photon-calcium-post-processing
This folder also includes the software used for producing deconvolved spike data from the DF/F data (used for event detection in Figure 2).

The 'Preprocessing and Analysis platforms' folder contains software for turning .mat files from the intial pipeline (above) into further processed and more organized 
data structures to be used for analysis. These processed data structures are the format in which the data for this study is provided, for ease of use. 
The 'Stability_masterscript' script demonstrates this process.

The 'Main figures' folder contains scripts for the analyses performed and figures generated for the main figures. 

The 'Other figures' folder contains scripts and functions for analysis done for supplementary figures and otherwise. 

The 'helper functions' folder contains various functions utilized by the code in the other folders.

