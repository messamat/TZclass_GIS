# TZclass_GIS

Purpose: These scripts include all steps used in producing the river classification for the Rufiji River basin and Tanzania as well as the 
		figures contained in the corresponding report and peer-reviewed publication.
		
All GIS analyses prior to these steps are available in messamat/TZclass_GIS and require an ESRI ArcGIS license (including the Spatial and Network Analyst extensions).

These scripts are annotated but could be challenging to follow. Please contact the author for comments and clarifications. 
Numbers cited in accompanying report and peer-reviewed publication can all be found in the scripts.

For the purpose of reproducing the full workflow used in this analysis, run the scripts in the following order:
1. Python_GIS/data_prep_TZ.py
2. Python_GIS/Rufiji_gages_envdata.py

Followed by the following scripts from the present directory, for the inductive classification: 
From https://github.com/messamat/TZclass
3. R/hydrodataprep.R
4. R/hydrodataprep2.R (optional)	
5. R/hydrodataprep3.R		
6. R/hydrodataprep4.R		
7. R/hydrodataprep5.R		
8. R/hydrostats_claspred.R
