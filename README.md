### README ###
Author:  Sam Zipper
Contact: samuelczipper@gmail.com

This repository is home to code for the High Resolution Mapping of 
EvapoTranspiration (HRMET; pronounced "hermit") model.

HRMET is introduced and described in the following publication:

    Zipper, S.C. & S.P. Loheide II (2014). Using evapotranspiration to
    assess drought sensitivity on a subfield scale with HRMET, a high
    resolution surface energy balance model. Agricultural & Forest
    Meteorology 197: 91-102. DOI: 10.1016/j.agrformet.2014.06.009

Please cite this paper if you use HRMET or any derivative thereof.

## Repository Contents ##
HRMET.R = this is the HRMET code

Publications.txt = a text file tracking publications that use HRMET; if you 
                   publish something, let me know by submitting a push request! (or send an email)

HRMET_HowTo/* = Example HRMET files showing how to run HRMET over a grid with 
                spatially variable inputs and estimate uncertainty in results; 
                however, these are with a no-longer-maintained MATLAB version of 
                the code (the version used in the original HRMET paper) and are 
                not complete. Regardless, they may help you understand how to use 
                the current, R version of the code at your own study site.

## Key Assumptions of HRMET ##
-HRMET calculated the 1D surface energy balance. However, it is typically
 applied over fields to produce raster maps of ET. In order to do this, 
 you simply have to define the relevant inputs at all locations you want
 to map ET, and then run HRMET at each grid point. 
 
-Thus, assumptions of spatial homogeneity of inputs should be made with care. 
 For example, in Zipper et al. (2014), we assume uniform meteorological conditions
 over our relatively small (~600 x 600 m) field. This assumption gets increasingly
 problematic as your spatial scale increases. HRMET was designed for precision-agriculture
 scale applications; however, the physical principles should work at larger scales, so
 long as the the input data is sufficiently high-resolution.

## Known Issues ##
-HRMET does not work well in extremely short canopies or deserts (h approaching 0 m).

-HRMET does not work well when the canopy height (h) exceeds the height of 
 temperature and wind speed measurements (Zair and Zu).
 
-The G_Tw coefficient (used in cloudiness estimation) takes a summer value by default; 
 future versions should automatically select this based on DOY.
 
## Bug Fix History ##
2017-06-07: Fix small bug in cloudiness fraction calculation (it was being inadvertently set to 0 under most conditions)
