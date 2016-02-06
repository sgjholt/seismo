# seismo
=========================================================================================================================================
Files for creating kik-net earthquake database.
=========================================================================================================================================
Important Scripts in the repository:
-----------------------------------------------------------------------------------------------------------------------------------------
**1)build_GMdb.py** 
  *Dependancies: process_kikhdr.py*

**2)GMPE_compare.py**

Function of main scripts:
-----------------------------------------------------------------------------------------------------------------------------------------
**1)** Builds a database from kik-net data of station locations, ground motions (PGA) and a range of distance metrics 
   (Rrup/Rjb/Rhyp/Repi etc) for a given earthqauke.
   Also builds a separate earthuake metadata file including location, depth, magnitude, strike, dip, fault type, crustal/slab. 
   *More parameters may be added in the future.* 

**2)** Plots data and a GMPE of choice  Note currently in very early stage of development. 
   Current functionality includes plotting zhao_06 GMPE versus data of a single component of seismograms.
   Calculate residuals values between GMPE of choice and perform statisitcal analysis.
   *Future functionality will include ability to plot multiple GMPES from different models and the ability to combine componants.
   Also to compare downhole with suface seismometer data.* 
   
   *Instructions for use to follow shortly.*
