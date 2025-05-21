# HydroFill-Example

Before investigating the software example, it is recommended that you look at the presentation slides "Eli_Robert_ASCE_2025_WE&WRC_Wide-Slide.pdf" and read the proceedings paper "ASCE_2025_WE&WRC_Paper_Robert_N_Eli.pdf". Also, this example repository is not static and will be updated with additional files and instructional information over time.

There are five MATLAB .m script files included in this example repository. Only three are required to preprocess the raw DEM to produce the fully drainable DEM for illustration purposes. It should be noted that the first .m script to be executed (West_Run_I68_MSB_DEM_Creation.m) combines the corrected original raw DEM with a Microsoft Buildings raster file to produce a modified raw DEM that includes the building structures. The second .m script to be executed (West_Run_I68_MSB_Culvert_Import.m) adds the culvert (or drain) inlets and outlets via a location table CIO and the raster file nbflag. The third .m script to be executed (HydroFill_West_Run_I68.m) produces the final fully drainable DEM that includes Microsoft Building structures and 35 culverts or drains.

It is mandatory that the .m scripts be executed in the proper sequence since intermediate files are generated that are required for input to the next .m script in the sequence. The remaining two .m script files provide illustration of the virtual ponding capability (HydroFill_Ponding_West_Run_I68.m) and the further processing of the fully drainable DEM to identify the drainage basin number bn and the accumulative drainage area di for all grid cells (HydroDrain_West_Run_I68.m). These latter two .m scripts can be executed after the required input files have been generated.

The following sequence lists each .m file and the corresponding input files required, and the resulting output files. Each .tif file is in fact a geoTIFF file that contains all the geospatial metadata necessary to view the file as a map layer in GIS software such as QGIS or ArcGIS.

1. West_Run_I68_MSB_DEM_Creation.m
   
   Input files:
   1. West_Run_I68_Fill_corrected.tif
   2. WV_MSB_WestRun_I68_clip.tif
  
   Output files:
   1. West_Run_I68_MSB_DEM.tif
   2. West_Run_I68_MSB_DEM_Creation.mat

2. West_Run_I68_MSB_Culvert_Import.m

   Input files:
   1. West_Run_I68_MSB_DEM_Creation.mat
   2. West_Run_I68_Culvert_IO.tif
  
   Output files:
   1. West_Run_I68_nbflag.tif
   2. West_Run_I68_MSB_Culvert_Import.mat

3. HydroFill_West_Run_I68.m

   Input files:
   1. West_Run_I68_MSB_DEM_Creation.mat
   2. West_Run_I68_MSB_Culvert_Import.mat
  
   Output files:
   1. HydroFill_West_Run_I68.mat
   2. HydroFill_West_Run_I68_wdepth.tif
   3. HydroFill_West_Run_I68_idem.tif

4. HydroFill_Ponding_West_Run_I68.m

   Input files:
   1. West_Run_I68_MSB_DEM_Creation.mat

   Output files:
   1. HydroFill_Ponding_West_Run_I68.mat
   2. HydroFill_Ponding_West_Run_I68_wdepth.tif
  
5. HydroDrain_West_Run_I68.m

   Input files:
   1. HydroFill_West_Run_I68.mat

   Output files:
   1. HydroDrain_West_Run_I68.mat
   2. HydroDrain_West_Run_I68_di.tif
   3. HydroDrain_West_Run_I68_bn.ti
