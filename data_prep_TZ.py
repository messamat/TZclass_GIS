__author__ = 'Mathis Messager'
#Creation date: 2018/03/08
#Last updated: 2018/26/08

#Project: USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#(Contract No. EDH-I-00-08-00023-00 621-TO 12)

#Purpose: Format spatial data for analysis

#Bugs to correct
    #Compute max statistics
    #Re-run analysis for Protected Areas: first dissolving PAs, then re-running whole accumulation script
    #Re-run analysis for forest loss (missing a small portion of Tanzania): tiles were added to the data folder so re-running the script as is shouldn't be an issue.

import arcpy
from arcpy import env
from arcpy.sa import *
import ArcHydroTools
import os
import sys
import numpy
import numpy.lib.recfunctions
import math
import time
import glob
from collections import defaultdict

arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = False
if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.AddMessage("Checking out Spatial")
    arcpy.CheckOutExtension("Spatial")
else:
    arcpy.AddError("Unable to get spatial analyst extension")
    arcpy.AddMessage(arcpy.GetMessages(0))
    sys.exit(0)
#TO UPDATE
rootdir="F:/Tanzania/Tanzania/"
datadir= rootdir + "data/"
wd = rootdir + "results/"
outhydro = wd + "outhydro.gdb/"
env.workspace = wd

#Usual variables
scatch = outhydro + "catchgrid118"
scatchpoly = outhydro + "catchpoly118"
sline = outhydro + "streamnet118"
spoint = outhydro+'point118'

##################################################################################################################################
# FUNCTIONS
##################################################################################################################################
def WGS84_pixelarea(in_ras, out_wd=wd):
    #Careful, only works for square grids
    #Require math, arcpy
    in_ras = Raster(in_ras)
    arcpy.env.extent = in_ras
    arcpy.env.cellSize = in_ras
    ###############################################################
    # Generate latitude grid (from https://community.esri.com/thread/43907)
    # Calculate $$NROWS and $$NCOLS from current environment
    cellSize = float(env.cellSize)
    nrows = int((env.extent.YMax - env.extent.YMin) / cellSize)
    ncols = int((env.extent.XMax - env.extent.XMin) / cellSize)
    # Bill Huber's method for $$XMAP and $$YMAP: "1" flows "right", "64" (63+1) flows "up"
    print("Create constant raster")
    tmpg = CreateConstantRaster(1)
    #xmap = (arcpy.sa.FlowAccumulation(tmpg) + 0.5) * cellSize + env.extent.XMin
    print("Compute flow accumulation")
    ymap = (arcpy.sa.FlowAccumulation(tmpg + 63) + 0.5) * cellSize + env.extent.YMin
    ###############################################################
    #Compute pixel area (from http://www.jennessent.com/downloads/Graphics_Shapes_Online.pdf p61-62)
    WGS84_radius = 6371000.79000915
    rad = math.pi/180
    print("Compute pixel area")
    pixelheight = WGS84_radius*2*ATan2(SquareRoot(((math.sin(math.radians(cellSize/2))**2)+(math.sin(math.radians(0))**2)*Cos(rad*(ymap+cellSize/2))*Cos(rad*(ymap-cellSize/2)))),
          SquareRoot(1-((math.sin(math.radians(cellSize/2))**2)+(math.sin(math.radians(0))**2)*Cos(rad*(ymap+cellSize/2))*Cos(rad*(ymap-cellSize/2)))))
    pixelwidth = WGS84_radius*2*ATan2(SquareRoot(((math.sin(math.radians(0))**2)+(math.sin(math.radians(cellSize/2))**2)*Cos(rad*(ymap+cellSize/2))*Cos(rad*(ymap-cellSize/2)))),
          SquareRoot(1-((math.sin(math.radians(0))**2)+(math.sin(math.radians(cellSize/2))**2)*Cos(rad*(ymap+cellSize/2))*Cos(rad*(ymap-cellSize/2)))))
    pixelarea = pixelheight*pixelwidth
    print("Save pixel area raster")
    pixelarea.save(os.path.join(out_wd,'pixelarea.tif'))
    arcpy.ClearEnvironment('extent')
    arcpy.ClearEnvironment('cellSize')
    arcpy.Delete_management(pixelheight)
    arcpy.Delete_management(pixelwidth)
    arcpy.Delete_management(ymap)

def ras_catcount(in_zone_data, in_class_data, output_wd, out_table, scratch_wd="C:\\temp", pixel_area=outhydro + 'pixelarea'): #class_data must be in integer format
    if pixel_area is None:
        if not arcpy.Exists(scratch_wd + '\\pixelarea.tif'):
            print('Computing pixel area')
            WGS84_pixelarea(in_zone_data, scratch_wd)
            pixel_area=scratch_wd + '\\pixelarea.tif'

    arcpy.MakeTableView_management(in_class_data,'zonetab')
    #Iterate over every category in class data
    print('Creating scratch GDB')
    scratchgdb = scratch_wd + '\\scratch.gdb'
    arcpy.CreateFileGDB_management(scratch_wd, out_name='scratch.gdb')
    try:
        with arcpy.da.SearchCursor('zonetab', ['Value']) as cursor:
            row = next(cursor)
            print(row[0])
            classbool = Con(Raster(in_class_data)==row[0],pixel_area,0) #Create a raster of pixel area in cells with that category
            print("Boolean done!")
            scratchtable=os.path.join(scratchgdb, 'catcount{}'.format(row[0]))
            ZonalStatisticsAsTable(in_zone_data,'Value',classbool,scratchtable,"DATA","SUM")
            print("Zonal stats done!")
            arcpy.AlterField_management(scratchtable, field="SUM", new_field_name='SUM_{}'.format(row[0]))
            tabjoin=arcpy.da.TableToNumPyArray(scratchtable, ('Value', 'SUM_{}'.format(row[0])))
            print("Create new numpy array done!")
            for row in cursor:
                print(row[0])
                classbool = Con(Raster(in_class_data)==row[0],pixel_area,0) #Create a raster of pixel area in cells with that category
                print("Boolean done!")
                scratchtable=os.path.join(scratchgdb, 'catcount{}'.format(row[0]))
                ZonalStatisticsAsTable(in_zone_data,'Value',classbool,scratchtable,"DATA","SUM")
                print("Zonal stats done!")
                arcpy.AlterField_management(scratchtable, field="SUM", new_field_name='SUM_{}'.format(row[0]))
                print("Alter field done!")
                array=arcpy.da.TableToNumPyArray(scratchtable, ('Value', 'SUM_{}'.format(row[0])))
                print("Numpy array done!")
                try:
                    tabjoin = numpy.lib.recfunctions.join_by('Value', tabjoin, array, jointype='inner', usemask=False)
                    print("Join done!")
                except: #Can run into MemoryError if array gets too big
                    print("Join didn't work, might have run out of memory, output array and write another one")
                    arcpy.da.NumPyArrayToTable(tabjoin, os.path.join(output_wd,out_table+str(row[0])))
                    print("Array to table done!")
                    del tabjoin
                    tabjoin=arcpy.da.TableToNumPyArray(scratchtable, ('Value', 'SUM_{}'.format(row[0])))
                    print("Create new numpy array done!")
            arcpy.da.NumPyArrayToTable(tabjoin, os.path.join(output_wd,out_table+str(row[0])))
            print("Array to table done!")
    except:
        del cursor
        print("Delete row and cursor done!")
        arcpy.Delete_management(scratchgdb)
        del tabjoin
        print("Delete gdb and array done!")

##################################################################################################################################
# FORMAT HYDROSHEDS DATA
##################################################################################################################################
#Mosaic hydrologically conditioned DEM
demconfolder = "F:/Tanzania/Tanzania/data/HydroSHEDS_201802/HydroSHEDS/DEMcon/"
demconlist = []
for (dirpath, dirnames, filenames) in os.walk(demconfolder):
    if len(dirnames)==2:
        demconlist.append(os.path.join(dirpath,dirnames[1]))
arcpy.MosaicToNewRaster_management(demconlist, output_location=demconfolder, raster_dataset_name_with_extension="DEMcon", number_of_bands=1)
#Mosaic direction raster
dirfolder = "F:/Tanzania/Tanzania/data/HydroSHEDS_201802/HydroSHEDS/dir/"
directionlist = []
for (dirpath, dirnames, filenames) in os.walk(dirfolder):
    if len(dirnames)==2:
        directionlist.append(os.path.join(dirpath,dirnames[1]))
arcpy.MosaicToNewRaster_management(directionlist, output_location=dirfolder, raster_dataset_name_with_extension="dir", number_of_bands=1)
#Create watershed polygon for Tanzania
tzwtshd = outhydro + "tzwtshd"
ArcHydroTools.BatchWatershedDelineationForPolygons("gadm1_lev0_Tanzania.shp", "dir", tzwtshd)
#Clip flow dir to Tanzania's drainage area
dirtz = outhydro+"dir_tz"
arcpy.Clip_management(dirfolder+"dir", arcpy.Describe(tzwtshd).extent, dirtz, in_template_dataset= tzwtshd,
                                clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Tends to crash
#Flow acc
flowacc = FlowAccumulation(dirtz, data_type="INTEGER")
flowaccras = outhydro+"flowacc"
flowacc.save(flowaccras)

#Delineate streams
#At 6.5 deg S, the approximate center of the raster, a degree of latitude is 110589 m, and a degree of long is 110609 m
#So 1 km2 is 1000000/(110589*110609*0.00083333333)= 118 cells
sras = outhydro+"streamras118"
streamras = Con(flowacc > 117, Raster(flowacc))
streamras.save(sras)
#Segment streams
sseg = outhydro+"streamseg118"
ArcHydroTools.StreamSegmentation(sras,dirtz, sseg)
#N.B: Important to have the double backward slash for path, otherwise tool won't run
#Generate streams polyline
sline = outhydro+"streamnet118"
ArcHydroTools.DrainageLineProcessing(sseg, dirtz, sline)
#Delineate subcatchments
scatch = outhydro + "catchgrid118"
ArcHydroTools.CatchmentGridDelineation(dirtz, sseg, scatch)
#Generate subcatchment polygons
scatchpoly = outhydro + "catchpoly118"
ArcHydroTools.CatchmentPolyProcessing(scatch, scatchpoly)
#A few polygons do not have corresponding streams, so this might be problematic -- to investigate down the line
#Generate drainage points for each subcatchments/streamline
spoint = outhydro + "point118"
arcpy.FeatureVerticesToPoints_management(sline, spoint, point_location='END')

#Compute flow length
"""Projecting the direction raster changes the values of too many cells and thus changes the flow direction overall.
However, flow length is not exact as it gives a measure in decimal degrees. If time allows, should create a weight
raster that will be equal to the distance across a cell given its latitude and the flow direction with raster
calculator"""

#########################################################################################################################
# FORMAT ELEVATION DATA
########################################################################################################################
SRTMfolder = "F:/Tanzania/Tanzania/data/SRTMGL1/unzipped/"
demlist = []
for (dirpath, dirnames, filenames) in os.walk(SRTMfolder):
    if len(filenames) != 0:
       demlist.append(os.path.join(dirpath,filenames[0]))
arcpy.MosaicToNewRaster_management(demlist, output_location="F:/Tanzania/Tanzania/data/SRTMGL1/", raster_dataset_name_with_extension="SRTMmosaic", number_of_bands=1)
#Clean erroneous values
SRTMmosaic= Raster("F:/Tanzania/Tanzania/data/SRTMGL1/SRTMmosaic")
SRTMmosaic_nodata = SetNull(SRTMmosaic, SRTMmosaic, "VALUE > 6000")
SRTMmosaic_nodata.save("F:/Tanzania/Tanzania/results/SRTMnodata")
SRTMmosaic_filled = EucAllocation(SRTMmosaic_nodata, maximum_distance="0.00416666665", out_distance_raster='SRTMfilleddis') #5-cells maximum fill-in
SRTMmosaic_filled.save("srtmfilled")
SRTMmosaic_clean = Con(IsNull(SRTMmosaic), SRTMmosaic, SRTMmosaic_filled) #To redefine coast as EucAlloc extended coast by 5 cells
SRTMmosaic_clean.save("srtmclean")
#Re-sample SRTM3 elevation at resolution and extent of catchments (conditioned DEM form hydrosheds is not suitable to compute
#slopes and elevation, aspect, etc. as many features are burnt and thus increase slope
cs = (arcpy.GetRasterProperties_management(dirtz, 'CELLSIZEX')).getOutput(0)
arcpy.env.snapRaster = scatch
arcpy.Resample_management(in_raster=wd+"srtmclean",out_raster=wd+'srtmclean90', cell_size=cs, resampling_type='BILINEAR')
#with Jess Jenness slope tool for WGS84 data: http://www.jennessent.com/arcgis/surface_area.htm, compute:
    # - Slope (Horn's Method): "/results/slopedeg
    # - Aspect (Horn's Method) : "/results/aspect"
    # - Profile curvature : "/results/profilcurvat"

########################################################################################################################
#Get catchment characteristics
########################################################################################################################
gdbname_cat = "catch_attri.gdb"
#arcpy.CreateFileGDB_management(wd, gdbname_cat)
env.workspace = wd + gdbname_cat

#########################################
# Drainage characteristics and topography
#Drainage area
arcpy.AddGeometryAttributes_management(scatchpoly, Geometry_Properties="AREA_GEODESIC", Area_Unit = "SQUARE_KILOMETERS")
arcpy.AlterField_management(scatchpoly, field="AREA_GEO", new_field_name="CatArea")
#Reach length
arcpy.AddGeometryAttributes_management(sline, Geometry_Properties="LENGTH_GEODESIC", Length_Unit = "KILOMETERS")
arcpy.AlterField_management(sline, field="LENGTH_GEO", new_field_name="CatLen")
#Flow accumulation
ZonalStatisticsAsTable(scatch, "Value", "flowacc",out_table="flowacc",ignore_nodata="DATA", statistics_type="MAXIMUM")
arcpy.AlterField_management("flowacc", "MAX", new_field_name="CatFlowAcc")
#Assign river order
ArcHydroTools.AssignRiverOrder(sline, "ReaOrd", 'Strahler', outhydro + "\\streamnet118_FS")
#Elevation
ZonalStatisticsAsTable(scatch, "Value", wd+"srtmclean90",out_table="elv",ignore_nodata="DATA", statistics_type="MIN_MAX_MEAN")
arcpy.AddField_management("elv", "CatElvRan", "SHORT")
arcpy.CalculateField_management("elv", "CatElvRan", "MAX-MIN")
arcpy.AlterField_management("elv","MEAN","CatElvAvg")
arcpy.AlterField_management("elv","MIN","CatElvMin")
arcpy.AlterField_management("elv","MAX","CatElvMax")
#Slope
ZonalStatisticsAsTable(scatch, "Value", wd+"slopedeg",out_table="slope",ignore_nodata="DATA", statistics_type="ALL")
arcpy.AlterField_management("slope","MEAN","CatSloAvg")
arcpy.AlterField_management("slope","STD","CatSloStd")
#Curvature
ZonalStatisticsAsTable(scatch, "Value", wd+"profilcurvat",out_table="profilcurvat",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("profilcurvat","MEAN","CatCurvAvg")
#Majority on flow dir (rather than aspect as an average on aspect would be meaningless)
ZonalStatisticsAsTable(scatch, "Value", dirtz,out_table="dir",ignore_nodata="DATA", statistics_type="MAJORITY")
arcpy.AlterField_management("dir","MAJORITY","CatDirMaj")
#Tabulate area on flow dir for accumulation over watersheds
ras_catcount(in_zone_data=scatch, in_class_data=dirtz, output_wd=wd+gdbname_cat, out_table='dirtabulate',
             scratch_wd="C:\\temp", pixel_area=outhydro + 'pixelarea')
for f in arcpy.ListFields('dirtabulate'):
    if 'SUM' in f.name:
       arcpy.AlterField_management("dirtabulate",f.name,f.name.replace("SUM","DirSum"))

#########################################
#Pekel data
#Water maximum extent
arcpy.ClearEnvironment("snapRaster")
maxextentras = "{0}{1}_20E_0N.tif;{0}{1}_30E_0N.tif;{0}{1}_30E_10N.tif;{0}{1}_30E_10S.tif".format(datadir + "Pekel2016\\","extent")
arcpy.MosaicToNewRaster_management(maxextentras, "F:\\Tanzania\\Tanzania\\results\\",wd+"pekelextent.tif", "",
                                   "8_BIT_UNSIGNED", "", "1", "LAST", "FIRST") #Export in GRID format didn't work
#Resample the catchment grid to match Pekel's resolution and alignment
pekeltemplate=wd+"pekelextent.tif"
cs = (arcpy.GetRasterProperties_management(pekeltemplate, 'CELLSIZEX')).getOutput(0)
arcpy.env.snapRaster = pekeltemplate
scatchpekel = outhydro+"\\catchgrid118_pekelresample"
arcpy.Resample_management(in_raster=scatch, out_raster=scatchpekel,cell_size=cs, resampling_type='NEAREST')
ZonalStatisticsAsTable(scatchpekel, "Value",wd+'pekelextent.tif', out_table="pekelextent",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("pekelextent","MEAN","CatWatExt")
arcpy.Delete_management(scatchpekel)
arcpy.Delete_management(wd+"pekelextent.tif")
#Water cover occurrence
maxchangeras = "{0}{1}_20E_0N.tif;{0}{1}_30E_0N.tif;{0}{1}_30E_10N.tif;{0}{1}_30E_10S.tif".format(datadir + "Pekel2016\\","occurrence")
arcpy.MosaicToNewRaster_management(maxchangeras, "F:\\Tanzania\\Tanzania\\results\\","pekeloccur.tif", "", "8_BIT_UNSIGNED", "", "1", "LAST", "FIRST")
ZonalStatisticsAsTable(scatchpekel, "Value",wd+'pekeloccur.tif', out_table="pekeloccur",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("pekeloccur","MEAN","CatWatOcc")
arcpy.Delete_management(wd+"pekeloccur.tif")
#Water cover change
maxchangeras = "{0}{1}_20E_0N.tif;{0}{1}_30E_0N.tif;{0}{1}_30E_10N.tif;{0}{1}_30E_10S.tif".format(datadir + "Pekel2016\\","change")
arcpy.MosaicToNewRaster_management(maxchangeras, "F:\\Tanzania\\Tanzania\\results\\","pekelchange.tif", "", "8_BIT_UNSIGNED", "", "1", "LAST", "FIRST")
pekelchange = Raster(wd+"pekelchange.tif")
pekelchange_format = Con(pekelchange < 253, pekelchange)
pekelchange_format.save(wd+"pekelchange_format.tif")
ZonalStatisticsAsTable(scatchpekel, "Value",wd+'pekelchange_format.tif', out_table="pekelchange",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("pekelchange","MEAN","CatWatCha")
arcpy.Delete_management(pekelchange)
arcpy.Delete_management(wd+"pekelchange_format.tif")
#Water cover seasonality
maxchangeras = "{0}{1}_20E_0N.tif;{0}{1}_30E_0N.tif;{0}{1}_30E_10N.tif;{0}{1}_30E_10S.tif".format(datadir + "Pekel2016\\","seasonality")
arcpy.MosaicToNewRaster_management(maxchangeras, "F:\\Tanzania\\Tanzania\\results\\","pekelseason.tif", "", "8_BIT_UNSIGNED", "", "1", "LAST", "FIRST")
ZonalStatisticsAsTable(scatchpekel, "Value",wd+'pekelseason.tif', out_table="pekelseason",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("pekelseason","MEAN","CatWatSea")
arcpy.Delete_management(wd+"pekelseason.tif")
arcpy.ClearEnvironment("snapRaster")

#########################################
#WorldClim and derived products data
#Bioclim
#scatch = outhydro + "\\catchgrid118"
bioclimdir = datadir + "WorldClim\\wc2_0_30s_bio"
bioclimlist = [filenames for (dirpath, dirnames, filenames) in os.walk(bioclimdir)][0]
biocliplist = [wd+"bio{}clip.tif".format(biolyr[-6:-4]) for biolyr in filenames]
for biolyr in bioclimlist:
    arcpy.Clip_management(os.path.join(bioclimdir, biolyr), arcpy.Describe(tzwtshd).extent, wd+"bio{}clip.tif".format(biolyr[-6:-4]),
                          in_template_dataset= tzwtshd,clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Tends to crash
arcpy.env.snapRaster = scatch
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
for biolyr in biocliplist:
    outbiolyr=wd+"bio{}_3s.tif".format(biolyr[-10:-8])
    if not arcpy.Exists(outbiolyr):
        arcpy.Resample_management(biolyr, outbiolyr, cs, resampling_type='CUBIC')
    print(outbiolyr)
    ZonalStatisticsAsTable(scatch, "Value", outbiolyr, "bio{}_3s".format(biolyr[-10:-8]),ignore_nodata="DATA",statistics_type="MEAN")
    arcpy.AlterField_management("bio{}_3s".format(biolyr[-10:-8]), "MEAN", "CatBio{}Av".format(biolyr[-10:-8]))
    print("stats done!")
    if arcpy.Exists(outbiolyr):
        arcpy.Delete_management(outbiolyr)

#Solar radiation
solardir = datadir + "WorldClim\\wc2_0_30s_srad"
solarlist = [os.path.join(solardir,files) for files in [filenames for (dirpath, dirnames, filenames) in os.walk(solardir)][0][1:]]
solarannual = CellStatistics(solarlist, "SUM")
solarannual.save(solardir + "\\solarannual.tif")
arcpy.Clip_management(solardir + "\\solarannual.tif", arcpy.Describe(tzwtshd).extent, wd+"solarannual.tif",
                      in_template_dataset= tzwtshd,clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Crashes
arcpy.Resample_management(wd + "solarannualclip.tif", wd + "solarannual3s.tif", cs, resampling_type='BILINEAR')
ZonalStatisticsAsTable(scatch, "Value", wd + "solarannual3s.tif", "solarannual",ignore_nodata="DATA",statistics_type="MEAN")
arcpy.AlterField_management("solarannual","MEAN","CatSolAvg")
arcpy.Delete_management(solardir + "\\solarannual.tif")
#Global PET
PET3s = wd + "pet3s.tif"
PET = datadir + "GlobalPET\\PET_he_annual\\pet_he_yr"
arcpy.Clip_management(PET, arcpy.Describe(tzwtshd).extent, wd+"petclip.tif",
                      in_template_dataset= tzwtshd,clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Crashes
arcpy.Resample_management(wd+"petclip.tif", PET3s, cs, resampling_type='NEAREST')
ZonalStatisticsAsTable(scatch, "Value", PET3s, "pet",ignore_nodata="DATA",statistics_type="MEAN")
arcpy.AlterField_management("pet","MEAN","CatPETAvg")
arcpy.Delete_management(PET3s)
#Global Aridity Index
AI3s = wd + "aridity.tif"
AI = datadir + "GlobalAridity\\AI_annual\\ai_yr"
arcpy.Clip_management(AI, arcpy.Describe(tzwtshd).extent, wd+"aiclip.tif",
                      in_template_dataset= tzwtshd,clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Crashes
arcpy.Resample_management(wd+"aiclip.tif", AI3s, cs, resampling_type='NEAREST')
ZonalStatisticsAsTable(scatch, "Value", AI3s, "ai",ignore_nodata="DATA",statistics_type="MEAN")
arcpy.AlterField_management("ai","MEAN","CatAIAvg")
arcpy.Delete_management(AI3s)
#Rainfall Erosivity
erosivity = wd + "erosivity.tif"
glo_erosivity = datadir + "GloREDa\\GlobalR_NoPol.tif"
arcpy.Clip_management(glo_erosivity , arcpy.Describe(tzwtshd).extent, wd+"erosivityclip.tif",
                      in_template_dataset= tzwtshd,clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Crashes
arcpy.Resample_management(wd+"erosivityclip.tif", erosivity, cs, resampling_type='NEAREST')
ZonalStatisticsAsTable(scatch, "Value", erosivity, "erosivity",ignore_nodata="DATA",statistics_type="MEAN")
arcpy.AlterField_management("erosivity","MEAN","CatEroAvg")
arcpy.Delete_management(erosivity)
arcpy.ClearEnvironment("snapRaster")

###############################################
#Percentage watershed drained from lake
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
arcpy.env.extent = arcpy.Describe(scatch).extent
hydrolakes = datadir + "HydroLAKES\\HydroLAKES_polys_v10.gdb\\HydroLAKES_polys_v10"
#Convert lake polygons to raster
lakeras = wd+"hydrolakes"
arcpy.FeatureToRaster_conversion(hydrolakes, field="Hylak_id", out_raster=lakeras, cell_size=cs)
#Crop flow accumulation to lake boundaries
lakeaccras = wd + "lakeacc"
lakeacc = Con(Raster(lakeras) > 0, Raster(flowaccras),0)
lakeacc.save(lakeaccras)
#Create catchment lake pourpoint raster
catchmaxlakacc = wd + "catmaxlakacc"
maxacc = ZonalStatistics(in_zone_data=scatch,zone_field="Value", in_value_raster=lakeaccras, statistics_type="MAXIMUM",ignore_nodata="DATA")
maxacc.save(catchmaxlakacc)
catlakpoint =wd + "catlakpoint"
pour = Con(Raster(catchmaxlakacc) == Raster(lakeaccras), Raster(lakeaccras), 0)
pour.save(catlakpoint)
catlakid =wd + "catlakid"
lakid = Con(Raster(catlakpoint)>0, Raster(lakeras))
lakid.save(catlakid)
#Extract values to table
ZonalStatisticsAsTable(in_zone_data=scatch, zone_field="Value", in_value_raster=catlakpoint, out_table=wd + 'catlakacc.dbf', ignore_nodata='DATA', statistics_type='MAXIMUM')
ZonalStatisticsAsTable(in_zone_data=scatch, zone_field="Value", in_value_raster=catlakid, out_table=wd + 'catlakid.dbf', ignore_nodata='DATA', statistics_type='MAJORITY')
#Output catchment CatLakInd
arcpy.CopyRows_management(wd + 'catlakacc.dbf', 'lakeacc')
arcpy.AlterField_management(in_table="lakeacc",field="MAX",new_field_name="CatLakInd",new_field_alias='CatLakInd')
#Join tables
arcpy.MakeTableView_management(wd + 'catlakacc.dbf', 'catlakacc')
arcpy.AddJoin_management('catlakacc', 'Value', wd + 'catlakid.dbf','Value')
catlakacid= wd +'catlakaccid.dbf'
arcpy.CopyRows_management('catlakacc', wd + 'catlakaccid.dbf')
arcpy.Delete_management(wd + 'catlakacc.dbf')
arcpy.Delete_management(wd + 'catlakid.dbf')
#Only keep records of downstream-most intersecting river segment for each lake
duplitab = wd + catlakacid[:-4]+'_identical.dbf'
arcpy.FindIdentical_management(catlakacid,duplitab,fields='MAJORITY',output_record_option='ONLY_DUPLICATES')
dupliID = [id[0] for id in arcpy.da.SearchCursor(duplitab, ['IN_FID'])]
[f.name for f in arcpy.ListFields(catlakacid)]
duplicates = [[row[0],row[1],row[2]] for row in arcpy.da.SearchCursor(catlakacid, ['OID','MAJORITY','MAX']) if row[0] in dupliID]
if len(duplicates)>0:
    d={}
    ldel=[]
    for sub in duplicates: #Inspired from https://stackoverflow.com/questions/34334381/removing-duplicates-from-a-list-of-lists-based-on-a-comparison-of-an-element-of
        k=sub[1]
        if k not in d.keys():
            d[k]=sub
        elif sub[2] > d[k][2]:
            ldel.append(d[k][0])
            d[k]=sub
        else:
            #print(sub[0])
            ldel.append(sub[0])
expr= 'NOT "OID" IN ' + str(tuple(ldel))
arcpy.MakeTableView_management(catlakacid, out_view='catlakacid_view', where_clause=expr)
arcpy.CopyRows_management('catlakacid_view',wd+'catlakac_sub.dbf')
arcpy.Delete_management(duplitab)
arcpy.Delete_management(lakeras)
arcpy.ClearEnvironment("extent")

#Percentage watershed drained from reservoir #WOULD HAVE TO REPRODUCE WORKFLOW ABOVE TO CONTINUE CODE BLOCK BELOW
# hydrolakes = datadir + "HydroLAKES\\HydroLAKES_polys_v10.gdb\\HydroLAKES_polys_v10"
# arcpy.MakeFeatureLayer_management(hydrolakes, 'reservoirlyr')
# arcpy.SelectLayerByAttribute_management('reservoirlyr', "NEW_SELECTION", 'Lake_type > 1')
# reservoiras = wd+"hydrores"
# arcpy.FeatureToRaster_conversion('reservoirlyr', field="Hylak_id", out_raster=reservoiras, cell_size=cs)
# resaccras = wd + "resacc.tif"
# resacc = Con(Raster(reservoiras) > 0, flowaccras,0)
# resacc.save(resaccras)

#####################################################
#Lithology
geology = datadir + "GMIS\\geology\\geol\\geology.shp"
geologyras= wd + 'georas.tif'
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
arcpy.env.extent = arcpy.Describe(scatch).extent
arcpy.PolygonToRaster_conversion(geology,"symbol_nr",geologyras,"MAXIMUM_AREA",priority_field=None,cellsize=cs) #crashes every time
#Majority of coverage within each catchment
ZonalStatisticsAsTable(scatch, "Value", geologyras,out_table="geology",ignore_nodata="DATA", statistics_type="MAJORITY")
arcpy.AlterField_management("geology", "MAJORITY", "CatGeoMaj")
#Tabulate area on geology for watershed calculations
ras_catcount(in_zone_data=scatch, in_class_data=geologyras, output_wd=wd+gdbname_cat, out_table='geologytabulate',
             scratch_wd="C:\\temp", pixel_area=outhydro + 'pixelarea')
for tab in arcpy.ListTables('geologytabulate*'):
    for f in arcpy.ListFields(tab):
        if 'SUM' in f.name:
            arcpy.AlterField_management(tab,f.name,f.name.replace("SUM","GeolSum"))
arcpy.Delete_management(geologyras)
arcpy.ClearEnvironment('extent')

#####################################################
#Soil
soilpoly= datadir + "Africa_soil_WRB\\africasoilmap.shp"
arcpy.MakeTableView_management(soilpoly,out_view='soiltabview')
[f.name for f in arcpy.ListFields('soiltabview')]
arcpy.Frequency_analysis('soiltabview',wd+'soilfreqafrica2.dbf',frequency_fields='SU_WRB1_PH') #Create a unique number ID for each soil type
arcpy.MakeFeatureLayer_management(soilpoly,'soilyr')
arcpy.AddJoin_management('soilyr', 'SU_WRB1_PH', wd+'soilfreqafrica2.dbf', join_field='SU_WRB1_PH')
soilpolyfreq=wd+'africasoilmap_freq.shp'
arcpy.CopyFeatures_management('soilyr', soilpolyfreq)
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
arcpy.env.extent = arcpy.Describe(scatch).extent
soilras= wd + 'soilras'
[f.name for f in arcpy.ListFields(soilpolyfreq)]
arcpy.PolygonToRaster_conversion(soilpolyfreq,"OID_",soilras,"MAXIMUM_AREA",cellsize=cs) #crashes every time
soilnowaterras = wd + 'soilnowater'
soilnowater = Con(Raster(soilras)>0,soilras) #Get rid of seawater
soilnowater.save(soilnowaterras)
#Majority of coverage within each catchment
ZonalStatisticsAsTable(scatch, "Value", soilnowaterras,out_table="soil",ignore_nodata="DATA", statistics_type="MAJORITY")
arcpy.AlterField_management("soil", "MAJORITY", "CatSoilMaj")
#Tabulate area on soil for watershed calculations
ras_catcount(in_zone_data=scatch, in_class_data=soilnowaterras, output_wd=wd+gdbname_cat, out_table='soiltabulate',
             scratch_wd="C:\\temp", pixel_area=outhydro + 'pixelarea')
for tab in arcpy.ListTables('soiltabulate*'):
    for f in arcpy.ListFields(tab):
        if 'SUM' in f.name:
            arcpy.AlterField_management(tab,f.name,f.name.replace("SUM","SoilSum"))
arcpy.Delete_management(soilras)
arcpy.Delete_management(soilpolyfreq)
arcpy.ClearEnvironment('extent')

#######################################################
#Depth to bedrock
depthrock = datadir+'SoilGrids\\BDTICM_M_250m.tif'
depthrock3s = wd + "BDTICM_M_3s.tif"
arcpy.Clip_management(depthrock, arcpy.Describe(tzwtshd).extent, wd+"BDTICM_M_250mclip.tif",
                      in_template_dataset= tzwtshd,clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Crashes
arcpy.env.snapRaster = scatch
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
arcpy.Resample_management(wd+"BDTICM_M_250mclip.tif", depthrock3s, cs, resampling_type='BILINEAR')
ZonalStatisticsAsTable(scatch, "Value", depthrock3s, "depthrock",ignore_nodata="DATA",statistics_type="MEAN")
arcpy.AlterField_management("depthrock", "MEAN", "CatDRocAvg")
arcpy.Delete_management(depthrock3s)

#######################################################
#Permeability and porosity
glhymps = datadir + "GLHYMPS\\GLHYMPS.gdb\\Final_GLHYMPS_Polygon"
glhympsproj = wd+'glhympsproj.shp'
sr84 = 4326 #Project to WGS84
arcpy.Project_management(glhymps , glhympsproj, sr84)
porosityras= wd + 'porosityras'
permeabras = wd + 'permeabras'
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
arcpy.env.extent = arcpy.Describe(scatch).extent
arcpy.PolygonToRaster_conversion(glhymps,"Porosity",porosityras,"MAXIMUM_AREA",cellsize=cs)
arcpy.PolygonToRaster_conversion(glhympsproj,"Permeability_no_permafrost",permeabras,"MAXIMUM_AREA",cellsize=cs)
ZonalStatisticsAsTable(scatch, "Value", porosityras,out_table="porosity",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("porosity", "MEAN", "CatPoroAvg")
ZonalStatisticsAsTable(scatch, "Value", permeabras,out_table="permeability",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("permeability", "MEAN", "CatPermAvg")
arcpy.Delete_management(porosityras)
arcpy.Delete_management(permeabras)
arcpy.ClearEnvironment('extent')

########################################################
#Land cover
LC = datadir + 'Sentinel2\\ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif'
cs = (arcpy.GetRasterProperties_management(LC, 'CELLSIZEX')).getOutput(0)
scatchLC = outhydro+"catchgrid118_LCresample"
arcpy.Resample_management(in_raster=scatch, out_raster=scatchLC,cell_size=cs, resampling_type='NEAREST')
LCclip = wd+'lcclip'
arcpy.Clip_management(LC, arcpy.Describe(scatchLC).extent, LCclip,
                      in_template_dataset= scatchLC,clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT") #Crashes
LCcleanras = wd+'lcclean.tif'
LCclean= Con(Raster(LCclip)<100,Raster(LCclip)) #Clean out glitchy pixels with a value of 200
LCclean.save(LCcleanras)
ras_catcount(in_zone_data=scatchLC, in_class_data=LCcleanras, output_wd=wd+gdbname_cat, out_table='landcovertabulate',
             scratch_wd="C:\\temp", pixel_area=outhydro + 'pixelarea')
for tab in arcpy.ListTables('landcovertabulate*'):
    for f in arcpy.ListFields(tab):
        if 'SUM' in f.name:
            arcpy.AlterField_management(tab,f.name,f.name.replace("SUM","LCSum"))
arcpy.Delete_management(LCclip)
arcpy.Delete_management(scatchLC)

########################################################
#Mines
mines = datadir + "GMIS\\mines\\mines.shp"
arcpy.MakeFeatureLayer_management(mines, 'mineslyr')
arcpy.SelectLayerByAttribute_management('mineslyr', "NEW_SELECTION",'NOT "miningexpl" = {}'.format("'4 - Prospect inactive'"))
#arcpy.CopyFeatures_management('mineslyr', wd+'minesedit.shp') #Started to satellite imagery to place points more precisely and estimate footprint but did not have enough time
minesbuffer = wd+'minesbuffer.shp'
arcpy.Buffer_analysis('mineslyr', minesbuffer, '500 meters')
minesbufras= wd+'minesbufras.tif'
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
arcpy.env.snapRaster = scatch
arcpy.PolygonToRaster_conversion(minesbuffer, 'FID', minesbufras, "MAXIMUM_AREA", cellsize=cs)
minesbusrasformat = Con(IsNull(Raster(minesbufras)), 0, Raster(minesbufras))
ZonalStatisticsAsTable(scatch, "Value",minesbusrasformat, out_table="mines",ignore_nodata="DATA", statistics_type="VARIETY")
arcpy.AlterField_management("mines", "VARIETY", "CatMineDen")
arcpy.ClearEnvironment("snapRaster")

########################################################
#Dams
damsTZ= datadir + "TZ_OpenDataPortal20180221\\XYalldamsconstructedbyninewaterbasinsboard\\XYalldamsconstructedbyninewaterbasinsboard.shp"
GranD = datadir + 'GranD_version_1_1\\GranD_dams_v1_1.shp'
arcpy.MakeFeatureLayer_management(GranD, 'grandlyr')
arcpy.SelectLayerByAttribute_management('grandlyr', 'NEW_SELECTION', "{0}=4042 OR {0}=4501 OR {0}=4502".format("GRAND_ID")) #Select only those dams in TZ that are not already in the TZ dataset
#Merge with GranD 2 dams
dams=wd+'damsjoin.shp'
arcpy.Merge_management([damsTZ, 'grandlyr'], dams)
arcpy.AddField_management(dams, 'Count', 'SHORT')
arcpy.CalculateField_management(dams, 'Count', 1)
damscatchinters=wd+'damscatchinters.shp'
arcpy.Intersect_analysis([dams, scatchpoly],damscatchinters)
arcpy.Statistics_analysis(damscatchinters, "dams", [["Count", "SUM"]],case_field="GridID")
arcpy.AlterField_management("dams", "SUM_COUNT", "CatDamDen")

#########################################################
#Population - AfriPop data
afripopdir=datadir+"Afripop\\"
afripoplist = [fn for fn in glob.glob(afripopdir+'*15*.tif') if fn not in glob.glob(afripopdir+'*pph*.tif') if fn not in glob.glob(afripopdir+'*adj*')]
arcpy.MosaicToNewRaster_management(afripoplist, output_location=afripopdir, raster_dataset_name_with_extension="afripopmos", number_of_bands=1)
afripopcor = Con(IsNull(Raster(afripopdir+"afripopmos")), 0,Raster(afripopdir+"afripopmos"))
afripopcor.save(afripopdir+'afripopcor')
cs = (arcpy.GetRasterProperties_management(scatch, 'CELLSIZEX')).getOutput(0)
pop3s = wd+"pop3s"
arcpy.env.snapRaster = scatch
arcpy.Resample_management(in_raster=afripopdir+"afripopcor", out_raster=pop3s,cell_size=cs, resampling_type='NEAREST')
popdens= 1000000*Raster(pop3s)/Raster(outhydro+'pixelarea')
popdens.save(wd+"popdens")
ZonalStatisticsAsTable(scatch, "Value",popdens, out_table="afripopdens2",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("afripopdens2", "MEAN", "CatPopDen")
arcpy.Delete_management(afripopdir+'afripopmos')
arcpy.Delete_management(afripopdir+'afripopcor')
arcpy.Delete_management(pop3s)
arcpy.ClearEnvironment('snapRaster')

#########################################################
#Forest cover loss
forestlossdir=datadir+'HansenGFC2016\\'
forestlosslist = glob.glob(forestlossdir+'*lossyear*.tif')
arcpy.MosaicToNewRaster_management(forestlosslist, wd,"forestloss.tif", "",
                                   "8_BIT_UNSIGNED", "", "1", "LAST", "FIRST") #Export in GRID format didn't work
forestlossbool=Con(Raster(wd+"forestloss.tif")>0,1,0)
forestlossboolras = wd+"forestlossbool.tif"
forestlossbool.save(forestlossboolras)
#Resample the catchment grid to match resolution and alignment
cs = (arcpy.GetRasterProperties_management(forestlossboolras, 'CELLSIZEX')).getOutput(0)
arcpy.env.snapRaster = forestlossboolras
scatchfl = outhydro+"catchgrid118_forestloss"
arcpy.Resample_management(in_raster=scatch, out_raster=scatchfl,cell_size=cs, resampling_type='NEAREST')
ras_catcount(scatchfl, forestlossboolras, output_wd=wd+gdbname_cat,out_table="forestloss",scratch_wd="C:\\temp",pixel_area=wd+"pixelarea.tif")
for f in arcpy.ListFields("forestloss"):
    if 'SUM' in f.name:
        arcpy.AlterField_management("forestloss",f.name,f.name.replace("SUM","CatFLosSum"))
arcpy.Delete_management(wd+"forestloss.tif")
arcpy.Delete_management(scatchfl)

###########################################################
#Roads
#Include all main roads and roads links types, except residential, as many are just dirt roads and likely more inconsistently mapped
osmdir = datadir+'osm'
roadlist = glob.glob(osmdir+"/*/*roads*")
railwaylist = glob.glob(osmdir+"/*/*railways*")
#Take out the dot in the shapefile names and only keep forward slash as otherwise cause error in Merge
for roadlyr in roadlist:
    if "gis." in roadlyr:
        arcpy.Rename_management(roadlyr, roadlyr.replace("gis.",""))
for raillyr in railwaylist:
    if "gis." in raillyr:
        arcpy.Rename_management(raillyr, raillyr.replace("gis.",""))
roadlist = [roadlyr.replace("\\","/") for roadlyr in glob.glob(os.path.join(osmdir,"*/*roads*.shp"))]
railwaylist = [raillyr.replace("\\","/") for raillyr in glob.glob(os.path.join(osmdir,"*/*railways*.shp"))]
arcpy.Merge_management(roadlist,os.path.join(osmdir,"roads.shp"))
arcpy.MakeFeatureLayer_management(os.path.join(osmdir,"roads.shp"), "roadslyr")
expr = "{0}='motorway' OR {0}='motorway_link' OR {0}= 'primary' OR {0}= 'primary_link' OR {0}= 'secondary' OR {0}= 'secondary_link' OR {0}= 'service' OR {0}= 'tertiary' OR \
{0}= 'tertiary_link' OR {0}= 'trunk' OR {0}= 'trunk_link' OR {0}= 'unclassified'".format("fclass") #Select only larger roads aside from residential streets which are very often just dirt roads and simply correlate to settlement
arcpy.SelectLayerByAttribute_management('roadslyr', 'NEW_SELECTION', expr)
arcpy.Merge_management(railwaylist,os.path.join(osmdir,"railways.shp"))
arcpy.MakeFeatureLayer_management(os.path.join(osmdir,"railways.shp"), "raillyr")
expr = "{0}='construction' OR {0}='disused' OR {0}= 'funicular' OR {0}= 'light_rail' OR {0}= 'miniature' OR {0}= 'monorail' OR {0}= 'narrow_gauge' OR {0}= 'preserved' OR \
{0}= 'rail' OR {0}= 'subway' OR {0}= 'tram' OR {0}= 'unclassified'".format("fclass") #Select only actual railways rather than features (e.g. stations, stops, etc.)
arcpy.SelectLayerByAttribute_management('raillyr', 'NEW_SELECTION', expr)
arcpy.Merge_management(['roadslyr', "raillyr"], wd+'roadsrailways.shp')
roadcatchinters= wd+"roadcatchinters.shp"
arcpy.Intersect_analysis([wd+'roadsrailways.shp', scatchpoly],roadcatchinters)
arcpy.AddGeometryAttributes_management(roadcatchinters, "LENGTH_GEODESIC", "kilometers")
arcpy.Statistics_analysis(roadcatchinters, "roadsrails", [["LENGTH_GEO", "SUM"]],case_field="GridID")
arcpy.AlterField_management("roadsrails","SUM_LENGTH_GEO","CatRoadDen")
arcpy.Delete_management(os.path.join(osmdir,"roads.shp"))
arcpy.Delete_management(os.path.join(osmdir,"railways.shp"))

###########################################################
#Protected areas
wdpa = datadir+'WDPA_Mar2018_Public\\WDPA_Mar2018_Public.gdb\\WDPA_poly_Mar2018'
arcpy.Intersect_analysis([scatchpoly,wdpa],wd+'wdpacatchinters.shp')
arcpy.AddGeometryAttributes_management(wd+'wdpacatchinters.shp', Geometry_Properties="AREA_GEODESIC", Area_Unit="SQUARE_KILOMETERS")
arcpy.Statistics_analysis(wd+'wdpacatchinters.shp', "PA", [["AREA_GEO", "SUM"]], case_field="GridID")
arcpy.AlterField_management("PA","SUM_AREA_GEO","CatPAPer")
arcpy.Delete_management(wd+'wdpacatchinters.shp')

##
arcpy.ClearEnvironment('workspace')

########################################################################################################################
##Make master table of catchment characteristics and calibrate for watershed accumulation
########################################################################################################################
arcpy.CopyRows_management(scatchpoly, 'master')
arcpy.DeleteField_management('master', ['Shape_Length','Shape_Area'])
arcpy.MakeTableView_management('master', 'masterview')
for tab in arcpy.ListTables():
    print(tab)
    if not tab=='master' or not tab=='masterdepthrock':
        try:
            if 'Value' in [f.name for f in arcpy.ListFields(tab)]:
                arcpy.AddJoin_management('masterview', 'GridID', tab, 'Value', join_type='KEEP_ALL')
            elif 'catchpoly188.GridID' in [f.name for f in arcpy.ListFields(tab)]:
                arcpy.AddJoin_management('masterview', 'GridID', tab, 'catchpoly188.GridID', join_type='KEEP_ALL')
            elif 'GRIDID' in [f.name for f in arcpy.ListFields(tab)]:
                arcpy.AddJoin_management('masterview', 'GridID', tab, 'GRIDID', join_type='KEEP_ALL')
        except:
            print('Could not join')
gdbname_ws = wd+'watershed_attri.gdb\\'
arcpy.CopyRows_management('masterview',gdbname_ws +'catchment_attributes') #CopyRows and TabletoTable after joining systematically crash.Not sure why. Do it in arcmap.
#In ArcMap, join all tables and take out useless fields (OBJECTID, etc.)
#Change <Null> values to 0 when relevant (road density, dam density, etc.)
attriformat = gdbname_ws +'catchment_attributes_format'
arcpy.CopyRows_management(gdbname_ws +'catchment_attributes', gdbname_ws +'catchment_attributes_format') #CopyRows and TabletoTable after joining systematically crash.Not syre why. Do it in arcmap.

#Multiply fields by catchment area so that can then be summed and divided by total watershed area so as to obtain area-weighted average for watershed
areaprodfields=['OBJECTID','CatArea','CatSolAvg','CatSloStd','CatSloAvg','CatPoroAvg','CatPETAvg','CatPermAvg','CatEroAvg','CatElvAvg','CatDRocAvg','CatBio01Av',
 'CatBio02Av','CatBio03Av','CatBio04Av','CatBio05Av','CatBio06Av','CatBio07Av','CatBio08Av','CatBio09Av','CatBio10Av','CatBio11Av','CatBio12Av','CatBio13Av',
 'CatBio14Av','CatBio15Av','CatBio16Av','CatBio17Av','CatBio18Av','CatBio19Av','CatAIAvg','CatPopDen','CatWatExt']
with arcpy.da.UpdateCursor(attriformat, areaprodfields) as cursor:
    for row in cursor:
        print(row[0])
        for f in range(2,len(areaprodfields)):
            if row[f]== None:
                print('pass')
            else:
                row[f]=row[f]*row[1]
                cursor.updateRow(row)
del row
del cursor

#Multiply pekel water extent indicators by the extent of water so that it is weighted by water extent rather than by area.
watareafields= ['OBJECTID','CatWatExt','CatWatSea','CatWatOcc','CatWatCha']
with arcpy.da.UpdateCursor(attriformat, watareafields) as cursor:
    for row in cursor:
        print(row[0])
        for f in range(2, len(watareafields)):
            if row[f]== None:
                print('pass')
            else:
                row[f]=row[f]*row[1]
                cursor.updateRow(row)
del row
del cursor


#Run, then delete, then run 703-752, then run "Compute watershed max statistics", then append river order>1 max statistics with river order=1, then append to final table, then copy to final dataset.


######################################################################################################################
#Route attributes to full watershed for each segment
#########################################################################################################################
#Join catchment attributes to drainage lines
arcpy.MakeFeatureLayer_management(sline,'slinelyr')
arcpy.AddJoin_management('slinelyr', 'GridID', attriformat, 'GridID', 'KEEP_ALL')
slineattri = outhydro+'streamnet118_attri'
arcpy.CopyFeatures_management('slinelyr',slineattri)
#Create feature dataset
feature_dat = "tznet"
pr = arcpy.Describe(slineattri).SpatialReference
arcpy.CreateFeatureDataset_management(outhydro, feature_dat, pr)
arcpy.CopyFeatures_management(slineattri,outhydro+"\\"+feature_dat+"\\lines")
#Create geo network
out_name = "tznet"
in_source_feature_classes = "lines SIMPLE_EDGE NO"
arcpy.CreateGeometricNetwork_management(outhydro + "\\" + feature_dat, out_name, in_source_feature_classes)
network = outhydro + feature_dat + "\\" + out_name
#Set flow direction of geometric network, a pre-requisite for tracing
arcpy.SetFlowDirection_management(network, flow_option = "WITH_DIGITIZED_DIRECTION")

###################################################################
# Create points for tracing
arcpy.FeatureVerticesToPoints_management(in_features=slineattri,out_feature_class=outhydro+'linemidpoints_attri', point_location='MID')
#Make list of variables
arcpy.MakeFeatureLayer_management(outhydro+'linemidpoints_attri', 'spointlyr',where_clause="ReaOrd>1") ####### UPDATE #########
catstatslist = [[f.name, 'SUM'] for f in arcpy.ListFields('spointlyr') if not f.name in ['CatLakInd', 'CatResInd','CatElvMax','CatElvMin','ReaOrd','FlowAcc']][12:-1] #Create list of all fields to sum
sumfields = [f.name for f in arcpy.ListFields('spointlyr') if not f.name in ['CatLakInd', 'CatResInd','CatElvMax','CatElvMin','ReaOrd','FlowAcc']][12:-1]
arcpy.env.workspace=wd+'watershed_attri.gdb'

########################################################################################################################
# Compute watershed SUM statistics
########################################################################################################################
#Create template table #Do not reiterate if re-running code and no copy has been made of data
sumout_tab = gdbname_ws+ 'Ws_Attri_Sum'
scratch_tab = gdbname_ws+ 'scratch'
sumtabfields = sumfields + ['GridID']
arcpy.MakeFeatureLayer_management(outhydro+'linemidpoints_attri', 'spointlyrreduce',where_clause="ReaOrd>1") #Create another temporary layer as numpy array cannot hold the entire table (MemoryError: cannot allocate array memory)
array_sum =arcpy.da.FeatureClassToNumPyArray('spointlyrreduce', sumtabfields)
arcpy.da.NumPyArrayToTable(array_sum, scratch_tab)
#arcpy.CreateTable_management(wd+'watershed_attri.gdb', 'Ws_Attri_Sum', scratch_tab)
arcpy.Delete_management(scratch_tab)
arcpy.Delete_management('spointlyrreduce')

#Create template table for appending data when cannot allocate memory to numpy array
memerror_tab = 'Ws_Attri_Sum_MemoryError'
#arcpy.Statistics_analysis('spointlyr', memerror_tab, catstatslist)
#arcpy.AddField_management(memerror_tab, field_name='GridID', field_type='LONG')
#arcpy.CalculateField_management(memerror_tab,field='GridID', expression='999999')

nms= array_sum.dtype.names[:-1]
array_types = [(i,array_sum[i].dtype.kind) for i in nms if array_sum[i].dtype.kind in ('i', 'f')] + [('GridID','i')]

#Make subset of reaches to get data for
tz_bd = wd+'gadm1_lev0_Tanzania.shp'
arcpy.MakeFeatureLayer_management(outhydro+'linemidpoints_attri', 'spointlyr',where_clause="ReaOrd>1") ####### UPDATE #########
arcpy.SelectLayerByLocation_management('spointlyr', 'INTERSECT', tz_bd, selection_type='SUBSET_SELECTION')
arcpy.GetCount_management('spointlyr')
start = time.time()
x=0
by_col = numpy.empty([0, len(array_types)], dtype=array_types)
with arcpy.da.SearchCursor('spointlyr', ["GridID"]) as cursor:
    for row in cursor:
            #print(x)
            #x=x+1
            print(row[0])
            if row[0]==294478:
                #print(row[0])
                subcatch_id = row[0]
                expr = '"GridID" = %s' %subcatch_id
                arcpy.SelectLayerByAttribute_management('spointlyr', 'NEW_SELECTION', expr)
                arcpy.TraceGeometricNetwork_management(in_geometric_network= network, out_network_layer = "up_trace_lyr", in_flags = "spointlyr",
                                                       in_trace_task_type= "TRACE_UPSTREAM")
                up_tr = "up_trace_lyr\\lines"
                #print('Trace done!')
                try:
                    array=arcpy.da.FeatureClassToNumPyArray(up_tr, sumfields)
                    #print('Feature class to numpy array done!')
                    if len(by_col)==0:
                        by_col = numpy.array([tuple([array[i].sum() for i in nms if array[i].dtype.kind in ('i', 'f')]+[row[0]])],
                                              dtype=array_types)
                    elif len(by_col)>0 and len(by_col)<500:
                        by_col = numpy.concatenate((by_col,numpy.array([tuple([array[i].sum() for i in nms if array[i].dtype.kind in ('i', 'f')]+[row[0]])],dtype=array_types)),axis=0)
                    elif len(by_col)==500:
                        by_col = numpy.concatenate((by_col,numpy.array([tuple([array[i].sum() for i in nms if array[i].dtype.kind in ('i', 'f')]+[row[0]])],dtype=array_types)),axis=0)
                        scratch=arcpy.da.NumPyArrayToTable(by_col, scratch_tab)
                        #print('Numpy array to table done!')
                        arcpy.Append_management(scratch_tab, sumout_tab, "TEST")
                        #print('Append done!')
                        arcpy.Delete_management(scratch_tab)
                        by_col = numpy.empty([0, len(array_types)], dtype=array_types)
                except:
                    print('Error with numpy array method, use arcpy.Statistics_analysis')
                    inmemtab = r"in_memory/scratch"
                    arcpy.Statistics_analysis(up_tr, inmemtab, catstatslist)
                    print("stats analysis failed")
                    arcpy.AddField_management(inmemtab, field_name='GridID', field_type='LONG')
                    #with arcpy.da.UpdateCursor(inmemtab, ['GridID']) as cursormem:
                    #    for rowmem in cursormem:
                    #        rowmem[0]=row[0]
                    #        cursormem.updateRow(rowmem)
                    #del rowmem
                    #del cursormem
                    arcpy.CalculateField_management(inmemtab,field='GridID', expression=row[0])
                    arcpy.Append_management(inmemtab, memerror_tab, "TEST")
                    arcpy.Delete_management(inmemtab)
    scratch=arcpy.da.NumPyArrayToTable(by_col, scratch_tab)
    #print('Numpy array to table done!')
    arcpy.Append_management(scratch_tab, sumout_tab, "TEST")
arcpy.Delete_management(scratch_tab)
arcpy.Delete_management('spointlyr')
del row
del cursor
end = time.time()
print(end - start)

########################################################################################################################
# Compute watershed MAX statistics
#   Edit 2018/06/29: did not process these statistics for all of Tanzania. 10-day processing times. Future iterations could
#   merge this loop with the SUM statistics loop. In the end, compute watershed lake and reservoir index with quicker method
#   see L XXX and this loop would only be useful to compute WsElvMax, which is not used later on. So commented this section
#   out.
########################################################################################################################

#Create template table #Do not reiterate if re-running code and no copy has been made of data
# maxout_tab = gdbname_ws+ 'Ws_Attri_Max'
# scratch_tab = gdbname_ws+ 'scratch'
# maxfields = ['CatLakInd', 'CatResInd','CatElvMax'] #Create list of all fields to find max
# maxtabfields = maxfields + ['GridID']
# arcpy.MakeFeatureLayer_management(outhydro+'linemidpoints_attri', 'spointlyrreduce',where_clause="ReaOrd>3") #Create another temporary layer as numpy array cannot hold the entire table (MemoryError: cannot allocate array memory)
# array_max =arcpy.da.FeatureClassToNumPyArray('spointlyrreduce', maxtabfields)
# #arcpy.da.NumPyArrayToTable(array_max, scratch_tab)
# #arcpy.CreateTable_management(wd+'watershed_attri.gdb', 'Ws_Attri_Max', scratch_tab)
# arcpy.Delete_management(scratch_tab)
# arcpy.Delete_management('spointlyrreduce')
#
# nms= array_max.dtype.names[:-1]
# array_types = [(i,array_max[i].dtype.kind) for i in nms if array_max[i].dtype.kind in ('i', 'f')] + [('GridID','i')]
# arcpy.MakeFeatureLayer_management(outhydro+'linemidpoints_attri', 'spointlyr',where_clause="ReaOrd>1") ####### UPDATE #########
# arcpy.GetCount_management('spointlyr')
# start = time.time()
# x=0
# by_col = numpy.empty([0, len(array_types)], dtype=array_types)
# with arcpy.da.SearchCursor('spointlyr', ["GridID"]) as cursor:
#     for row in cursor:
#             print(x)
#             x=x+1
#             #if row[0]==450238:
#             #print(row[0])
#             subcatch_id = row[0]
#             expr = '"GridID" = %s' %subcatch_id
#             arcpy.SelectLayerByAttribute_management('spointlyr', 'NEW_SELECTION', expr)
#             arcpy.TraceGeometricNetwork_management(in_geometric_network= network, out_network_layer = "up_trace_lyr", in_flags = "spointlyr",
#                                                    in_trace_task_type= "TRACE_UPSTREAM")
#             up_tr = "up_trace_lyr\\lines"
#             #print('Trace done!')
#             array=arcpy.da.FeatureClassToNumPyArray(up_tr, maxfields)
#             #print('Feature class to numpy array done!')
#             if len(by_col)==0:
#                 by_col = numpy.array([tuple([array[i].max() for i in nms if array[i].dtype.kind in ('i', 'f')]+[row[0]])],
#                                       dtype=array_types)
#             elif len(by_col)>0 and len(by_col)<500:
#                 by_col = numpy.concatenate((by_col,numpy.array([tuple([array[i].max() for i in nms if array[i].dtype.kind in ('i', 'f')]+[row[0]])],dtype=array_types)),axis=0)
#             elif len(by_col)==500:
#                 by_col = numpy.concatenate((by_col,numpy.array([tuple([array[i].max() for i in nms if array[i].dtype.kind in ('i', 'f')]+[row[0]])],dtype=array_types)),axis=0)
#                 scratch=arcpy.da.NumPyArrayToTable(by_col, scratch_tab)
#                 #print('Numpy array to table done!')
#                 arcpy.Append_management(scratch_tab, maxout_tab, "TEST")
#                 #print('Append done!')
#                 arcpy.Delete_management(scratch_tab)
#                 by_col = numpy.empty([0, len(array_types)], dtype=array_types)
#     scratch=arcpy.da.NumPyArrayToTable(by_col, scratch_tab)
#     #print('Numpy array to table done!')
#     arcpy.Append_management(scratch_tab, maxout_tab, "TEST")
# arcpy.Delete_management(scratch_tab)
# arcpy.Delete_management('spointlyr')
# del row
# del cursor
# end = time.time()
# print(end - start)

########################################################################################################################
# Compute watershed lake and reservoir indices
#Order from low to high accumulation
#Loop through each segment with max accumulation, select all segments downstream
#Assign accumulation to all segments downstream
#Make subset of reaches to get data for
arcpy.MakeFeatureLayer_management(outhydro+'linemidpoints_attri', 'spointlyr',where_clause="ReaOrd>=1")
start = time.time()
x=0
maxaccdic = defaultdict(list)
errorlist = []
with arcpy.da.SearchCursor(wd+'catlakac_sub.dbf', ["VALUE", "MAX", "MAJORITY"]) as Outcursor:
    for Outrow in Outcursor:
            print(x)
            x=x+1
            print(Outrow[0])
            try:
                subcatch_id = Outrow[0]
                expr = '"GridID" = %s' %subcatch_id
                arcpy.SelectLayerByAttribute_management('spointlyr', 'NEW_SELECTION', expr)
                arcpy.TraceGeometricNetwork_management(in_geometric_network= network, out_network_layer = "down_trace_lyr", in_flags = "spointlyr",
                                                       in_trace_task_type= "TRACE_DOWNSTREAM")
                up_tr = "down_trace_lyr\\lines"
                with arcpy.da.SearchCursor(up_tr, ['GridID']) as Incursor:
                    for Inrow in Incursor:
                        maxaccdic[Inrow[0]].append([Outrow[1], Outrow[2]])
            except:
                print('No flag found')
                #These are catchments < 118 pixels without corresponding line segment. Here only occur in endorheic lakes.
                #and does not impact downstream routing as there are no downstream segments of endorheic point of max accumulation.
arcpy.Delete_management('spointlyr')
del Inrow
del Outrow
del Incursor
del Outcursor
end = time.time()
print(end - start)

#Only keep maximum accumulation value for each segment downstream of one or more lakes
maxaccdicsub = defaultdict(list)
for k, v in maxaccdic.iteritems(): #For each segment downstream of a lake
    for lst in v: #For each upstream lake
        try:
            if lst[0]>maxaccdicsub[k][0]: #If lake accumulation value > previous lake acc value for that segment
                maxaccdicsub[k]=lst #Replace the value in the output dic
        except: #If this segment does not have a lake acc value yet
            maxaccdicsub[k]=lst #Add the current lake acc value to the output dic

#Make table of watershed lake index
wslaktab =  gdbname_ws+'Ws_lakeindex'
arcpy.MakeTableView_management(sline, 'sline_view')
arcpy.AddJoin_management('sline_view', 'GridID', wd+gdbname_cat+'/lakeacc', 'Value')
arcpy.CopyRows_management('sline_view', wslaktab)

#######TO RUN, THEN MERGE WITH MAIN ATTRIBUTE TABLE, COMPUTE PROPER LAKE INDEX, RE-RUN DEDUCTIVE CLASSIFICATION
arcpy.AddField_management(in_table=wslaktab, field_name='MaxLakAcc', field_type='LONG')
arcpy.CalculateField_management(wslaktab, 'MaxLakAcc', expression=0)
arcpy.AddField_management(in_table=wslaktab, field_name='Hylak_id', field_type='LONG')

[f.name for f in arcpy.ListFields(wslaktab)]

with arcpy.da.UpdateCursor(wslaktab, ['GridID', 'Value','MaxLakAcc', 'MAX','Hylak_id','MAJORITY']) as cursor:
    for row in cursor:
        if row[0]==row[1]: #If segment intersects lake, get accumulation value from Catchment lakeacc table
            row[2]=row[3]
            row[4]=row[5]
        if row[0] in maxaccdicsub.keys(): #IF segment is downstream from a lake
            row[2]=maxaccdicsub[row[0]][0]  #Get accumulation from the dictionary
            row[2]=maxaccdicsub[row[0]][5]
        cursor.updateRow(row)
del row
del cursor

#################################################
#Make a copy and merge tables with SUM statistics
arcpy.Copy_management(in_data=memerror_tab,out_data='Ws_Attri_Sum_MemoryError_20180531')
arcpy.Copy_management(in_data=sumout_tab,out_data=sumout_tab+'_20180531')
#Delete dummy first line in memerror_tab
with arcpy.da.UpdateCursor(memerror_tab, ['GridID']) as cursor:
    for row in cursor:
        if row[0] == 999999:
            cursor.deleteRow()
del row, cursor
#Edit field names and types for both tables to match
for fd in arcpy.ListFields(sumout_tab):
    if 'Cat' in fd.name: arcpy.AlterField_management(sumout_tab, fd.name, new_field_name=fd.name.replace('Cat','Ws'))
    if 'Cat' in fd.aliasName: arcpy.AlterField_management(sumout_tab, fd.name.replace('Cat','Ws'), new_field_alias=fd.aliasName.replace('Cat','Ws'))
tab1_fields=[f.name for f in arcpy.ListFields(sumout_tab)]
for fd in arcpy.ListFields(memerror_tab):
    if 'Cat' in fd.name or 'SUM_' in fd.name: arcpy.AlterField_management(memerror_tab, fd.name, new_field_name=fd.name.replace('Cat','Ws').replace('SUM_',''))
    if 'Cat' in fd.aliasName or 'SUM_' in fd.aliasName: arcpy.AlterField_management(memerror_tab, fd.name.replace('Cat','Ws').replace('SUM_',''), new_field_alias=fd.name.replace('Cat','Ws').replace('SUM_',''))
try:
    arcpy.DeleteField_management(memerror_tab,'FREQUENCY')
except Exception:
    e = sys.exc_info()[1]
    print(e.args[0])
#Modify 'WsMineden' field type
arcpy.AddField_management(sumout_tab, field_name='WsMineDen2', field_type=[f.type for f in arcpy.ListFields(memerror_tab) if f.name=='WsMineDen'][0])
arcpy.CalculateField_management(sumout_tab, field='WsMineDen2', expression='!WsMineDen!',expression_type="PYTHON")
arcpy.DeleteField_management(sumout_tab, 'WsMineDen')
arcpy.AlterField_management(sumout_tab, field='WsMineDen2', new_field_name='WsMineDen',new_field_alias='WsMineDen')
#Modify 'WsFlowAcc' field type
arcpy.AddField_management(sumout_tab, field_name='WsFlowAcc2', field_type=[f.type for f in arcpy.ListFields(memerror_tab) if f.name=='WsFlowAcc'][0] )
arcpy.CalculateField_management(sumout_tab, field='WsFlowAcc2', expression='!WsFlowAcc!',expression_type="PYTHON")
arcpy.DeleteField_management(sumout_tab, 'WsFlowAcc')
arcpy.AlterField_management(sumout_tab, field='WsFlowAcc2', new_field_name='WsFlowAcc',new_field_alias='WsFlowAcc')
#Check for differences among tables
comparetabs=arcpy.TableCompare_management(sumout_tab, memerror_tab, sort_field='GridID', compare_type='SCHEMA_ONLY',continue_compare='CONTINUE_COMPARE',
                                          out_compare_file=wd+'compare')
#Append tables
arcpy.Append_management(memerror_tab, sumout_tab, 'TEST')
arcpy.Delete_management(wd+'compare')
#Check for duplicates and keep the one with the maximum watershed area
arcpy.FindIdentical_management(sumout_tab,sumout_tab+'_identical',fields='GridID',output_record_option='ONLY_DUPLICATES')
dupliID = [id[0] for id in arcpy.da.SearchCursor(sumout_tab+'_identical', ['IN_FID'])]
duplicates = [[row[0],row[1],row[2]] for row in arcpy.da.SearchCursor(sumout_tab, ['OBJECTID','GridID','WsArea']) if row[0] in dupliID]
if len(duplicates)>0:
    d={}
    ldel=[]
    for sub in duplicates: #Inspired from https://stackoverflow.com/questions/34334381/removing-duplicates-from-a-list-of-lists-based-on-a-comparison-of-an-element-of
        k=sub[1]
        if k not in d:
            d[k]=sub
        elif sub[2] > d[k][2]:
            ldel.append(d[k][0])
            d[k]=sub
        else:
            #print(sub[0])
            ldel.append(sub[0])
expr= 'NOT "OBJECTID" IN ' + str(tuple(ldel))
arcpy.MakeTableView_management(sumout_tab, out_view='Ws_Attri_Sum_view', where_clause=expr)
arcpy.CopyRows_management('Ws_Attri_Sum_view','Ws_Attri_Sum_nodupli')
arcpy.Delete_management(sumout_tab+'_identical')
arcpy.Delete_management(wd+'compare')

###################################################
# Join SUM and MAX statistics table
# for fd in arcpy.ListFields(maxout_tab):
#     if 'Cat' in fd.name: arcpy.AlterField_management(maxout_tab, fd.name, new_field_name=fd.name.replace('Cat','Ws'))
#     if 'Cat' in fd.aliasName: arcpy.AlterField_management(maxout_tab, fd.name.replace('Cat','Ws'), new_field_alias=fd.aliasName.replace('Cat','Ws'))
# arcpy.MakeTableView_management('Ws_Attri_Sum_nodupli', 'Ws_Attri_Sum_noduplilyr')
# arcpy.AddJoin_management('Ws_Attri_Sum_noduplilyr', 'GridID', maxout_tab, 'GridID')
# ws_tab=gdbname_ws+'Ws_Attri'
# arcpy.CopyRows_management('Ws_Attri_Sum_noduplilyr',ws_tab)
# arcpy.Delete_management('Ws_Attri_Sum_noduplilyr')
# arcpy.DeleteField_management(ws_tab, ['OBJECTID_1','GridID_1'])

###################################################
# Compute final variables for watersheds
areadivlist = [f.name for f in arcpy.ListFields(ws_tab) if not f.name in ['OBJECTID','WsWatOcc','WsWatcha','WsWatSea','WsLen_1','GridID',
                                                                              'WsRoadDen','WsDamDen','WsMineDen','WsLakInd', 'WsResInd','WsElvMax']] #Keep area as first field and don't divide for dam, roads, mines, and PA density
watextdivlist = ['WsWatExt','WsWatOcc','WsWatcha','WsWatSea']
with arcpy.da.UpdateCursor(ws_tab,watextdivlist) as cursor:
    for row in cursor:
        if row[0]>0:
            for i in range(1,len(watextdivlist)): #Iterate over every selected field after area and divide it by water area
                row[i]=float(row[i])/float(row[0])
        cursor.updateRow(row)
del row
del cursor

with arcpy.da.UpdateCursor(ws_tab,areadivlist) as cursor:
    for row in cursor:
        for i in range(1,len(areadivlist)): #Iterate over every selected field after area and divide it by watershed area
            try:
                row[i]=float(row[i])/float(row[0])
            except:
                print('Error with field: '+areadivlist[i])
        cursor.updateRow(row)
del row
del cursor

#For land cover, direction, geology,and soil, compute percentage (for each row, dividing by sum over all fields) and identify field which is majority
for var in ['Dir','Geol','LC','Soil']:
    print(var)
    if 'Ws{}Maj'.format(var) not in [f.name for f in arcpy.ListFields(ws_tab)]:
        arcpy.AddField_management(in_table=ws_tab,field_name='Ws{}Maj'.format(var),field_type='TEXT')
    varfields = [f.name for f in arcpy.ListFields(ws_tab, '{}Sum*'.format(var))]
    varfields.insert(0, 'Ws{}Maj'.format(var))
    with arcpy.da.UpdateCursor(ws_tab, varfields) as cursor:
        for row in cursor:
            if row[1] is not None: denom=sum(row[1:]) #Sometimes geology was not computed for some watershed because of small boundaries of geology data
            for i in range(1,len(varfields)):
                if row[i] == max(row):
                    row[0]=str(varfields[i][len(var)+4:]) #Find the category ID that covers the most area in the watershed
                if row[i] is not None:
                    row[i]=row[i]/denom #Compute percentage area for each category
            cursor.updateRow(row)
    del cursor, row
#For forest loss, compute percentage and get rid of WsFLosSum_0 field
with arcpy.da.UpdateCursor(ws_tab, ['WsFLosSum_0','WsFLosSum_1']) as cursor:
    for row in cursor:
        if row[1] is not None:
            row[1]=row[1]/sum(row)
            cursor.updateRow(row)
    del cursor, row
arcpy.DeleteField_management(ws_tab,'WsFLosSum_0')

########################################################################################################################
# Get remaining reach characteristics
########################################################################################################################
gdbname_reach = "reach_attri.gdb"
if not arcpy.Exists(wd+gdbname_reach):
    arcpy.CreateFileGDB_management(wd,gdbname_reach)
env.workspace = wd + gdbname_reach

#Elevation
ZonalStatisticsAsTable(sseg, "Value", wd+"srtmclean90",out_table="elv",ignore_nodata="DATA", statistics_type="MIN_MAX_MEAN")
arcpy.AlterField_management("elv","MEAN","ReaElvAvg")
arcpy.AlterField_management("elv","MIN","ReaElvMin")
arcpy.AlterField_management("elv","MAX","ReaElvMax")
#Direction
ZonalStatisticsAsTable(sseg, "Value", dirtz,out_table="dir",ignore_nodata="DATA", statistics_type="MAJORITY")
arcpy.AlterField_management("dir","MAJORITY","ReaDirMaj")
#Curvature
ZonalStatisticsAsTable(sseg, "Value", wd+"profilcurvat",out_table="profilcurvat",ignore_nodata="DATA", statistics_type="MEAN")
arcpy.AlterField_management("profilcurvat","MEAN","ReaCurvAvg")
#Protected areas
wdpa = datadir+'WDPA_Mar2018_Public\\WDPA_Mar2018_Public.gdb\\WDPA_poly_Mar2018'
arcpy.Dissolve_management(wdpa, wd+'wdpadissolve.shp')
arcpy.Intersect_analysis([sline,wdpa],wd+'wdpalineinters.shp')
arcpy.AddGeometryAttributes_management(wd+'wdpalineinters.shp', Geometry_Properties="LENGTH_GEODESIC", Length_Unit="KILOMETERS")
arcpy.Statistics_analysis(wd+'wdpalineinters.shp', "PA", [["LENGTH_GEO", "SUM"]], case_field="GridID")
arcpy.AlterField_management("PA","SUM_LENGTH_GEO",new_field_name="ReaPAPer",new_field_alias='ReaPAPer')
arcpy.Delete_management(wd+'wdpalineinters.shp')

arcpy.MakeTableView_management("elv","elvview")
arcpy.AddJoin_management("elvview", 'Value', "dir", 'Value')
arcpy.AddJoin_management("elvview", 'Value', "profilcurvat", 'Value')
arcpy.AddJoin_management("elvview", 'Value', "PA", 'GridID')
arcpy.CopyRows_management('elvview','reach_attri') #Crashes so do it in arcmap + remove redundant fields

with arcpy.da.UpdateCursor('reach_attri', ['ReaPAPer']) as cursor:
    for row in cursor:
        if row[0] is None:
            row[0]=0
            cursor.updateRow(row)
    del cursor, row

arcpy.ClearEnvironment('workspace')
########################################################################################################################
# Generate final dataset
########################################################################################################################
finalgdb = wd+'TZ_final.gdb\\'
if not arcpy.Exists(finalgdb):
    arcpy.CreateFileGDB_management(wd, 'TZ_final.gdb')
arcpy.MakeFeatureLayer_management(sline, 'slineattrilyr')
#Join catchment attributes
arcpy.AddJoin_management('slineattrilyr', 'GridID', gdbname_ws +'catchment_attributes','GridID')
#Add profile curvature
arcpy.AddJoin_management('slineattrilyr', 'GridID', wd + gdbname_cat +'\\profilcurvat','Value')
#Create subnetwork for Rufiji
slinecat = finalgdb +'streamnet118_cat'
arcpy.CopyFeatures_management('slineattrilyr', slinecat) #Sometimes crashes with join.
arcpy.Delete_management('slineattrilyr')
arcpy.DeleteField_management(slinecat, ['OBJECTID_1','HydroID_1','GridID_1','OBJECTID_12','Value','COUNT','AREA'])

#For land cover, direction, geology,and soil, compute percentage (for each row, dividing by sum over all fields) and identify field which is majority
for var in ['Dir','Geol','LC','Soil']:
    arcpy.AddField_management(in_table=slinecat,field_name='Cat{}Maj'.format(var),field_type='TEXT')
    varfields = [f.name for f in arcpy.ListFields(slinecat, '{}Sum*'.format(var))]
    varfields.insert(0, 'Cat{}Maj'.format(var))
    with arcpy.da.UpdateCursor(slinecat, varfields) as cursor:
        for row in cursor:
            if row[1] is not None: denom=sum(row[1:]) #Sometimes geology was not computed for some watershed because of small boundaries of geology data
            for i in range(1,len(varfields)):
                if row[i] == max(row):
                    row[0]=str(varfields[i][len(var)+4:]) #Find the category ID that covers the most area in the watershed
                if row[i] is not None:
                    row[i]=row[i]/denom #Compute percentage area for each category
            cursor.updateRow(row)
    del row, cursor
#For forest loss, compute percentage and get rid of WsFLosSum_0 field
with arcpy.da.UpdateCursor(slinecat, ['CatFLosSum_0','CatFLosSum_1']) as cursor:
    for row in cursor:
        if row[1] is not None:
            row[1]=row[1]/sum(row)
        cursor.updateRow(row)
arcpy.DeleteField_management(slinecat, 'CatFLosSum_0')

##########################################################################################
# Fill in watershed attributes for streams of first order (i.e. append catchment attributes)
arcpy.MakeFeatureLayer_management(slinecat, 'slinecatlyr')
arcpy.SelectLayerByAttribute_management('slinecatlyr', 'NEW_SELECTION', 'ReaOrd=1')
order1=finalgdb +'catattri_order1'
arcpy.CopyRows_management('slinecatlyr',order1)

#Format catchment attribute table into Ws table
appfields=[f.name for f in arcpy.ListFields(order1) if f.name not in ['Shape_Length','CatFlowAcc','CatElvMin','CatCurvAvg']][9:]
appfields.insert(0,'OBJECTID')
appfields.append('GridID')
len(appfields)
len([f.name for f in arcpy.ListFields(ws_tab)])
arcpy.DeleteField_management(order1, [f.name for f in arcpy.ListFields(order1) if f.name not in appfields])
for fd in arcpy.ListFields(order1): #Edit field names and types for both tables to match
    if 'Cat' in fd.name: arcpy.AlterField_management(order1, fd.name, new_field_name=fd.name.replace('Cat','Ws'))
    if 'Cat' in fd.aliasName: arcpy.AlterField_management(order1, fd.name.replace('Cat','Ws'), new_field_alias=fd.aliasName.replace('Cat','Ws'))

arcpy.Merge_management(inputs=[ws_tab,order1],output=finalgdb+'wsattri_all')
arcpy.Delete_management(order1)

#Join watershed attributes
arcpy.MakeFeatureLayer_management(slinecat, 'slinecatlyr')
arcpy.AddJoin_management('slinecatlyr', 'GridID', finalgdb+'wsattri_all', 'GridID')
slinecatws = finalgdb +'streamnet118_catws'
arcpy.CopyFeatures_management('slinecatlyr', slinecatws)
arcpy.DeleteField_management(slinecatws, ['OBJECTID_1','GridID_1'])
[f.name for f in arcpy.ListFields(slinecatws)]

#Compute mines, dams, road, PA, and drainage densities for catchment and watershed
arcpy.AddField_management(slinecatws, 'CatDen', field_type='DOUBLE')
arcpy.AddField_management(slinecatws, 'WsDen', field_type='DOUBLE')
with arcpy.da.UpdateCursor(slinecatws, ['CatArea','CatLen_1','CatDen','CatRoadDen','CatDamDen','CatMineDen','CatPAPer']) as cursor:
    for row in cursor:
        row[2]=row[1]/row[0]
        row[3]=row[3]/row[0]
        row[4]=row[4]/row[0]
        row[5]=row[5]/row[0]
        cursor.updateRow(row)
    del row, cursor
with arcpy.da.UpdateCursor(slinecatws, ['WsArea','WsLen_1','WsDen','WsRoadDen','WsDamDen','WsMineDen','WsPAPer']) as cursor:
    for row in cursor:
        row[2]=row[1]/row[0]
        row[3]=row[3]/row[0]
        row[4]=row[4]/row[0]
        row[5]=row[5]/row[0]
        cursor.updateRow(row)
    del row, cursor

#Compute lake index and reservoir index
with arcpy.da.UpdateCursor(slinecatws, ['CatFlowAcc','WsLakInd', 'WsResInd']) as cursor:
    for row in cursor:
        row[1]=row[1]/row[0]
        row[2]=row[2]/row[0]
        cursor.updateRow(row)
    del row, cursor

###############################
# Add reach characteristics
arcpy.MakeFeatureLayer_management(slinecatws, 'slinecatws_lyr')
arcpy.AddJoin_management('slinecatws_lyr', 'GridID', wd + gdbname_reach + '\\reach_attri','Value_')
slinefinal = finalgdb +'streamnet118_final'
arcpy.CopyFeatures_management('slinecatws_lyr', slinefinal)
arcpy.DeleteField_management(slinefinal, ['OBJECTID_1','Value_'])

#Compute reach slope
arcpy.AddField_management(slinefinal, field_name='ReaSloAvg',field_type='FLOAT')
arcpy.CalculateField_management(slinefinal, 'ReaSloAvg', "math.degrees(math.atan((!ReaElvMax!-!ReaElvMin!)/(1000*!CatLen!)))", expression_type='PYTHON')

################################EDIT FIELD NAMES ########################################
[f.name for f in arcpy.ListFields(slinefinal)]
#Rename geology fields
for f in arcpy.ListFields(slinefinal):
    if 'GeolSum' in f.name:
        print(f.name)
        if len(f.name[8:])<3:
            if float(f.name[8:])<67:
                arcpy.AlterField_management(slinefinal, f.name, 'CatGeol'+f.name[8:], new_field_alias='CatGeol'+f.name[8:])
            else:
                arcpy.AlterField_management(slinefinal, f.name, 'WsGeol'+f.name[8:], new_field_alias='WsGeol'+f.name[8:])
        else:
            arcpy.AlterField_management(slinefinal, f.name, 'WsGeol'+f.name[8:10], new_field_alias='WsGeol'+f.name[8:10])
#Rename soil fields
flist= [f.name for f in arcpy.ListFields(slinefinal)]
split_index = flist.index('CatSoilMaj')
for i in range(0,len(flist)):
    if 'SoilSum' in flist[i]:
        if i < split_index:
            arcpy.AlterField_management(slinefinal, flist[i], 'CatSoil'+flist[i][8:], new_field_alias='CatSoil'+flist[i][8:])
        if i > split_index:
            arcpy.AlterField_management(slinefinal, flist[i], 'WsSoil'+flist[i][8:], new_field_alias='WsSoil'+flist[i][8:])

#Rename LCSum fields
arcpy.AlterField_management(slinefinal, 'LCSum_1', 'CatTreePer', new_field_alias='CatTreePer')
arcpy.AlterField_management(slinefinal, 'LCSum_2', 'CatShruPer', new_field_alias='CatShruPer')
arcpy.AlterField_management(slinefinal, 'LCSum_3', 'CatGrasPer', new_field_alias='CatGrasPer')
arcpy.AlterField_management(slinefinal, 'LCSum_4', 'CatCropPer', new_field_alias='CatCropPer')
arcpy.AlterField_management(slinefinal, 'LCSum_5', 'CatAqVegPer', new_field_alias='CatAqVegPer')
arcpy.AlterField_management(slinefinal, 'LCSum_6', 'CatSparPer', new_field_alias='CatSparPer')
arcpy.AlterField_management(slinefinal, 'LCSum_7', 'CatBarePer', new_field_alias='CatBarePer')
arcpy.AlterField_management(slinefinal, 'LCSum_8', 'CatUrbPer', new_field_alias='CatUrbPer')
arcpy.AlterField_management(slinefinal, 'LCSum_10', 'CatWatPer', new_field_alias='CatWatPer')

arcpy.AlterField_management(slinefinal, 'LCSum_12', 'WsTreePer', new_field_alias='WsTreePer')
arcpy.AlterField_management(slinefinal, 'LCSum_23', 'WsShruPer', new_field_alias='WsShruPer')
arcpy.AlterField_management(slinefinal, 'LCSum_34', 'WsGrasPer', new_field_alias='WsGrasPer')
arcpy.AlterField_management(slinefinal, 'LCSum_45', 'WsCropPer', new_field_alias='WsCropPer')
arcpy.AlterField_management(slinefinal, 'LCSum_56', 'WsAqVegPer', new_field_alias='WsAqVegPer')
arcpy.AlterField_management(slinefinal, 'LCSum_67', 'WsSparPer', new_field_alias='WsSparPer')
arcpy.AlterField_management(slinefinal, 'LCSum_78', 'WsBarePer', new_field_alias='WsBarePer')
arcpy.AlterField_management(slinefinal, 'LCSum_89', 'WsUrbPer', new_field_alias='WsUrbPer')
arcpy.AlterField_management(slinefinal, 'LCSum_10_11', 'WsWatPer', new_field_alias='WsWatPer')

#Rename/Delete extra length field
arcpy.DeleteField_management(slinefinal, 'CatLen_1')
arcpy.AlterField_management(slinefinal, 'WsLen_1', 'WsLen', new_field_alias='WsLen')

#Rename forest loss field
arcpy.AlterField_management(slinefinal, 'CatFLosSum_1', 'CatFLosPer', new_field_alias='CatFLosPer')
arcpy.AlterField_management(slinefinal, 'WsFLosSum_1', 'WsFLosPer', new_field_alias='WsFLosPer')

#Export environmental attribute table to dbf
arcpy.CopyRows_management(slinefinal, wd+'streamnet118_final.csv')