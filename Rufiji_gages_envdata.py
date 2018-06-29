__author__ = 'Mathis Messager'
#Creation date: 2018/03/08
#Last updated: 2018/26/08

#Project: USAID/TANZANIA
#Technical Assistance to Support the Development of Irrigation and Rural Roads Infrastructure Project (IRRIP2)
#(Contract No. EDH-I-00-08-00023-00 621-TO 12)

#Purpose: Link gages to GIS river network and export environmental characteristics
import arcpy
from arcpy import env
from arcpy.sa import *
import os
import sys

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
finalgdb = wd+'rufiji_final.gdb\\'
env.workspace = wd
#Usual variables
gages= datadir+'sharepoint20180316\gis\hydro\\rufiji\RRB_RBWB_RGS_edit.shp'
slinerufifinal = finalgdb +'streamnet118_rufiji_final'

#################################
#Snap gages to river network
#################################
gagegdb = wd+'rufiji_gages.gdb'
if not arcpy.Exists(gagegdb):
        arcpy.CreateFileGDB_management(wd, 'rufiji_gages.gdb')
env.workspace=gagegdb
gagesnap = 'gages'
arcpy.CopyFeatures_management(gages, gagesnap)
[f.name for f in arcpy.ListFields(gagesnap)]

#Project both river network and gages
linesproj = 'streamnet118_rufiji_final_proj'
gagesproj = 'gagesnapproj'
sr_afeq = 102023 #Project to Africa Equidistant Conic
arcpy.Project_management(slinerufifinal, linesproj, sr_afeq)
arcpy.Project_management(gagesnap, gagesproj, sr_afeq)

#Compute coordinates
arcpy.AddGeometryAttributes_management(gagesproj, 'POINT_X_Y_Z_M')
#Check distance between gage and network
in_features = gagesproj
near_features = slinerufifinal
search_radius = "500 meters"
location = "LOCATION"
angle = "NO_ANGLE"
method = "GEODESIC"
arcpy.Near_analysis(in_features, near_features, search_radius, location, angle)

#Snap gages to nearest river reach within 50 meters, otherwise do not change position
snap_env = [linesproj, "EDGE", "50 meters"]
arcpy.Snap_edit(gagesproj, [snap_env])
#Inspect those gages using satellite imagery and alternative river networks,verify that those snapped make sense, manually snap the others — in ArcGIS
#40/58 are snapped manually
#arcpy.CopyFeatures_management(gagesproj, "snap_edit")
#arcpy.AddField_management("snap_edit", 'snap_manual', 'TEXT')
#arcpy.AddField_management("snap_edit", 'Comment', 'TEXT',field_length=200)

#####################################################################
# SNAP FROM HERE IF RIVER ENVIRONMENTAL ATTRIBUTES HAVE BEEN CHANGED
#Join to river network by location
gages_netjoin = 'gages_netjoin'
arcpy.SpatialJoin_analysis("snap_edit",linesproj, gages_netjoin, match_option='INTERSECT')
#Export gages environmental attribute table to dbf
arcpy.CopyRows_management(gages_netjoin, wd+'gages_netjoin.csv')
#Export rest of environmental data to dbf
arcpy.CopyRows_management(slinerufifinal, wd+'streamnet118_rufiji_finaltab.csv')