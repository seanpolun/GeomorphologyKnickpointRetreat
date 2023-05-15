import arcpy
import csv
import math
import pandas as pd
import numpy as np
import os
import sys
import scipy.io as spio

## ISSUES
#
# 6-13-18
# Geometries are sometimes read out of sequence. Not entire sure how to fix. -sgp
# Fixed: Make sure geometries are singlepart: Use Multipart to Singlepart, then integrate to remove
# small gaps, then unsplit lines. - sgp
#

waterShed = arcpy.GetParameterAsText(0)
inputLine = arcpy.GetParameterAsText(1)
inputFac = arcpy.GetParameterAsText(2)
inputDEM = arcpy.GetParameterAsText(3)
outputDirectory = arcpy.GetParameterAsText(4)
outputName = arcpy.GetParameterAsText(5)

arcpy.env.workspace = "in_memory"

WshedTable = arcpy.da.FeatureClassToNumPyArray(waterShed,("Shape_Area","HydroID","Name"))
WShedID = WshedTable["HydroID"]
WShedNames = WshedTable["Name"]

FacLines = "in_memory\\FacLines"
ElevLines = "in_memory\\ElevLines"
FacArea = 30*30



arcpy.InterpolateShape_3d(inputFac,inputLine,FacLines,z_factor=FacArea,vertices_only="VERTICES_ONLY")
arcpy.InterpolateShape_3d(inputDEM,inputLine,ElevLines,z_factor=1,vertices_only="VERTICES_ONLY")
FacCursor = arcpy.da.SearchCursor(FacLines,["SHAPE@"])
Basins = []
outputData = {}
with arcpy.da.SearchCursor(ElevLines,["OID@","SHAPE@","DrainID"]) as LineSearch:
    for row in LineSearch:
#        try:
        FacCursor.next()
        geoElev = row[1]
        geoFac = FacCursor[0]
        #TotalLeng = row[2]
        BasinID = row[2]
        Basins.append(BasinID)
        idd = np.where(WShedID == BasinID)[0][0]
        BasinName = WShedNames[idd]

        start = geoFac.firstPoint
        x2 = start.X
        y2 = start.Y
        dist2 = 0
        part_z = []
        partDist = []
        for part in geoElev.getPart():
            for pnt in part:
                part_z.append(pnt.Z)
                x1 = pnt.X
                y1 = pnt.Y
                dist1 = math.sqrt((x1 - x2)**2 + (y1 - y2)**2) + dist2
                partDist.append(dist1)
                x2 = x1
                y2 = y1
                dist2 = dist1
        Elev = np.asarray(part_z)
        Dist = np.asarray(partDist)
        partFac = []
        for part in geoFac.getPart():
            for pnt in part:
                partFac.append(pnt.Z)
        Fac = np.asarray(partFac)
        outDict = {"MDist":Dist,"Elev":Elev,"FlowArea":Fac}
        outputData[BasinName] = outDict

#        except:
#            arcpy.AddWarning("an error occurred")
#            continue


matout_name = outputName + ".mat"
mat_file = outputDirectory + "\\" + matout_name
spio.savemat(mat_file,outputData)
