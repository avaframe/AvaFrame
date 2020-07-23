#!/usr/bin/env python
# coding: utf-8
from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.PyQt.QtGui import QFont, QColor
from qgis.core import (
    QgsField,
    QgsFeature,
    QgsFeatureSink,
    QgsFeatureRequest,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSink,
    QgsProcessingFeedback,
    QgsVectorLayer,
    QgsMarkerSymbol,
    QgsGeometry,
    QgsProject,
    QgsField,
    QgsFields,
    QgsMessageLog,
    QgsPoint,
    Qgis,
    QgsSingleSymbolRenderer,
    QgsVectorFileWriter,
    QgsUnitTypes,
    QgsRasterLayer
)


import processing


import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from qgis.core import *
import qgis
from shapely.geometry import Point, LineString
import time
import datetime

import matplotlib as mpl
colors = ["#393955","#8A8A9B","#E9E940"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

class AlphaBeta(QgsProcessingAlgorithm):
    DEM = 'DEM'
    PROFILES = 'PROFILES'
    SPLITPOINTS = 'SPLITPOINTS'
    SMALLAVA = 'SMALLAVA'
    OUTPUT = 'OUTPUT'
    OUTPUT_LAYER = 'OUTPUT_LAYER'
    OUTPUTDATA = 'OUTPUTDATA'
    START = 'START'

    def __init__(self):
        super().__init__()

    def name(self):
        return "AlphaBeta"

    def tr(self, text):
        return QCoreApplication.translate("AlphaBeta", text)

    def displayName(self):
        return self.tr("AlphaBeta")

    def group(self):
        return self.tr("WLV Tools")

    def groupId(self):
        return "scripts"

    def shortHelpString(self):
        hstring = "AlphaBeta 4.1 \n\
        Textbaustein (erweiterbar) für Bericht, bitte kopieren und die nicht benötigten Punkte löschen: \n\
        Aufgrund unten genannter Kriterien wird eine Verschiebung des Alpha \
        Punktes mit X Standardabweichungen gewählt. \n\
        Punkte die für einen längeren Auslauf sprechen: \n\
        - Anbruchgebiet größer 6 ha \n\
        - Auslauf liegt über 1600hm SH \n\
        - Starke Einwehung \n\
        - Hohes Entrainment \n\
        Punkte die für einen kürzeren Auslauf sprechen: \n\
        - Anbruchgebiet kleiner 2ha \n\
        - Starke Auswehung \n\
        - Markante Richtungsänderungen \n\
        - Massenverlust über 20% \n\
        - Starke Bewaldung "

        return self.tr(hstring)

    def helpUrl(self):
        return "https://qgis.org"

    def createInstance(self):
        return type(self)()

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.DEM,
            self.tr("DEM layer")))
            # [QgsProcessing.TypeRaster]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.PROFILES,
            self.tr("Profile layer"),
            [QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.SPLITPOINTS,
            self.tr("Splitpoint layer"),
            [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterBoolean(
            self.SMALLAVA,
            self.tr("Small avalanche")
            ))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUTDATA,
            self.tr("AlphaBetaData"),
            QgsProcessing.TypeVectorAnyGeometry))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT,
            self.tr("AlphaBeta"),
            QgsProcessing.TypeVectorAnyGeometry))


    def processAlgorithm(self, parameters, context, feedback):

        sourceDEM = self.parameterAsRasterLayer(parameters, self.DEM, context)
        sourceProfiles = self.parameterAsVectorLayer(parameters, self.PROFILES, context)
        sourceSplitPoints = self.parameterAsVectorLayer(parameters,
                                                        self.SPLITPOINTS, context)

        # parameterAsBoolean is only available sind QGIS 3.8
        # smallAva = self.parameterAsBoolean(parameters, self.SMALLAVA, context)
        smallAva = self.parameterAsBool(parameters, self.SMALLAVA, context)


        QgsMessageLog.logMessage(sourceDEM.source(),
                                 'AlphaBeta', level=Qgis.Info)
        QgsMessageLog.logMessage(sourceProfiles.source(),
                                 'AlphaBeta', level=Qgis.Info)

        outPath = os.path.dirname(sourceProfiles.source())

        prepLine, SplitPoints, LineSource = self.PrepareLine(sourceProfiles,
                                                             sourceSplitPoints)

        procLine = self.ProcessLine(prepLine,sourceDEM.source(), LineSource)

        ablayer = self.AlphaBetaMain(procLine, SplitPoints, outPath, smallAva, feedback)

        dest_id=0
        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            ablayer.fields(),
            ablayer.wkbType(),
            ablayer.sourceCrs())

        features = ablayer.getFeatures(QgsFeatureRequest())
        sink.addFeatures(features, QgsFeatureSink.FastInsert)

        (sinkdata, destdata_id) = self.parameterAsSink(
            parameters,
            self.OUTPUTDATA,
            context,
            procLine.fields(),
            procLine.wkbType(),
            procLine.sourceCrs()
        )

        features = procLine.getFeatures(QgsFeatureRequest())
        sinkdata.addFeatures(features)


        # Set DataLayer visibility
        # Works, but only AFTER the script.... no solution yet...
        # root = QgsProject.instance().layerTreeRoot()
        # for child in root.children():
        #     if isinstance(child, QgsLayerTreeGroup):
        #         get_group_layers(child)
        #     elif isinstance(child, QgsLayerTreeLayer):
        #         print ('- layer: ' + child.name())
        #         if child.name() == 'AlphaBetaData' :
        #             child.setItemVisibilityChecked(False)

        #TODO
        # for feat in features:
        #     out_feat = QgsFeature()
        #     out_feat.setGeometry(feat.geometry())
        #     out_feat.setAttributes(feat.attributes())
        #     sink.addFeature(out_feat, QgsFeatureSink.FastInsert)
        # Add labels (maybe option?)

        # layer_settings  = QgsPalLayerSettings()
        # text_format = QgsTextFormat()

        # text_format.setFont(QFont("Arial", 8))
        # text_format.setSize(8)

        # buffer_settings = QgsTextBufferSettings()
        # buffer_settings.setEnabled(True)
        # buffer_settings.setSize(0.10)
        # buffer_settings.setColor(QColor("black"))

        # text_format.setBuffer(buffer_settings)
        # layer_settings.setFormat(text_format)

        # layer_settings.fieldName = "Name"
        # layer_settings.placement = 4

        # layer_settings.enabled = True

        # layer_settings = QgsVectorLayerSimpleLabeling(layer_settings)

        # sink.setLabelsEnabled(True)
        # sink.setLabeling(layer_settings)
        # sink.triggerRepaint()

        return {self.OUTPUT: dest_id}

    def AlphaBetaMain(self, layer,SplitPoints, outPath, smallAva, feedback):

        # Paramters for alpha beta
        abVersion = '4.1'
        feedback.pushInfo('Running Version: '+abVersion)

        if smallAva:
            feedback.pushInfo('KLEINLAWINEN-Setup\n\n')
            k1 = 0.933
            k2 = 0.0
            k3 = 0.0088
            k4 = -5.02
            SD = 2.36

            ParameterSet = "Kleinlawinen"
            LayerShortAppendix = "SM"

        else:
            feedback.pushInfo('Standard\n\n')
            k1 = 1.05
            k2 = -3130.0
            k3 = 0.0
            k4 = -2.38
            SD = 1.25

            ParameterSet = "Standard"
            LayerShortAppendix = "STD"



        #NOTE: Pandas would be way more convienient, but requires the
        #user to install it, as it is not included in a standard QGIS install...
        # get info about features
        # layer = QgsProject.instance().mapLayersByName(DataLayer)[0]
        features = layer.getFeatures()
        FeatDict = dict()
        field_names = [field.name() for field in layer.fields() ]

        crs = layer.crs().authid()

        ablayer = QgsVectorLayer('Point?crs=%s' % crs, 'AlphaBeta', 'memory')

        # Set the provider to accept the data source
        prov = ablayer.dataProvider()

        feats = []

        ablayer.startEditing()   # actually writes attributes

        prov.addAttributes([QgsField("fid", QVariant.Int),
                            QgsField("Name", QVariant.String),
                            QgsField("Angle", QVariant.Double),
                            QgsField("Distance", QVariant.Double)
        ])
        ablayer.commitChanges()

        # get fields for propagation to result layer
        fields = ablayer.fields()

        # fetch attributes
        # and split in featuers. i.e. in profiles
        for feature in features:
            fieldDict = dict(zip(field_names, feature.attributes()))
            FeatDict.setdefault(fieldDict['fid'],[]).append(fieldDict)

        # Sort data based on distance => c_meters
        for cukey in FeatDict.keys():
            CuFeat = sorted(FeatDict[cukey], key=lambda k: k['distance'])

            # get the matching splitpoint
            CuSplit = next(item for item in SplitPoints if item["fid"] == cukey)

            #Stackoverflow... nuff said...
            # convert to lists (dict are not necessarily sorted)
            Data = {k: [j[x] for j in CuFeat for x in j if x == k] for i in CuFeat for k in i}

            # # PANDAS Version
            # infile = "./LineTestSingle/LineTestSingle_ChRastGeomExport.csv"
            # data = pds.read_csv(infile,delimiter = ',')
            # data['dx'] = np.abs(data['xcoord'] - data['xcoord'].shift(1))
            # data['dy'] = np.abs(data['ycoord'] - data['ycoord'].shift(1))
            # data['dz'] = np.abs(data['DGM5mp'] - data['DGM5mp'].shift(1))
            # data['l'] = np.sqrt(data['dx']**2 + data['dy']**2)
            # data['Angle'] = np.rad2deg(np.arctan2(data['dz'], data['l']))
            # data.fillna(0,inplace=True)
            # data['lsum'] = data['l'].cumsum()

            # NON Pandas Version
            # id,x,z,xcoord,ycoord = np.genfromtxt(infile,skip_header=1,delimiter=',',unpack=True)
            xcoord = np.array(Data['xcoord'])
            ycoord = np.array(Data['ycoord'])
            x = np.array(Data['distance'])
            z = np.array(Data['z'])


            # Sanity check if first element of z is highest:
            # if not, flip all arrays
            if z[-1] > z[0]:
                print('Uho, reverse profile')
                z = np.flip(z)
                x = x[-1] - np.flip(x)
                xcoord = np.flip(xcoord)
                ycoord = np.flip(ycoord)

            dx = np.abs(x - np.roll(x,1))
            dz = np.abs(z - np.roll(z,1))
            dx[0] = 0.0
            dz[0] = 0.0
            angle = np.rad2deg(np.arctan2(dz, dx))


            #TODO SPLIT POINT READING
            # get all values where Angle < 10 but >0
            # get index of first occurance and go one back to get previous value
            # (i.e. last value above 10 deg)
            # tmp = x[(angle < 10.0) & (angle > 0.0) & (x > 450)]
            tmp = np.where((angle < 10.0) & (angle > 0.0) & (x > CuSplit['distance']))
            idx_10Point = tmp[0][0] - 1

            # Do a quadtratic fit and get the polynom for 2nd derivative later
            zQuad = np.polyfit(x,z,2)
            poly = np.poly1d(zQuad)

            # Get H0: max - min for parabola
            H0 = max(poly(x)) - min(poly(x))

            # get beta
            dz_beta = z[0] - z[idx_10Point]
            beta = np.rad2deg(np.arctan2(dz_beta, x[idx_10Point]))

            # get Alpha
            alpha = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4

            # get Alpha standard deviations
            SDs = [SD,-1*SD, -2*SD]
            alphaSD = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4 + SDs

            # Line down to alpha
            f = z[0] + np.tan(np.deg2rad(-alpha)) * x
            fplus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[0])) * x
            fminus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[1])) * x
            fminus2SD = z[0] + np.tan(np.deg2rad(-alphaSD[2])) * x

            # First it calculates f - g and the corresponding signs
            # using np.sign. Applying np.diff reveals all
            # the positions, where the sign changes (e.g. the lines cross).
            idx_alpha = np.argwhere(np.diff(np.sign(f - z))).flatten()
            idx_alphaP1SD = np.argwhere(np.diff(np.sign(fplus1SD - z))).flatten()
            idx_alphaM1SD = np.argwhere(np.diff(np.sign(fminus1SD - z))).flatten()
            idx_alphaM2SD = np.argwhere(np.diff(np.sign(fminus2SD - z))).flatten()

            # Only get the first index past the splitpoint
            try:
                idx_alpha = idx_alpha[x[idx_alpha]>CuSplit['distance']][0]
            except:
                feedback.pushInfo('Alpha out of profile')
                idx_alpha = None

            try:
                idx_alphaP1SD = idx_alphaP1SD[x[idx_alphaP1SD]>CuSplit['distance']][0]
            except:
                feedback.pushInfo('+1 SD above beta point')
                idx_alphaP1SD = None

            try:
                idx_alphaM1SD = idx_alphaM1SD[x[idx_alphaM1SD]>CuSplit['distance']][0]
            except:
                feedback.pushInfo('-1 SD out of profile')
                idx_alphaM1SD = None

            try:
                idx_alphaM2SD = idx_alphaM2SD[x[idx_alphaM2SD]>CuSplit['distance']][0]
            except:
                feedback.pushInfo('-2 SD out of profile')
                idx_alphaM2SD = None

            # Plot the whole shebang
            plt.close("all")
            fig = plt.figure(int(cukey),figsize=(10,6))
            titleText = 'Profil: '+CuFeat[0]['Name'] #+ \
                # '  (' + datetime.datetime.now().strftime("%d.%m.%y")  + \
                # '; ' + 'AlphaBeta ' + abVersion + ')'
            plt.title(titleText)

            xlabelText = 'Distance [m]\nBeta: '+str(round(beta,1))+ '$^\circ$' + \
                '  Alpha: '+str(round(alpha,1)) + '$^\circ$'
            plt.xlabel(xlabelText,multialignment='center')

            plt.ylabel('Height [m]')
            # plt.plot(x,z,'-', color=colors[0], label = 'DEM')
            plt.plot(x,z,'-', label = 'DEM')
            plt.plot(x,poly(x),':', label = 'QuadFit')
            plt.axvline(x=x[idx_10Point], color='0.8', \
                        linewidth=1, linestyle='-.',label='Beta')

            plt.plot(x, f, '-', label = 'AlphaLine')
            if idx_alpha:
                plt.plot(x[idx_alpha], z[idx_alpha], 'o', label='Alpha')
            if idx_alphaP1SD:
                plt.plot(x[idx_alphaP1SD], z[idx_alphaP1SD], 'd',markersize=5,
                        label='Alpha + 1SD')
            if idx_alphaM1SD:
                plt.plot(x[idx_alphaM1SD], z[idx_alphaM1SD], 'x',markersize=8,
                        label='Alpha - 1SD')
            if idx_alphaM2SD:
                plt.plot(x[idx_alphaM2SD], z[idx_alphaM2SD], 'x', markersize=8,
                        label='Alpha - 2SD')

            ax = plt.gca()
            # plt.text(0, 0, 'matplotlib', horizontalalignment='center', \
            #          verticalalignment='center', transform=ax.transAxes)
            versionText =  datetime.datetime.now().strftime("%d.%m.%y")  + \
                '; ' + 'AlphaBeta ' + abVersion + ' ' + ParameterSet
            plt.text(0, 0, versionText, fontsize=8 , verticalalignment='bottom', \
                     horizontalalignment='left', transform=ax.transAxes, \
                     color='0.5')
            # plt.text(-0.2, 0, 'matplotlib -2', \
            #          verticalalignment='center', transform=ax.transAxes)

            plt.gca().set_aspect('equal', adjustable='box')
            plt.grid(linestyle=':',color='0.9')
            plt.legend(frameon=False)
            plt.draw()

            # save_file = outPath + CuFeat[0]['Name'] + '.pdf'
            save_file = os.path.join(outPath, 'AlphaBeta_' + CuFeat[0]['Name'] + '_' + LayerShortAppendix+ '.pdf')
            plt.savefig(save_file)
            plt.close(fig)
            plt.close("all")
            # time.sleep(2)
            # plt.show()


            if idx_alpha:
                pointAlpha = QgsPointXY(xcoord[idx_alpha],ycoord[idx_alpha])
            if idx_alphaP1SD:
                pointAlphaP1SD = QgsPointXY(xcoord[idx_alphaP1SD],ycoord[idx_alphaP1SD])
            if idx_alphaM1SD:
                pointAlphaM1SD = QgsPointXY(xcoord[idx_alphaM1SD],ycoord[idx_alphaM1SD])
            if idx_alphaM2SD:
                pointAlphaM2SD = QgsPointXY(xcoord[idx_alphaM2SD],ycoord[idx_alphaM2SD])
            pointBeta = QgsPointXY(xcoord[idx_10Point],ycoord[idx_10Point])

            # Add a new feature and assign the geometry
            if idx_alpha:
                featAlpha = QgsFeature(fields)
                featAlpha['fid'] = int(cukey)
                featAlpha['Angle'] = float(round(alpha,2))
                featAlpha['Name'] =  CuFeat[0]['Name'] + '_Alpha'
                featAlpha['Distance'] = float(x[idx_alpha])
                featAlpha.setGeometry(QgsGeometry.fromPointXY(pointAlpha))
                feats.append(featAlpha)

            featBeta = QgsFeature(fields)
            featBeta['fid'] = int(cukey)
            featBeta['Angle'] = float(round(beta,2))
            featBeta['Name'] =  CuFeat[0]['Name'] + '_Beta'
            featBeta['Distance'] = float(x[idx_10Point])
            featBeta.setGeometry(QgsGeometry.fromPointXY(pointBeta))
            feats.append(featBeta)

            # Add a new feature and assign the geometry
            if idx_alphaP1SD:
                featAlphaP1SD = QgsFeature(fields)
                featAlphaP1SD['fid'] = int(cukey)
                featAlphaP1SD['Angle'] = float(round(alphaSD[0],2))
                featAlphaP1SD['Name'] =  CuFeat[0]['Name'] + '_AlPlus1SD'
                featAlphaP1SD['Distance'] = float(x[idx_alpha])
                featAlphaP1SD.setGeometry(QgsGeometry.fromPointXY(pointAlphaP1SD))
                feats.append(featAlphaP1SD)

            # Add a new feature and assign the geometry
            if idx_alphaM1SD:
                featAlphaM1SD = QgsFeature(fields)
                featAlphaM1SD['fid'] = int(cukey)
                featAlphaM1SD['Angle'] = float(round(alphaSD[1],2))
                featAlphaM1SD['Name'] =  CuFeat[0]['Name'] + '_AlMinus1SD'
                featAlphaM1SD['Distance'] = float(x[idx_alpha])
                featAlphaM1SD.setGeometry(QgsGeometry.fromPointXY(pointAlphaM1SD))
                feats.append(featAlphaM1SD)

            if idx_alphaM2SD:
                # Add a new feature and assign the geometry
                featAlphaM2SD = QgsFeature(fields)
                featAlphaM2SD['fid'] = int(cukey)
                featAlphaM2SD['Angle'] = float(round(alphaSD[2],2))
                featAlphaM2SD['Name'] =  CuFeat[0]['Name'] + '_AlMinus2SD'
                featAlphaM2SD['Distance'] = float(x[idx_alpha])
                featAlphaM2SD.setGeometry(QgsGeometry.fromPointXY(pointAlphaM2SD))
                feats.append(featAlphaM2SD)

        prov.addFeatures(feats)


        ablayer.updateExtents()
        ablayer.commitChanges()
        # Add the layer to the Layers panel
        # QgsProject.instance().addMapLayers([ablayer])

        return(ablayer)


    def create_points_at(self, distance,
                        fid,
                        name,
                        geom,
                        ):
        """Creating Points at coordinates along the line
        """

        if distance <= 0:
            distance = geom.length()

        length = geom.length()

        current_distance = distance
        startpoint = 0
        pn = 0

        feats = []

        # set the first point at startpoint
        point = geom.interpolate(startpoint)
        # convert 3D geometry to 2D geometry as OGR seems to have problems with this
        point = QgsGeometry.fromPointXY(point.asPoint())

        field_id = QgsField(name="fid", type=QVariant.Int)
        field_name = QgsField(name="Name", type=QVariant.String)
        field_pointnum = QgsField(name="PointNum", type=QVariant.Int)
        field_dist = QgsField(name="dist", type=QVariant.Double)
        fields = QgsFields()

        fields.append(field_id)
        fields.append(field_name)
        fields.append(field_pointnum)
        fields.append(field_dist)
        # fields.append(field)

        feature = QgsFeature(fields)
        feature['fid'] = fid
        feature['Name'] = name
        feature['PointNum'] = pn
        feature['dist'] = startpoint

        feature.setGeometry(point)
        feats.append(feature)

        while startpoint + current_distance <= length:
            pn = pn + 1
            # Get a point along the line at the current distance
            point = geom.interpolate(startpoint + current_distance)
            # Create a new QgsFeature and assign it the new geometry
            feature = QgsFeature(fields)
            feature['fid'] = fid
            feature['name'] = name
            feature['PointNum'] = pn
            feature['dist'] = (startpoint + current_distance)
            feature.setGeometry(point)
            feats.append(feature)
            # Increase the distance
            current_distance = current_distance + distance

        return feats


    def PrepareLine(self,layer,
                    splayer,
                    distance=10):

        # layer = QgsProject.instance().mapLayersByName(inlayer)[0]
        # layer = inlayer
        # splayer = QgsProject.instance().mapLayersByName(splitPointLayer)[0]

        crs = layer.crs().authid()

        # TODO check for virtual or shapelayer and set virt_layer according to it
        layer_type = "memory"

        layerout = "AB_LinePrep"
        virt_layer = QgsVectorLayer("Point?crs=%s" % crs,
                                    layerout,
                                    layer_type)

        provider = virt_layer.dataProvider()
        virt_layer.startEditing()   # actually writes attributes

        provider.addAttributes([QgsField("fid", QVariant.Int),
                                QgsField("Name", QVariant.String),
                                QgsField("PointNum", QVariant.Int),
                                QgsField("distance", QVariant.Double)
        ])

        def get_features():
                """Getting the features
                """
                if selected_only:
                    return layer.selectedFeatures()
                else:
                    return layer.getFeatures()

        QgsMessageLog.logMessage('Gonna Prepare da shit out of these lines',
                                 'AlphaBeta', level=Qgis.Info)

        field_names = [field.name() for field in layer.fields() ]

        # Loop through all (selected) features
        #     # Maybe use the get_feature def from above?
        for feature in layer.getFeatures():
            geom = feature.geometry()
            # Add feature ID of selected feature
            fid = feature.id()

            if not geom:
                QgsMessageLog.logMessage("No geometry", "QChainage")
                continue
            fieldDict = dict(zip(field_names, feature.attributes()))

            # TODO: flip line if lower part first

            #TODO: set generic name if none is given
            if 'Name' in fieldDict:
                AvaName = fieldDict['Name']
            else:
                AvaName = 'Profile_' + str(fid)

            QgsMessageLog.logMessage(AvaName,
                                 'AlphaBeta', level=Qgis.Info)

            features = self.create_points_at(distance,
                                        fid,
                                        AvaName,
                                        geom)
            provider.addFeatures(features)
            virt_layer.updateExtents()
            # line=geom
            # line = LineString(geom)
            # print(geom)


        #######################
        # Split Point handling
        #######################
        #Prepare SpatialIndex for SplitPoint analysis
        spIndex = QgsSpatialIndex() #create spatial index object

        feat = QgsFeature()
        fit = provider.getFeatures() #gets all features in layer

        # insert features to index
        while fit.nextFeature(feat):
            spIndex.insertFeature(feat)

        SplitPoints = list()
        field_names = [field.name() for field in provider.fields() ]
        # loop throug splitpoint and get closest point
        # from points generated above with create_points_at
        # to get the distance in the attributes
        for feat in splayer.getFeatures():
            cuSplitPoint = feat.geometry().asPoint()

            # QgsSpatialIndex.nearestNeighbor (QgsPoint point, int neighbors)
            nearestIds = spIndex.nearestNeighbor(cuSplitPoint,1) # we need only one neighbour
            featureId = nearestIds[0]
            fit2 = provider.getFeatures(QgsFeatureRequest().setFilterFid(featureId))
            mf = list(fit2)
            # SplitPoints.append(mf[0].attributes())
            fieldDict = dict(zip(field_names, mf[0].attributes()))
            SplitPoints.append(fieldDict)

        virt_layer.commitChanges()

        # # Uncomment if you want to have layer in project
        # proj = QgsProject.instance()
        # proj.addMapLayers([virt_layer])
        # virt_layer.reload()

        return(virt_layer, SplitPoints, layer.source())


    def ProcessLine(self,layer,dgm, LineSource):

        params = {
            'INPUT': layer,
            'CALC_METHOD' : 0,
            'OUTPUT': 'memory:MyAddedGeom'
        }

        proc = processing.run('qgis:exportaddgeometrycolumns', params)

        # add z values to points
        params = {
            'COLUMN_PREFIX' : 'z',
            'INPUT' : proc['OUTPUT'],
            'OUTPUT' : 'TEMPORARY_OUTPUT',
            'RASTERCOPY' : dgm
        }

        withZ = processing.run('qgis:rastersampling', params)['OUTPUT']

        # rename output column to z 
        fields = withZ.fields()
        for field in fields:
            if field.name() == 'z_1':
                with edit(withZ):
                    idx = fields.indexFromName(field.name())
                    withZ.renameAttribute(idx, 'z')


        return(withZ)


# https://docs.qgis.org/testing/en/docs/pyqgis_developer_cookbook/vector.html

# Used to develop together with plugin SCRIPT RUNNER
def run_script(iface):
    print("Script")
    ProfileLayer = "LineTest"
    DGMSource = 'D:/Current/AlphaBetaQgis/AB_Test/2_Originaldaten/Dgm/DGM_5m_p.asc'

    prepLine, SplitPoints, LineSource = PrepareLine(ProfileLayer)

    procLine = ProcessLine(prepLine,DGMSource, LineSource)

    AlphaBetaMain(procLine, SplitPoints)
