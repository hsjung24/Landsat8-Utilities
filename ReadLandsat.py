#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 19:38:52 2023

@author: junghs
"""

import os
import numpy as np
from osgeo import gdal
import pyproj
import warnings

class TOA_REFLECTANCE(object):
    
    ''' Initialize Class by Specifying
        1. Filename
        2. Band Number: 1~9
        3. Mask Confidence: 'high', 'medium', 'low'
    '''
    def __init__(self, Band_Number, region, input_dir = './input', loc_cloud = './cloud_label'):
        self.Band_Number = Band_Number
        self.confidence = None
        self.region = region
        self.input_dir = input_dir # For Me: '/home/junghs/Sea_Ice/data/LANDSAT
        self.filename = os.path.basename(input_dir)
        self.loc_cloud = loc_cloud # For Me: '/home/junghs/Sea_Ice/data/CloudMask_validation/eye_comparison
        
    ''' Will Return
        1. Spectral TOA Reflectance: [nx, ny]
        2. Longitude: [nx, ny]
        3. Latitude: [nx, ny] 
    '''
    # Retrieve Landsat Band TOA Reflectance
    def ReadLandsat(self):
        # Raw Band Information
        print("Calculating TOA Reflectance for Band {}...".format(self.Band_Number))
        dataset = gdal.Open(os.path.join(self.input_dir, self.filename +'_B{}.TIF'.format(self.Band_Number)))
        band = dataset.GetRasterBand(1)
        raw_band = band.ReadAsArray()
        # Quality Pixel & Radiometric Saturation Mask
        fill_mask = self.GetMask('fill')
        dilated_cloud_mask = self.GetMask('dilated_cloud')
        b5_saturation_mask = self.GetMask('band5', Type = 'RADSAT')
        b6_saturation_mask = self.GetMask('band6', Type = "RADSAT")
        terrain_occlusion_mask = self.GetMask('terrain_occlusion', Type = 'RADSAT')
        if self.confidence == 'medium':
            cloud_mask = self.GetMask('medium_confidence_cloud')
            cloud_shadow_mask = self.GetMask('reserved_cloud_shadow')
            cirrus_mask = self.GetMask('reserved_cirrus')
        elif self.confidence == 'high':
            cloud_mask = self.GetMask('high_confidence_cloud')
            cloud_shadow_mask = self.GetMask('high_confidence_cloud_shadow')
            cirrus_mask = self.GetMask('high_confidence_cirrus')
        elif self.confidence == 'low':
            cloud_mask = self.GetMask('low_confidence_cloud')
            cloud_shadow_mask = self.GetMask('low_confidence_cloud_shadow')
            cirrus_mask = self.GetMask('low_confidence_cirrus')
        total_mask = fill_mask*cloud_mask*cloud_shadow_mask*dilated_cloud_mask*cirrus_mask*b5_saturation_mask*b6_saturation_mask*terrain_occlusion_mask
        BAND_DN = raw_band * total_mask
        # Calculation of TOA Reflectance (ZANTER 2019)
        theta_sol = self.readmdata('SUN_ELEVATION')
        M_sf = self.readmdata('REFLECTANCE_MULT_BAND_{}'.format(self.Band_Number))
        A_sf = self.readmdata('REFLECTANCE_ADD_BAND_{}'.format(self.Band_Number))
        TOA_REFLECTANCE = (M_sf*BAND_DN+A_sf)/np.sin(np.pi/180*theta_sol)
        return TOA_REFLECTANCE


    # Read Metadata
    def readmdata(self, name_meta):
        for name in os.listdir(self.input_dir):
            if name.endswith('MTL.txt'):
                path = os.path.join(self.input_dir, name)
        f = open(path, 'r')
        for line in f:
            l = line.strip()
            v = l.split('=')
            if len(v) > 1:
                header = v[0].strip()
                contant = v[1].strip()
                if header == name_meta:
                    output = contant
                    break
        f.close()
        return float(output)
    

    # Convert UTM Grids into Longitude and Latitude
    def utm2lonlat(self):
        warnings.filterwarnings('ignore')
        p1 = pyproj.Proj(proj = 'utm', zone = self.readmdata('UTM_ZONE'), datum = 'WGS84', south = False)
        p2 = pyproj.Proj(proj = 'latlong', datum = 'WGS84')
        dataset = gdal.Open(os.path.join(self.input_dir, self.filename + '_B{}.TIF'.format(self.Band_Number)))
        ulx, xscale, _, uly, _, yscale = dataset.GetGeoTransform()
        utmx = []; utmy = []; nx = dataset.RasterXSize; ny = dataset.RasterYSize
        for i in range(nx):
            utmx.append(ulx + xscale*(i-1))
        for i in range(ny):
            utmy.append(uly + yscale*(i-1))
        utmx = np.array(utmx)
        utmy = np.array(utmy)
        xx, yy = np.meshgrid(utmx, utmy); x = xx.reshape(nx*ny); y = yy.reshape(nx*ny)
        lon_landsat, lat_landsat = pyproj.transform(p1, p2, x, y); lon_l = lon_landsat.reshape(ny,nx); lat_l = lat_landsat.reshape(ny,nx)
        # Quality Pixel & Radiometric Saturation Mask
        fill_mask = self.GetMask('fill')
        dilated_cloud_mask = self.GetMask('dilated_cloud')
        b5_saturation_mask = self.GetMask('band5', Type = 'RADSAT')
        b6_saturation_mask = self.GetMask('band6', Type = "RADSAT")
        terrain_occlusion_mask = self.GetMask('terrain_occlusion', Type = 'RADSAT')
        if self.confidence == 'medium':
            cloud_mask = self.GetMask('medium_confidence_cloud')
            cloud_shadow_mask = self.GetMask('reserved_cloud_shadow')
            cirrus_mask = self.GetMask('reserved_cirrus')
        elif self.confidence == 'high':
            cloud_mask = self.GetMask('high_confidence_cloud')
            cloud_shadow_mask = self.GetMask('high_confidence_cloud_shadow')
            cirrus_mask = self.GetMask('high_confidence_cirrus')
        elif self.confidence == 'low':
            cloud_mask = self.GetMask('low_confidence_cloud')
            cloud_shadow_mask = self.GetMask('low_confidence_cloud_shadow')
            cirrus_mask = self.GetMask('low_confidence_cirrus')
        total_mask = fill_mask*cloud_mask*cloud_shadow_mask*dilated_cloud_mask*cirrus_mask*b5_saturation_mask*b6_saturation_mask*terrain_occlusion_mask
        return lon_l*total_mask, lat_l*total_mask
    

    # Quality Acessment Mask Array   
    def GetMask(self, MaskType, Type = 'PIXEL'):
        dataset = gdal.Open(os.path.join(self.input_dir, self.filename + '_QA_{}.TIF'.format(Type)))
        band = dataset.GetRasterBand(1)
        raw_band = band.ReadAsArray()
        Q = np.array(raw_band)
        bits = {'fill': 0, 'dilated_cloud': 1, 
                'high_confidence_cirrus': 2, 'low_confidence_cirrus': 14, 'reserved_cirrus': 15,
                'high_confidence_cloud': 3, 'low_confidence_cloud': 8, 'medium_confidence_cloud': 9, 
                'high_confidence_cloud_shadow': 4, 'low_confidence_cloud_shadow': 10, 'reserved_cloud_shadow': 11, 
                'high_confidence_snow': 5, 'low_confidence_snow': 12, 'reserved_snow': 13, 
                'water': 7,
                'band5': 4, 'band6': 5, 'terrain_occlusion': 11}
        bit = bits[MaskType]
        val = 1*(2**bit)
        quotient = Q//val
        MaskArray = np.where(quotient%2 == 1, np.nan, 1)       
        if MaskType == 'snow' or MaskType == 'water': # 1 - Water/Snow | 0 - not Water/not Snow
            MaskArray_s = np.where(quotient%2 == 1, 1, 0)
            return MaskArray_s
        else:    
            return MaskArray
    

    # Cloud Contamination Label
    def readcloudlabel(self):
        # Cloud Contamination Labels 1
        loc_cloud = os.path.join(self.loc_cloud, '{}_Cloud_Contamination_Label.txt'.format(self.region))
        ccl = open(loc_cloud, 'r').read().split('\n')
        ccl = ccl[5:]
        underestimated = []
        for item in ccl:
            if len(item) != 0:
                if item[41] == str(1):
                    underestimated.append(item[:40])
        underestimated = np.array(underestimated)
        # Cloud Contamination Labels 2
        loc_cloud2 = os.path.join(self.loc_cloud, '{}_Cloud_Contamination_Label2.txt'.format(self.region))
        ccl2 = open(loc_cloud2, 'r').read().split('\n')
        overestimated = []
        for item in ccl2:
            if len(item) != 0:
                overestimated.append(item)
        overestimated = np.array(overestimated)

        idx_under = np.where(underestimated == self.filename, True, False)
        idx_over = np.where(overestimated == self.filename, True, False)
        if all(~idx_under):
            criteria = True
            if all(~idx_over):
                self.confidence = 'medium'
            else:
                self.confidence = 'high'
        else:
            criteria = False
        return criteria

