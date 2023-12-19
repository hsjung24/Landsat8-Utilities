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
    def __init__(self, filename, Band_Number, confidence = 'medium'):
        self.filename = filename
        self.Band_Number = Band_Number
        self.confidence = confidence
        self.input_dir = str('./data/'+self.filename+'/')
        self.output_dir = str('./output/'+self.filename+'/')
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
    
    ''' Will Return
        1. Spectral TOA Reflectance: [nx, ny]
        2. Longitude: [nx, ny]
        3. Latitude: [nx, ny]
    '''
    # Retrieve Landsat Band TOA Reflectance
    def ReadLandsat(self):
        # Raw Band Information
        dataset = gdal.Open(self.input_dir + self.filename + '_B{}.TIF'.format(self.Band_Number))
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
        lon, lat = self.utm2lonlat()
        # Save the Files in ./output/Filename
        with open(self.output_dir + '/TOA_Reflectance.npy', 'wb') as f:
            np.save(f, TOA_REFLECTANCE)
        with open(self.output_dir + '/longitude.npy', 'wb') as f:
            np.save(f, lon)
        with open(self.output_dir + '/latitude.npy', 'wb') as f:
            np.save(f, lat)
    

    # Read Metadata
    def readmdata(self, name_meta):
        for name in os.listdir(self.input_dir):
            if name.endswith('MTL.txt'):
                path = self.input_dir + name
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
        dataset = gdal.Open(self.input_dir + self.filename + '_B{}.TIF'.format(self.Band_Number))
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
        return lon_l, lat_l
        
    # Quality Acessment Mask Array   
    def GetMask(self, MaskType, Type = 'PIXEL'):
        dataset = gdal.Open(self.input_dir + self.filename + '_QA_{}.TIF'.format(Type))
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
    
        

