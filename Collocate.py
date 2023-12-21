from unittest import result
import numpy as np
import os
from tqdm import tqdm
import pyproj
import warnings
from scipy import stats
import pickle


class Collocation(object):

    ''' Initialize Class by Specifying
        1. region
        2. resolution: 6.25, 12.5, 25 km
        3. ncritical: 15000 for 6.25 km resolution / 100000 for 25 km resolution
        3. output_dir1: Output Directory for SIC_Landsat.npy
        4. output_dir2: Output Directory for Confidence_Interval.pkl
        5. loc_grid: Directory of the Polar Stereographic grid system
    '''
    def __init__(self, region, resolution, ncritical = 15000, output_dir1 = './output1', output_dir2 = './output2', loc_grid = './grid'):
        self.region = region
        self.resolution = resolution
        self.ncritical = ncritical
        self.NY = int(448*(25/self.resolution))
        self.NX = int(304*(25/self.resolution))
        self.output_dir1 = output_dir1 # For me: '/home/junghs/Sea_Ice/data/LANDSAT_gridded625'
        self.output_dir2 = output_dir2 # For me: '/home/junghs/Sea_Ice/data/LANDSAT_stats625'
        self.loc_grid = loc_grid # For Me: '/home/junghs/Sea_Ice/data/GRID'

    def collocate(self, lonpsn, latpsn, lon_l, lat_l, SI_landsat, filename, date):
        # 1. Collocation --------------------------------------------------------------------------------------------------------
        print('Collocating on PS Grid System...')
        nx, ny = SI_landsat.shape
        lon_flat = lon_l.reshape(nx*ny); lat_flat = lat_l.reshape(nx*ny); SI_flat = SI_landsat.reshape(nx*ny)
        idx_nan = np.isnan(SI_flat)
        lon_Landsat = lon_flat[~idx_nan]; lat_Landsat = lat_flat[~idx_nan]
        SI_Landsat = SI_flat[~idx_nan]
        p1 = pyproj.Proj('proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs')
        p2 = pyproj.Proj(proj='latlong', datum='WGS84')
        warnings.filterwarnings('ignore')
        x25, y25 = pyproj.transform(p2, p1, lonpsn, latpsn)
        warnings.filterwarnings('ignore')
        x_l, y_l = pyproj.transform(p2, p1, lon_Landsat, lat_Landsat)
        SIC_Landsat = np.zeros([self.NY, self.NX]); SIC_Landsat[:,:] = np.nan
        Samples = []
        for i in tqdm(range(self.NX)):
            if i == 1:
                x_low = x25[0,i]; x_high = 0.5*(x25[0,i] + x25[0,i+1])
            elif i == self.NX-1:
                x_low = 0.5*(x25[0,i] + x25[0,i-1]); x_high = x25[0,i]
            else:
                x_low = 0.5*(x25[0,i-1] + x25[0,i]); x_high = 0.5*(x25[0,i] + x25[0,i+1])
            idx = np.where(x_l >= x_low, True, False) * np.where(x_l < x_high, True, False)
            if any(idx):
                yl_contain = y_l[idx]
                yl_max = np.max(yl_contain); yl_min = np.min(yl_contain)
                y_dif_high = np.squeeze(y25[:,i]) - yl_max
                y_dif_low = np.squeeze(y25[:,i]) - yl_min
                # Y-dir Boundary Index
                idy_max = np.where(np.diff(np.sign(y_dif_high)))[0][0]
                idy_min = np.where(np.diff(np.sign(y_dif_low)))[0][0] + 1
                for j in np.arange(idy_max, idy_min+1, 1):
                    if j == 1:
                        y_high = y25[j,0]
                        y_low = (y_high + y25[j+1,0])/2
                    elif j == self.NY-1:
                        y_low = y25[j,i]
                        y_high = (y_low + y25[j-1,0])/2
                    else:
                        y_high = (y25[j-1,0] + y25[j,0])/2
                        y_low = (y25[j,0] + y25[j+1,0])/2
                    idy = np.where(y_l >= y_low, True, False)*np.where(y_l < y_high, True, False)
                    id_SI = idx*idy
                    N = np.count_nonzero(id_SI)
                    if N >= 1:
                            sample = {}
                            sample['size'] = N
                            sample['distribution'] = SI_Landsat[id_SI]
                            sample['pixel_x_index'] = i
                            sample['pixel_y_index'] = j
                            Samples.append(sample)
                    if N >= self.ncritical: # Point of Saturation for the Bootstrap Confidence Interval
                        SI_contained = SI_Landsat[id_SI]
                        Ns = np.count_nonzero(SI_contained)
                        SIC_Landsat[j,i] = Ns/N*100
                    else:
                        continue       
            else:
                continue

        # 2. Bootstrap Confidence Interval Calculation ------------------------------------------------------------------
        print('Calculaing BS Confidence Interval...')
        ds = {}
        NUM = np.array([item['size'] for item in Samples])
        X_IDX = np.array([item['pixel_x_index'] for item in Samples])
        Y_IDX = np.array([item['pixel_y_index'] for item in Samples])
        DISTRIBUTION = [item['distribution'] for item in Samples]
        num_subsample = 100
        MEAN_SAMPLE = []
        MEAN_RESAMPLE = []
        DELTA025 = []
        DELTA975 = []
        LENGTH_CI = []
        for i in tqdm(range(len(Samples))):
            sample = DISTRIBUTION[i]
            mean_sample = np.mean(sample)
            resample = np.zeros([num_subsample, len(sample)])
            for j in range(num_subsample):
                rs = np.random.choice(sample, len(sample), replace = True)
                resample[j,:] = rs
            mean_resample = np.mean(resample, axis = 1)
            delta_resample = mean_resample - mean_sample
            delta_resample.sort()
            delta_resample025 = np.percentile(delta_resample, 2.5)
            delta_resample975 = np.percentile(delta_resample, 97.5)
            length_CI = delta_resample975 - delta_resample025
            MEAN_SAMPLE.append(mean_sample) # True Mean [1]
            MEAN_RESAMPLE.append(mean_resample) # Distribution of Resampled Mean [100]
            DELTA025.append(delta_resample025)
            DELTA975.append(delta_resample975)
            LENGTH_CI.append(length_CI)
        MEAN_SAMPLE = np.array(MEAN_SAMPLE)*100
        MEAN_RESAMPLE = np.array(MEAN_RESAMPLE)*100
        DELTA025 = np.array(DELTA025)*100
        DELTA975 = np.array(DELTA975)*100
        LENGTH_CI = np.array(LENGTH_CI)*100
        
        ds['filename'] = filename # Landsat Filename
        ds['date'] = date # Time
        ds['pixel_x_index'] = X_IDX # SSMIS Pixel x-Index
        ds['pixel_y_index'] = Y_IDX # SSMIS Pixel y-Index
        ds['true_mean'] = MEAN_SAMPLE # True Mean of Pixel
        ds['confidence_interval'] = LENGTH_CI # Bootstrap 95% CI Length
        ds['sample_size'] = NUM # Sample Size
        ds['ci_low'] = DELTA025 # Lower Bound of CI
        ds['ci_high'] = DELTA975 # Upper Bound of CI
        ds['distribution'] = MEAN_RESAMPLE # Bootstrap Resample Mean Distribution

        # 3. Save SIC as SIC_Landsat_{filename}.npy --------------------------------------------------------------------------
        if not os.path.isdir(os.path.join(self.output_dir1, self.region)):
            path_region = os.path.join(self.output_dir1, self.region)
            os.mkdir(path_region)
        if not os.path.isdir(os.path.join(path_region, date)):
            path_date = os.path.join(path_region, date)
            os.mkdir(path_date)
        with open(os.path.join(path_date, 'SIC_Landsat_{}.npy'.format(filename)), 'wb') as f1:
            np.save(f1, SIC_Landsat)

        # 4. Save Confidence Interval as Confidence_Interval_{filename}.pkl --------------------------------------------------
        if not os.path.isdir(os.path.join(self.output_dir2, self.region)):
            path_region2 = os.path.join(self.output_dir2, self.region)
            os.mkdir(path_region2)
        if not os.path.isdir(os.path.join(path_region2, date)):
            path_date2 = os.path.join(path_region2, date)
            os.mkdir(path_date2)
        with open(os.path.join(path_date2, 'Confidence_Interval_{}.pkl'.format(filename)), 'wb') as f2:
            pickle.dump(ds, f2)

    def readgrid(self, type):
        f = open(os.path.join(self.loc_grid, 'psn06{}_v3.dat'.format(type)), 'rb')
        data = np.fromfile(f, dtype = 'int32')
        data = data.reshape(self.NY, self.NX)
        data = data * 1e-5
        return data