import argparse
import os
import numpy as np
from tqdm import tqdm
from ReadLandsat import TOA_REFLECTANCE
from Collocate import Collocation


def run(input_dir, output_dir1, output_dir2, loc_grid, loc_cloud, region, resolution, ncritical):
    print("GENERATION OF LANSAT SIC | REGION: {} | RESOLUTION: {} | N_CRITICAL: {}".format(region, resolution, ncritical))
    # =============================================================================
    # Read the Polar Stereographic Grid System: The system we want to project to
    # =============================================================================
    cls_col = Collocation(region, resolution, ncritical, output_dir1, output_dir2, loc_grid)
    lonpsn = cls_col.readgrid(type = 'lons')
    latpsn = cls_col.readgrid(type = 'lats')

    # =============================================================================
    # Read the Landsat Band Files: The system we want to project
    # =============================================================================
    loc_region = os.path.join(input_dir, region)
    loc_dates = [os.path.join(loc_region, date) for date in os.listdir(loc_region)]
    path_list = []
    for n in range(len(loc_dates)):
        files = os.listdir(loc_dates[n])
        files.sort()
        path_file = []
        for file in files:
            path_file.append(os.path.join(loc_dates[n], file))
        path_list.append(path_file)
    # Total Number of Files
    total_size = 0
    for n in range(len(path_list)):
        total_size += len(path_list[n])
    
    count = 1
    for n in range(len(path_list)):
        for m in range(len(path_list[n])):
            path = path_list[n][m]
            filename = os.path.basename(path)
            cls5 = TOA_REFLECTANCE(5, region, path, loc_cloud)
            cls6 = TOA_REFLECTANCE(6, region, path, loc_cloud)
            criteria = cls5.readcloudlabel()
            _ = cls6.readcloudlabel()
            if criteria:
    # =============================================================================
    # Classifying Ice / Water using NDSI
    # =============================================================================
                R5 = cls5.ReadLandsat()
                R6 = cls6.ReadLandsat()
                lon_l, lat_l = cls5.utm2lonlat()
                print('Classifying Ice / Water...')
                NDSI = (R5-R6)/(R5+R6)
                SI_landsat = np.zeros_like(NDSI)
                for i in tqdm(range(SI_landsat.shape[0])):
                    for j in range(SI_landsat.shape[1]):
                        if (NDSI[i,j] > 0.45) and (R5[i,j] > 0.08):
                            SI_landsat[i,j] = 1
                idx_nan = np.where(np.isnan(R5), True, False)
                SI_landsat[idx_nan] = np.nan
    # =============================================================================
    # Calculating Sea Ice Concentration
    # =============================================================================
                cls_col.collocate(lonpsn, latpsn, lon_l, lat_l, SI_landsat, filename, filename[17:25])
                print('[{} / {}] GENERATED SIC for {}'.format(count, total_size, filename))
                count += 1      
    print("GENERATED LANSAT SIC for | REGION: {} | RESOLUTION: {} | N_CRITICAL: {}".format(region, resolution, ncritical))

def parse_opt():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type = str, default = './input', help = 'Directory of Input Files') # For Me: '/home/junghs/Sea_Ice/data/LANDSAT'
    parser.add_argument('--output_dir1', type = str, default = '/output1', help = 'Directory of Output Files: SIC') # For me: '/home/junghs/Sea_Ice/data/LANDSAT_gridded625'
    parser.add_argument('--output_dir2', type = str, default = '/output2', help = 'Directory of Output Files: Statistics') # For me: '/home/junghs/Sea_Ice/data/LANDSAT_stats625'
    parser.add_argument('--loc_grid', type = str, default = './grid', help = 'Directory of Grid system') # For Me: '/home/junghs/Sea_Ice/data/GRID'
    parser.add_argument('--loc_cloud', type = str, default = '/cloud_label', help = 'Directory of Cloud Contamination Label') # For Me: '/home/junghs/Sea_Ice/data/CloudMask_validation/eye_comparison
    parser.add_argument('--region', type = str, required = True, help = 'Region Name')
    parser.add_argument('--resolution', type = float, required = True, help = 'Grid Resolution: 6.25, 12.5, 25 km')
    parser.add_argument('--ncritical', type = int, default = 15000, help = 'Critical Sample Size') # 15000 for 6.25 km resolution
    opt = parser.parse_args()
    return opt

def main(opt):
    run(**vars(opt))

if __name__ == '__main__':
    opt = parse_opt()
    main(opt)
    



    

