# Landsat 8 Sea Ice Concentration

Code to read Landsat 8 band images 
and calculate the Sea Ice Concentration (SIC) in 6.25, 12.5, and 25 km Polar Stereographic Format
using methods suggested in Tanaka and Lu (2023)

## Environment Setup

- Python 3.9.13
- anaconda 4.12.0
- Please install the dependencies in **requirements.txt** with the following line
```
conda install -r requirements.txt
```

## How to Run the Code

- Clone the repository
```
git clone https://github.com/hsjung24/Landsat8-Utilities.git
```
- Go into the cloned folder
```
cd Landsat8_Sea_Ice_Concentration
```
- Place the Landsat-8 Level 1 Data files in the **input** folder
```
# Directory Structure should be the same as shown
  1. region: Name of the sub-regions of the Arctic Sea
  2. date: Acquisition date in YYYYMMDD Format

./input/region/date/Landsat_8_Filename
```
- Run **main.py** with specified arguments
```
# For selection of resolution
python main.py --region "region" --resolution 6.25   # 6.25 km resolution
python main.py --region "region" --resolution 25 --ncritical 100000    # 25 km resolution

# For selection of the range of files you want to process
python main.py --region "region" --resolution 6.25 --count_start 50 --count_end 100   # Will Generate SIC for files number 50 ~ 100
```
- Output values will be saved in the **output1** and **output2** folders

  - ***./output1/SIC_Landsat_{filename}.npy***: 2D Numpy array of Landsat-8 Sea Ice Concentration
  - ***./output2/Confidence_Interval_{}.pkl***: Dictionary containing the Confidence Interval Length for each calculated Landsat-8 Sea Ice Concentration



