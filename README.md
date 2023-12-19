# Landsat 8 Utilities

Code to read Landsat 8 band images and return the top-of-atmosphere (TOA) reflectance values in Numpy array format

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
git clone 
```
- Go into the cloned folder
```
cd Landsat8-Utilities
```
- Place the Landsat-8 Level 1 Data files in the **data** folder
```
./data/Landsat_8_Filename
```
- Run **main.py** with specified arguments
```
# For default Cloud Mask Confidence (Medium)
python main.py --filename "filename" --Band_Number "Band_Number"

# For selection of Cloud Mask Confidence
python main.py --filename "filename" --Band_Number "Band_Number" --confidence "low"   #low confidence
python main.py --filename "filename" --Band_Number "Band_Number" --confidence "high"  #high confidence
```
- Output Numpy arrays will be saved in the **output** folder under the same name as the input file containing...
  - ***TOA_Reflectance.npy***: Spectral op-of-atmosphere reflectance
  - ***longitude.npy***: Longitude coordinate of each pixel
  - ***Latitude.npy**: Latitude coordinate of each pixel
```
./output/Landsat_8_Filename
```



