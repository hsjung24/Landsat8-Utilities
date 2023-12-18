import argparse
from ReadLandsat import TOA_REFLECTANCE

def parse_opt():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename', type=str, required = True, help='Name of Landsat File')
    parser.add_argument('--Band_Number', type=str, required = True, help='Band Number: 1~9')
    parser.add_argument('--confidence', type=str, default='medium', help='(optional) Confidence Level of Mask: high, medium, low')
    opt = parser.parse_args()
    return opt

def main(opt):
    cls = TOA_REFLECTANCE(**vars(opt))
    cls.ReadLandsat()

if __name__ == '__main__':
    opt = parse_opt()
    main(opt)

