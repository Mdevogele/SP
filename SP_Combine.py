#!/usr/bin/env python

""" SP_Combine - Combine Extracted spectra
    v1.0: 2018-04-19, mdevogele@lowell.edu
"""

import argparse, shlex
import numpy as np
import _SP_conf


def Combine_Spectra(filenames,OutFile,Verbose):

    SpecT = [] 
    for elem in filenames:
        with open(elem,'r') as f:
            SpecT.append(f.read().splitlines())
    SpecT = np.array(SpecT).astype(float)
    SpecT[SpecT<=0] = np.nan
    SpecT[SpecT>10] = np.nan
    for idx,elem in enumerate(SpecT):
        SpecT[idx] = SpecT[idx]/np.nanmedian(SpecT[idx][1300:1500])
    SpecTT = np.nanmedian(SpecT,axis=0)
    f = open(OutFile,'w')
    for item in SpecTT:
      f.write("%s\n" % item)
    f.close()
    
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract spectrum')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
#    parser.add_argument('-target', help='primary targetname override',
#                        default=None)
#    parser.add_argument('-m', help='add flats information to diagnostic.htlm file',
#                        choices=['auto','range'])
#    parser.add_argument('-s', help='If there is several series of flat \n || Options: none: Only one serie \n || index: split according to the index of the files \n || target: split according to the target name in fits headers \n || pointing: split according to telescope pointing according to fits headers',
#                        choices=['none','index','target','pointing'],
#                        default = 'None')
#    parser.add_argument('-b',
#                        help='Name of the master bias to use \n || Can use None if no bias are available',
#                        default='MasterBias.fits')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-o',
                        help='Name of the combined spectrum',
                        default = 'SpecOut.spec')    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
#    man_targetname = args.target
#    Method = args.m
#    Series = args.s
#    MasterBias = args.b
    Verbose = args.v
    OutFile = args.o
    filenames = args.images  
    
    print(filenames)

    
    
    # call run_the_pipeline only on filenames
    Combine_Spectra(filenames,OutFile,Verbose)
    pass
    