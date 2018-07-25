#!/usr/bin/env python

""" SP_WavCal - Perform the wavelength calibration
    v1.0: 2018-04-19, mdevogele@lowell.edu
"""

import argparse, shlex
import numpy as np
import SP_Toolbox as SP
import SP_Deveny_Toolbox as SP_Dev_Tool
from astropy.io import fits

from SP_CheckInstrument import CheckInstrument


def WavCal(filenames,ArcsFile,OutFile,Verbose):
    
    
    Arcs = []
    for elem in ArcsFile:
        hdulist = fits.open(elem)
        Arcs.append(hdulist[0].data)
       
    Arcs = np.median(Arcs,axis=0)
    
    telescope, obsparam = CheckInstrument([ArcsFile[0]])  
    
    print(telescope)
    
    if telescope == 'GMOSS' or telescope == 'GMOSN':
        Detector = hdulist[0].header['DETECTOR']
        Gratting = hdulist[0].header[obsparam['grating']]
        Binning = hdulist[1].header['CCDSUM'][0]
        DetecFlags = {'Instrument':telescope, 'Binning' : Binning, 'Gratting':Gratting, 'Detector':Detector}  
        print DetecFlags
    if telescope == 'DEVENY':
        DetecFlags = {'Instrument':'Deveny'}  

       
    print(OutFile)
    Wav = SP.Wav_Cal2(Arcs,**DetecFlags)
    Wav = np.array(Wav)
    Wav = Wav/10000
    Wav = np.flip(Wav,axis=0)
    print(Wav)
    
      
    with open(filenames[0],'r') as f:
        SpecA = f.read().splitlines()  
    
    SpecA = np.flip(SpecA,axis=0)
    
    f = open(OutFile,'w')
    for Wave,Spec in zip(Wav,SpecA):
      f.write("{} \t {} \n".format(Wave,Spec))
    
    
    return None


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
    parser.add_argument('-a',
                        help='List of arcs files',
                        nargs='+')
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
    ArcsFile = args.a
    Verbose = args.v
    OutFile = args.o
    filenames = args.images  
    
    print(ArcsFile)
    print(OutFile)

    
    
    # call run_the_pipeline only on filenames
    WavCal(filenames,ArcsFile,OutFile,Verbose)
    pass