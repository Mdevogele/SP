#!/usr/bin/env python

""" SP_WavCal - Perform the wavelength calibration
    v1.0: 2018-04-19, mdevogele@lowell.edu
"""

import argparse, shlex
import numpy as np
import SP_Toolbox as SP
from astropy.io import fits
import pickle
import operator
from itertools import combinations
import os
import _SP_conf

from SP_CheckInstrument import CheckInstrument


def WavCal(filenames,ArcsFile,OutFile,Verbose,Method,Line,Full):
    
#    logging.info('***************************************')
#    logging.info('****** Start of SP_Wavcal script ******')
#    logging.info('***************************************')    
    
    
    Pipe_Path = _SP_conf.rootpath
    Arcs = []
    for elem in ArcsFile:
        hdulist = fits.open(elem)
        Arcs.append(hdulist[0].data)
       
    Arcs = np.median(Arcs,axis=0)
    
    telescope, obsparam = CheckInstrument([ArcsFile[0]])  
    
    print(telescope)
    
    SpecA = np.loadtxt(filenames[0]).transpose()
    SpecA = np.array(SpecA)
    
    if Method == 'template':
    
    
        if telescope == 'GMOSS' or telescope == 'GMOSN':
            Detector = hdulist[0].header['DETECTOR']
            print(Detector)
            Gratting = hdulist[0].header[obsparam['grating']]
            Binning = hdulist[1].header['CCDSUM'][0]
            DetecFlags = {'Instrument':telescope, 'Binning' : Binning, 'Gratting':Gratting, 'Detector':Detector}  
            print DetecFlags
        if telescope == 'DEVENY':
            DetecFlags = {'Instrument':'Deveny'}  
        if telescope == 'SOAR':
            DetecFlags = {'Instrument':'Soar'}      
        if telescope == 'NOT':
            Wav = ((10.93/2)*np.array(range(np.shape(Arcs)[1]))+3609.8)/10000
            print(np.shape(Arcs)[1])
            np.savetxt(OutFile,np.array([Wav,SpecA[0,:],SpecA[1,:]]).transpose())

            return None
    
           
        print(DetecFlags)
        Wav = SP.Wav_Cal2(Arcs,**DetecFlags)
        Wav = np.array(Wav)
        Wav = Wav/10000
#        Wav = np.flip(Wav,axis=0)
        print(Wav)
        
          
#        with open(filenames[0],'r') as f:
#            SpecA = f.read().splitlines()  
        
#        SpecA = np.loadtxt(filenames[0]).transpose()
#        SpecA = np.array(SpecA)
#        SpecA = np.flip(SpecA,axis=0)
        
#        f = open(OutFile,'w')
#        for Wave,Spec in zip(Wav,SpecA):
#          f.write("{} \t {} \n".format(Wave,Spec))
          
        np.savetxt(OutFile,np.array([Wav,SpecA[0,:],SpecA[1,:]]).transpose())
          
    if Method == 'auto':
        Dim = []

        Dim.append(np.size(Arcs,0))
        Dim.append(np.size(Arcs,1))
    
        Arcs_L = Arcs[Line,:] 
        Arcs_Loc = SP.Auto_Detect_Lines(Arcs_L, Tresh_Det = 1.35, Tresh_Arcs = [2, 12] )
        
        
        if telescope == 'DEVENY':
            f = open(Pipe_Path +'/Pre_Comp_Deveny')
            Pre = pickle.load(f)
            f.close()
            WV = [7503.9,7635.1,7723.8,7948.2,8014.8,8115.3,8264.5,8424.6,8521.4,9122.9,9224.5,9657.8]
#            WV = [3261.05, 3610.51, 3650.15, 4046.56, 4358.33, 4678.16, 4799.92, 5085.82, 5460.74, 5769.6, 5790.7, 6965.5, 7067.2, 7272.9, 7384.0]
        else:
            f = open(Pipe_Path + '/Wav_Precomp','r')
            Pre = pickle.load(f)
            f.close()
            WV = [7067.2,7147.0, 7272.9, 7384.0 ,7503.9, 7635.1, 7723.8, 7948.2, 8014.8, 8115.3, 8264.5, 8408.2, 8521.4, 8667.9, 9122.9, 9224.5,9354.2,9657.8, 9784.5] 
        
        One = Arcs_Loc[:9]
        One = np.array(One).astype(float)
        O = One-np.max(One)
        O = O/np.min(O)

        D1_one = np.sort(O)[1]
        D2_one = np.sort(O)[2]
        D3_one = np.sort(O)[3]
        D4_one = np.sort(O)[4]
        D5_one = np.sort(O)[5]
        D6_one = np.sort(O)[6]
        D7_one = np.sort(O)[7]
        
        dist = []
        for elem1,elem2,elem3,elem4,elem5,elem6,elem7 in zip(Pre[0],Pre[1],Pre[2],Pre[3],Pre[4],Pre[5],Pre[6]):
            dist.append((elem1-D1_one)**2 + (elem2-D2_one)**2 + (elem3-D3_one)**2 + (elem4-D4_one)**2 + (elem5-D5_one)**2 + (elem6-D6_one)**2 + (elem7-D7_one)**2)

        AS = np.argsort(dist)
        max_index, max_value = min(enumerate(dist), key=operator.itemgetter(1))

#        WV = [7067.2,7147.0, 7272.9, 7384.0 ,7503.9, 7635.1, 7723.8, 7948.2, 8014.8, 8115.3, 8264.5, 8408.2, 8521.4, 8667.9, 9122.9, 9224.5,9354.2,9657.8, 9784.5] 

        comb = combinations(WV, 9)

        WV_Comb = []

        for elem in comb:
            WV_Comb.append(elem)


        a,b = np.polyfit(WV_Comb[AS[0]],np.sort(One)[::-1],1)

#        Red = a*np.array(WV)+b
#        WV_Sol = []
#        for elem in Arcs_Loc[:9]: 
#            Dist = (Red-elem)**2
#            SD = np.argsort(Dist)
#            WV_Sol.append(WV[SD[0]])
#
#        aa= np.polyfit(WV_Sol,Arcs_Loc[:9],1,full=True)
#        
#        print(aa[1])

        with open(filenames[0],'r') as f:
            SpecA = f.read().splitlines()  
        
        SpecA = np.flip(SpecA,axis=0)

        print('bla')
        Wav = (np.array(range(0,len(SpecA))) - b)/a
        Wav = Wav/10000
        Wavel = Wav[::-1]
        
#        Wav = (np.array(range(0,len(SpecA))) - aa[0][1])/aa[0][0]
#        Wav = Wav/10000
#        Wavel = Wav[::-1]
        
        f = open(OutFile,'w')
        for Wave,Spec in zip(Wavel,SpecA):
          f.write("{} \t {} \n".format(Wave,Spec))
        
        
#        
#    logging.info('*************************************')
#    logging.info('****** End of SP_Wavcal script ******')
#    logging.info('*************************************')
#        
        
    return None


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract spectrum')

    parser.add_argument('-a',
                        help='List of arcs files',
                        nargs='+')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    
    parser.add_argument('-o',
                        help='Name of the combined spectrum',
                        default = 'SpecOut.spec')   
    
    parser.add_argument('-m',
                        help='Method to use for calibration: auto or template',
                        default = 'auto')   
    
    parser.add_argument('-l',
                        help='Line to use for wavelength calibration',
                        default = 250)

    parser.add_argument('-f',
                        help='Do full frame wavelength calibration',
                        action="store_true",
                        default = False)      
    
    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
    ArcsFile = args.a
    Verbose = args.v
    OutFile = args.o
    Method = args.m
    filenames = args.images
    full = args.f
    
    Line = int(args.l)
    
    print(ArcsFile)
    print(OutFile)

    
    WavCal(filenames,ArcsFile,OutFile,Verbose,Method,Line,full)
    pass