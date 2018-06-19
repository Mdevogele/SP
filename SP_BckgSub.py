#!/usr/bin/env python

""" SP_BckgSub - Substract background to science data  
    v1.0: 2018-04-19, mdevogele@lowell.edu
  
"""
import argparse, shlex

from astropy.io import fits
import SP_Toolbox as SP

from SP_CheckInstrument import CheckInstrument

import numpy as np

def BckgSub(FileName, Verbose, Method ,Suffix,Area = [250,350]):
   
    telescope, obsparam = CheckInstrument([FileName[0]])
    print telescope
    for elem in FileName:
        
        hdulist = fits.open(elem)
        image = hdulist[0].data
        Mask = np.zeros_like(image)
        Mask[:,:] =  0
        
        if Method == 'auto':
            print telescope
            if telescope == 'Deveny' or 'DEVENY':
                DetecFlags = {'Instrument':'Deveny'}
            if telescope == 'GMOSS' or telescope == 'GMOSN':
                DetecFlags = {'Instrument':'GMOS'}
                
            Center = SP.Detect_Spectra(image,Bin =2,**DetecFlags)
            print(Center)
            Area[0] = int(Center-70)
            Area[1] = int(Center+70)                
        
        print(Area)
        
        Mask[0:Area[0],:] = 1
        Mask[Area[1]:,:] = 1
        Mask = Mask.astype(bool)
        
        image2 = np.ma.masked_array(image,mask = Mask)
        
        X = range(image.shape[0])
        X = np.array(X)
        
        
        for i in range(image.shape[1]):
            index = np.argwhere(np.isnan(image2[:,i]))
            image2[index,i]=0
            image3 = SP.Sigma_Clip(image2[:,i])
            XX = np.ma.masked_array(X,image3.mask)
            z = np.ma.polyfit(XX,image3 , 1)
            p = np.poly1d(z)
            image[:,i] = image[:,i] - p(X)
        
        hdulist.writeto(elem.replace('.fits','').replace('_Procc','') + '_' + Suffix + '.fits')



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master flats')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
    parser.add_argument('-s',
                        help='Suffix to add to processed files',
                        default='Bckg')
    parser.add_argument('-m', help='Method to use for the selection of region to be considered for the evaluation of the background',
                        choices=['auto','range'],
                        default = 'auto')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-r',
                        help='Range of pixels to use for background subtraction. 2 arguments, the pixel center and the number of pixels to consider',
                        default = '300 200',
                        nargs=2)    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
    Suffix = args.s
    Method = args.m
    Verbose = args.v
    Range = args.r
    filenames = args.images    
    
    Ran = []
    Ran.append(int(Range[0])-int(Range[1])/2)
    Ran.append(int(Range[0])+int(Range[1])/2)
    
    
    # call run_the_pipeline only on filenames
    BckgSub(filenames,Verbose,Method,Suffix,Ran)
    pass