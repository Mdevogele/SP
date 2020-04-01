#!/usr/bin/env python

""" SP_Bias - wrapper for creating bias 
    v1.0: 2018-04-17, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import numpy as np

from astropy.io import fits
from astropy.time import Time


from SP_CheckInstrument import CheckInstrument



import SP_diagnostics as diag

def Create_Bias(filenames,MasterName,Verbose,Diagnostic,Method):
    
    
    if Verbose:
        print('Beginning bias processing')
        print('Processing files:')
        print('index \t filename')
        for idx,elem in enumerate(filenames):
            print('{} \t {}'.format(idx+1,elem))

  
    if Verbose:
        print('Creating the master bias')
           
    Bias = []
    times  = []

    for image in filenames:
        toopen = image
        hdulist = fits.open(toopen)
        times.append(hdulist[0].header['DATE-OBS'])        
        Bias.append(hdulist[0].data)

    t = Time(times, format='isot', scale='utc')
    MidTime = t[0] + (t[-1]-t[0])/2
        
    hdulist[0].data
    
    if Method == 'Median':
        MasterBias = np.median(Bias,axis = 0 )
    elif Method == 'Mean':
        MasterBias = np.mean(Bias,axis = 0 )
    else:
        print('Method ' + Method + ' was not recognized, Median used by default')
        MasterBias = np.median(Bias,axis = 0 )
    
    hdulist[0].data = MasterBias

    print(MasterName)
    
    
    ## Header updates 
    hdulist[0].header['MIDTIME'] = (MidTime.isot, 'SP: midtime of the individual biases')
    hdulist[0].header['CREATIME'] = (Time.now().isot, 'SP: Creation time of the Master bias')
    hdulist[0].header['NUMDARKS'] = (len(filenames), 'SP: Number of individual biases used')
    hdulist[0].header['PIPELINE'] = ('Spectroscopic Pipeline', 'SP: pipeline used to created this file')
    hdulist[0].header['PROCTYPE'] = ('MASTER BIAS', 'SP: Type of processed file')
    hdulist[0].header['FILENAME'] = (MasterName, 'SP: Name of the file')
    hdulist[0].header['METHOD'] = (Method,'Method used to combine the individual biases')
    
    hdulist[0].header['HISTORY'] =  Method + ' master bias created from ' +  str(len(filenames)) + ' individual biases'
    hdulist[0].header['HISTORY'] =  'on ' + str(Time.now().isot) + ' (YYYY-MM-DD hh:mm:ss UT) from'
    for elem in filenames:
        hdulist[0].header['HISTORY'] = elem
        
   
    
    
    hdulist.writeto(MasterName, overwrite = True)
    hdulist.close()
    
    if Verbose:
        print('Master bias save to {}'.format(MasterName))
        hdulist = fits.open(MasterName)
        data = hdulist[0].data
        print('Statistics of the Master bias: ')
        print('Mean: {} \t Median: {} \t std: {}'.format(np.nanmean(data), np.nanmedian(data), np.nanstd(data)))
        print('End of bias processing')

    if Diagnostic:
        diag.create_website('Bias_Log.html')
        diag.add_BiasSummary(filenames,MasterName,'Bias_Log.html')
        diag.add_BiasList(filenames,'Bias_Log.html')        
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master bias')

    parser.add_argument('-lv',
                        help='decrease verbosity',
                        action="store_false",
                        default = True)    
    parser.add_argument('-o',
                        help='Name of the master bias file',
                        default='MasterBias.fits')  
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)   
    parser.add_argument('-m',
                        help='Method used to combine data',
                        default = 'Median')   
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()

    Verbose = args.lv
    MasterName = args.o
    filenames = args.images    
    Diagnostic = args.d
    Method = args.m

    SP.Check(filenames)

    Create_Bias(filenames,MasterName,Verbose,Diagnostic,Method)
    pass
