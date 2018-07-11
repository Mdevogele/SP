#!/usr/bin/env python

""" SP_Prepare2 - Substract background to science data  
    v1.0: 2018-04-19, mdevogele@lowell.edu
  
"""


import argparse, shlex

import os
import logging
import datetime
from astropy.io import ascii
from astropy.io import fits

import SP_Toolbox as toolbox
from SP_CheckObsType import CheckObsType
import _SP_conf
from SP_CheckInstrument import CheckInstrument

# setup logging
logging.basicConfig(filename = _SP_conf.log_filename,
                    level    = _SP_conf.log_level,
                    format   = _SP_conf.log_formatline,
                    datefmt  = _SP_conf.log_datefmt)


def Prepare(filenames,Verbose):
    
    
    # Get current directory 
    
    Directory = os.getcwd()
    
    # Removes the './' if present in filenames
    
    for idx, filename in enumerate(filenames):
        filenames[idx] = filename.replace('./','')
    
    # start logging
    logging.info('Preparing data with parameters: %s',filenames)
    print('Preparing data with parameters: %s',filenames)

    # Create the procc folder 

#    now = datetime.datetime.now()

#    Dir_To_Create = 'Procc_' + str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '_' + str(now.hour) + '_' + str(now.minute)           + '_' + str(now.second)
#    os.mkdir(Dir_To_Create)
#    Proc_Dir = Dir_To_Create
    
    
    #  Copy raw data to a /raw directory for safety
    
    if not os.path.exists('raw'):
        os.makedirs('raw')
        for idx, filename in enumerate(filenames):
            shutil.copy(filename, './' + 'raw/' + filename)


    telescope, obsparam = CheckInstrument(filenames)
    # Copy raw file to Procc folders (These files will be modified)
    
#    for idx, filename in enumerate(filenames):
#        shutil.copy('./' + 'raw/' + filename, './' + Dir_To_Create + '/' + filename)
        
    # Move the working directory to the processing directory
    
#    os.chdir(Directory + '/' + Dir_To_Create)
    
    # change FITS file extensions to .fits
    for idx, filename in enumerate(filenames):
        if filename.split('.')[-1] in ['fts', 'FTS', 'FITS', 'fit', 'FIT']:
            os.rename(filename, '.'.join(filename.split('.')[:-1])+'.fits')
            filenames[idx] = '.'.join(filename.split('.')[:-1])+'.fits'
            logging.info('change filename from "%s" to "%s"' %
                         (filename, filenames[idx]))
            
    # identify keywords for GENERIC telescopes

    # open one sample image file
    hdulist = fits.open(filenames[0], verify='ignore',
                        ignore_missing_end='True')
    header = hdulist[0].header
    
    # keywords that have to be implanted into each image
    implants = {}
    
    # prepare image headers for spectroscopic pipeline
    
    
    CheckObsType(filenames)
    
#    for filename in filenames:

#        if display:
#            print('preparing', filename)
            
        # open image file
#        hdulist = fits.open(filename, mode='update', verify='silentfix',
#                            ignore_missing_end=True)
#        header = hdulist[0].header

        # add other headers, if available
        
#        if len(hdulist) > 1:
#            for i in range(len(hdulist)):
#                try:
#                    header += hdulist[i].header
#                except:
#                    pass
                
        # read image data
#        imdata = hdulist[0].data
    
        # read out image binning mode
#        binning = toolbox.get_binning(header, obsparam)
    
     # create diagnostics
#    if diagnostics:
#        diag.create_index(filenames, os.getcwd(), obsparam, display)

#    logging.info('Done! -----------------------------------------------------')

    return None   


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Prepare the files for the spectroscopic pipeline')
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
#    parser.add_argument('-r',
#                        help='Range of pixels to use for background subtraction',
#                        nargs=2)    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
#    man_targetname = args.target
#    Method = args.m
#    Series = args.s
#    MasterBias = args.b
    Verbose = args.v
#    Range = args.r
    filenames = args.images  
    
    print(filenames)

    
    
    # call run_the_pipeline only on filenames
    Prepare(filenames,Verbose)
    pass





