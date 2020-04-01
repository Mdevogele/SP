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

import shutil

import SP_Toolbox as toolbox
from SP_CheckObsType import CheckObsType
import _SP_conf
from SP_CheckInstrument import CheckInstrument


def Prepare(filenames,Verbose=True):
    
    logging.info('****************************************')
    logging.info('****** Start of SP_Prepare script ******')
    logging.info('****************************************')

    # Get current directory 
    
    Directory = os.getcwd()
    
    # Removes the './' if present in filenames
    
    for idx, filename in enumerate(filenames):
        filenames[idx] = filename.replace('./','')
    
    # start logging
#    logging.info('Preparing data with parameters: %s',filenames)

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
            logging.info('All data have been backup to the ./raw directory')


    telescope, obsparam = CheckInstrument(filenames)


    # change FITS file extensions to .fits
    for idx, filename in enumerate(filenames):
        if filename.split('.')[-1] in ['fts', 'FTS', 'FITS', 'fit', 'FIT']:
            os.rename(filename, '.'.join(filename.split('.')[:-1])+'.fits')
            filenames[idx] = '.'.join(filename.split('.')[:-1])+'.fits'
            logging.info('change filename from "%s" to "%s"' %
                         (filename, filenames[idx]))
#        if telescope == 'SOAR':
#            hdulist = fits.open(filename)
#            hdulist[0].data = hdulist[0].data[200:-200,:]
#            hdulist.writeto(filename,overwrite = True,checksum=False, output_verify='ignore')



    
    CheckObsType(filenames,telescope,obsparam)


    logging.info('**************************************')
    logging.info('****** end of SP_Prepare script ******')
    logging.info('**************************************')


    return None   


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Prepare the files for the spectroscopic pipeline')

    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()

    Verbose = args.v
    filenames = args.images  
    
    Prepare(filenames,Verbose)
    pass





