#!/usr/bin/env python

""" SP_MaskSpectrum - Mask a spectrum close to the asteroid spectrum  
    v1.0: 2018-06-20, mdevogele@lowell.edu
  
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

import random
# setup logging
logging.basicConfig(filename = _SP_conf.log_filename,
                    level    = _SP_conf.log_level,
                    format   = _SP_conf.log_formatline,
                    datefmt  = _SP_conf.log_datefmt)


def Spec_Mask(filenames,loc,outfile,repl):
    
    
    hdulist = fits.open(filenames[0])
    
    data = hdulist[0].data
    
    data[int(loc[0])-int(loc[1]):int(loc[0])+int(loc[1]),:] = data[int(repl),:] 
    
    hdulist[0].data = data
    
    hdulist.writeto(outfile)

 
    return None   


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Mask a spectrum too close to the target')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
#    parser.add_argument('-target', help='primary targetname override',
#                        default=None)
#    parser.add_argument('-m', help='add flats information to diagnostic.htlm file',
#                        choices=['auto','range'])
#    parser.add_argument('-s', help='If there is several series of flat \n || Options: none: Only one serie \n || index: split according to the index of the files \n || target: split according to the target name in fits headers \n || pointing: split according to telescope pointing according to fits headers',
#                        choices=['none','index','target','pointing'],
#                        default = 'None')
    parser.add_argument('-o',
                        help='Name of the outputfile')
    parser.add_argument('-l',
                        help='Possition of the spectrum and number of pixels to mask',
                        nargs=2)
    parser.add_argument('-r',
                        help='Value to replace the pixel for',
                        default = 100)
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
    loc = args.l
    repl = args.r
    outfile = args.o
    filenames = args.images  
    
    print(filenames)

    
    
    # call run_the_pipeline only on filenames
    Spec_Mask(filenames,loc,outfile,repl)
    pass





