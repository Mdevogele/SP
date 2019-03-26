#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:57:17 2019

@author: maximedevogele
"""

import SP_Prepare2
import SP_Log
import SP_Bias
import SP_Flat

import glob

def Run(filenames): 
    
SP_Prepare2.Prepare(filenames)

Bias = glob('*Bias*')
SP_Bias.Create_Bias(Bias)






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
    parser.add_argument('-g',
                        help='Generic name of the offset and spectra location')
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
    Suffix = args.s
    Method = args.m
    Verbose = args.v
    Range = args.r
    filenames = args.images    
    Spec_loc = args.g
    
    Ran = []
    Ran.append(int(Range[0])-int(Range[1])/2)
    Ran.append(int(Range[0])+int(Range[1])/2)
    
    
    # call run_the_pipeline only on filenames
    Run(filenames)
    pass