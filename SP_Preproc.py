#!/usr/bin/env python

""" SP_Preproc - Apply the flat and bias to science data
    v1.0: 2018-04-19, mdevogele@lowell.edu
        
"""
import argparse, shlex

import SP_Toolbox as SP
import logging

import SP_diagnostics as diag


def Preproc(filenames,MasterBias,MasterFlat,Verbose,Suffix,Diagnostic):

    logging.info('****************************************')
    logging.info('****** Start of SP_Preproc script ******')
    logging.info('****************************************')
    
    
    
    Object_out = []
    Object_in = []
    compteur = 0
    for elem in filenames:
        compteur += 1
        Object_out.append(elem.replace('.fits','') + '_' + Suffix)
        Object_in.append(elem[:-5])


    print(Object_in)

    ReducFlags = {'Bias': MasterBias, 
                      'Flat':MasterFlat,'OutFilelist': Object_out,'OverWrite':True, 'Verbose' : True, 'IsGMOS': False, 'RawPath':''}   
       
    
    SP.Reduce(Object_in,**ReducFlags)


    if Diagnostic: 
        diag.create_website('Pre-processing_Log.html')
        diag.add_PreProc(Object_out,MasterBias,MasterFlat,'Pre-processing_Log.html')


    logging.info('**************************************')
    logging.info('****** End of SP_Preproc script ******')
    logging.info('**************************************')        



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master flats')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
    parser.add_argument('-s',
                        help='Suffix to add to processed files',
                        default='Procc')
    parser.add_argument('-b',
                        help='Name of the master bias to use \n || Can use None if no bias are available',
                        default='MasterBias.fits')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-f',
                        help='Name of the master flat to use \n || Can use None if no bias are available',
                        default='MasterFlat.fits')    
    
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  
    
    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
    Suffix = args.s
    MasterBias = args.b
    Verbose = args.v
    MasterFlat = args.f
    filenames = args.images    
    Diagnostic = args.d
    
    print(filenames)
    
    Preproc(filenames,MasterBias,MasterFlat,Verbose,Suffix,Diagnostic)
    pass