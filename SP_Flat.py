#!/usr/bin/env python

""" SP_Flat - wrapper for creating flat(s) 
    v1.0: 2018-04-17, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import numpy as np

from astropy.io import fits

from itertools import groupby
from operator import itemgetter

import SP_diagnostics as diag
   

def Create_Flat(filenames,MasterName,Verbose,Bias,Series,Diagnostic):

     

    if Verbose:
        print('Beginning flat processing')
        print('Processing files:')
        print('index \t filename')
        for idx,elem in enumerate(filenames):
            print('{} \t {}'.format(idx+1,elem))
        print('Is using {} as master bias'.format(Bias))
        
    if Series == 'none':
        List_Idx_Series = [range(len(filenames))]
        File_Idx_Series = filenames
    
    if Series == 'index':
        index = []
        for elem in filenames:
            index.append(int(elem.split('_')[1]))
        File_Idx_Series = []
        List_Idx_Series = []
        for k, g in groupby(enumerate(index), lambda x: x[1]-x[0]):
            File_Idx_Series.append(map(itemgetter(1), g))
        for k, g in itertools.groupby(enumerate(index), lambda x: x[1]-x[0]):
            List_Idx_Series.append(map(itemgetter(0), g))  
    
        
        if Verbose:    
            print('SP_Flat detected {} serie(s) of flats based on files indexes'.format(len(File_Idx_Series)))
            for idx, elem in enumerate(File_Idx_Series):
                print('Serie {} contains indexes {}'.format(idx+1,elem))

    for idx, elem in enumerate(List_Idx_Series):
        
        FlatName = MasterName.replace('.fits','_' + str(idx + 1)) + '.fits'
        
        if Verbose:
            print('Creating the master flat {}'.format(FlatName))  


        FlatFlags = {
            'LogFile':'FlatLog.txt','RawPath':'','WriteFile': FlatName,
            'Verbose':False, 'Bias': './' + Bias, 'OverWrite' : True,'AddFits': False, 'IsGMOS': False
            }
        
        FileList = []
        for filesidx in elem:
            FileList.append(filenames[filesidx])
            
            
            
        SP.Create_Flat(FileList,**FlatFlags)

        if Diagnostic:
            diag.create_website('Flat_Log.html')
            diag.add_BiasSummary(filenames,FlatName,'Flat_Log.html')
            diag.add_FlatList(filenames,'Flat_Log.html')

        if Verbose:
            print('Master flat save to {}'.format(FlatName))
            hdulist = fits.open(FlatName)
            data = hdulist[0].data
            print('Statistics of the master flat')
            print('Mean: {} \t Median: {} \t std: {}'.format(np.mean(data), np.median(data), np.std(data)))
            print('End of flat processing')
            hdulist.close()
            

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master flats')

    parser.add_argument('-m', help='If there is several series of flat \n || Options: none: Only one serie \n || index: split according to the index of the files \n || target: split according to the target name in fits headers \n || pointing: split according to telescope pointing according to fits headers',
                        choices=['none','index','target','pointing'],
                        default = 'none')
    parser.add_argument('-b',
                        help='Name of the master bias to use \n || Can use None if no bias are available',
                        default='MasterBias.fits')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-o',
                        help='Prefix of the master flats files',
                        default='MasterFlat.fits')  
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
#    man_targetname = args.target
    Series = args.m
    Bias = args.b
    Verbose = args.v
    MasterName = args.o
    filenames = args.images    
    Diagnostic = args.d

    
    # call run_the_pipeline only on filenames
    Create_Flat(filenames,MasterName,Verbose,Bias,Series,Diagnostic)
    pass
