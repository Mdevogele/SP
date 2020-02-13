#!/usr/bin/env python

""" SP_Asteroid - wrapper for asteroid(s) processing 
    v1.0: 2018-04-18, mdevogele@lowell.edu
"""

import argparse, shlex

def Parse_Asteroids(filenames,Verbose):
    
    Asteroid_List = []
    for idx, elem in enumerate(filenames):
        Asteroid_List.append(elem.split('_')[4])
    
    Asteroid_List = set(Asteroid_List)   
    
    if Verbose:
        print('{} asteroids have been found'.format(len(Asteroid_List)))
        for elem in Asteroid_List:
            print(elem)

    return Asteroid_List


def Get_IndivualList(filenames,Ast):
    
    Asteroid_List = []
    for idx, elem in enumerate(filenames):
        if elem.split('_')[4] == Ast:
            Asteroid_List.append(elem)

    return Asteroid_List

def Get_IndivualList(filenames,Ast):
    
    Asteroid_List = []
    for idx, elem in enumerate(filenames):
        if elem.split('_')[4] == Ast:
            Asteroid_List.append(elem)

    return Asteroid_List




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master flats')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
#    parser.add_argument('-target', help='primary targetname override',
#                        default=None)
    parser.add_argument('-d', help='add flats information to diagnostic.htlm file',
                        action="store_true")
    parser.add_argument('-s', help='If there is several series of flat \n || Options: none: Only one serie \n || index: split according to the index of the files \n || target: split according to the target name in fits headers \n || pointing: split according to telescope pointing according to fits headers',
                        choices=['none','index','target','pointing'],
                        default = 'None')
    parser.add_argument('-b',
                        help='Name of the master bias to use \n || Can use None if no bias are available',
                        default='MasterBias.fits')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-o',
                        help='Prefix of the master flats files',
                        default='MasterFlat.fits')    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
#    man_targetname = args.target
    Diagnostic = args.d
    Series = args.s
    Bias = args.b
    Verbose = args.v
    MasterName = args.o
    filenames = args.images    
    
    # call run_the_pipeline only on filenames
    Asteroid_List = Parse_Asteroids(filenames,Verbose)
    
    for idx, elem in enumerate(Asteroid_List):
        FileList = Get_IndivualList(filenames,elem)
        print(FileList)

    
    pass
