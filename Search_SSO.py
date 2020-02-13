#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 23:06:13 2019

@author: maximedevogele
"""

from astropy.io import fits
import astroquery

import requests
from astropy.io import ascii

import argparse, shlex


def Search(filenames):
    
    for elem in filenames:

        hdulist = fits.open(elem)
        
        RA = hdulist[0].header['RA']
        DEC = hdulist[0].header['DEC']
        print(RA)
        print(DEC)
        
        ra_deg = (float(RA[0:2])+float(RA[3:5])/60+float(RA[6:8])/3600)*15
        if '-' in DEC[0]:
            dec_deg = float(DEC[0:3])-float(DEC[4:6])/60-float(DEC[7:9])/3600
        else:
            dec_deg = float(DEC[0:2])+float(DEC[3:5])/60+float(DEC[6:8])/3600
        
        
        Time = hdulist[0].header['DATE-OBS']
        
        server = 'http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php'
        
        
        r = requests.get(server, 
                             params= {'RA': ra_deg, 'DEC': dec_deg, 
                                      'SR': 0.12, 'EPOCH': Time,
                                      '-output': 'object',
                                      '-mime': 'text'},
                             timeout=180)
 
        if 'No solar system object was found in the requested FOV' in r.text:
            print('No Solar System object found')
            f = open(elem.replace('.fits','_SSO.txt'),'w')
            f.write('No Solar System object found')
            f.close()
        else:

       
        
            results = ascii.read(r.text, delimiter='|',
                                     names=('number', 'name', 'ra', 'dec', 'type',
                                            'V', 'posunc', 'd'))
            
            ascii.write(results,elem.replace('.fits','_SSO.txt'))
            
            print(results)


if __name__ == '__main__':

    # command line arguments
    parser = argparse.ArgumentParser(description='Search for SSO in the images')

    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    filenames = args.images


    Search(filenames)
