#!/usr/bin/env

""" SP_ObsTypeCheck - Script that check the observation type for each file and correct it if necessary
    v2.0: 2018-12-14, mdevogele@lowell.edu
"""
from __future__ import print_function

# Spectroscopic Pipeline
# Copyright (C) 2018  Maxime Devogele, mdevogele@lowell.edu

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

import _SP_conf
import sys
import logging
import numpy as np

from astropy.io import fits

from SP_CheckInstrument import CheckInstrument

# setup logging
logging.basicConfig(filename = _pp_conf.log_filename,
                    level    = _pp_conf.log_level,
                    format   = _pp_conf.log_formatline,
                    datefmt  = _pp_conf.log_datefmt)


def CheckObsType(filenames):
    
    telescope, obsparam = CheckInstrument(filenames)
    
    CollFoc = [];
    for idx, filename in enumerate(filenames):
        try:
            hdulist = fits.open(filename, mode='update' ,ignore_missing_end=True)
        except IOError:
            logging.error('cannot open file %s' % filename)
            print('ERROR: cannot open file %s' % filename)
            filenames.pop(idx)
            continue

        data = hdulist[0].data
        Med = np.median(data[10:-10,10:-10])
        std0 = np.std(np.median(data[10:-10,10:-10],axis=0))
        std1 = np.std(np.median(data[10:-10,10:-10],axis=1))
        
        if Med > 10000:
            print('%s changed to %s' % (hdulist[0].header[obsparam['obstype']], 'FLAT'))
            hdulist[0].header[obsparam['obstype']] = 'FLAT'
        else:
            if std0 + std1 < 10:
                print('%s changed to %s' % (hdulist[0].header[obsparam['obstype']], 'BIAS'))
                hdulist[0].header[obsparam['obstype']] = 'BIAS'
            else:
                if std0/std1 > 100:
                    if telescope == 'DEVENY':
                        CollFoc.append((int(filename.split('.')[1]),hdulist[0].header['COLLFOC'],idx))                                
                        print('%s changed to %s' % (hdulist[0].header[obsparam['obstype']], 'ARCS/FOCUS'))
                        hdulist[0].header[obsparam['obstype']] = 'ARCS/FOCUS'                                   
                    else:
                        print('%s changed to %s' % (hdulist[0].header[obsparam['obstype']], 'ARCS'))
                        hdulist[0].header[obsparam['obstype']] = 'ARCS'
                else:
                    print('%s changed to %s' % (hdulist[0].header[obsparam['obstype']], 'OBJECT'))
                    hdulist[0].header[obsparam['obstype']] = 'OBJECT'
        
        hdulist.close()

                 
    Cons = toolbox.Get_Consecutive(np.array(CollFoc)[:,0].astype(int))
    Focus = []
    Series = []
    for idx, elem in enumerate(Cons):
        for idx2, elem2 in enumerate(elem):
            Focus.append(CollFoc[Cons[idx][idx2][0]][1])
        if len(set(Focus)) == 1:
            Series.append('ARCS')
        else:
            Series.append('FOCUS')  
        print(Focus)
        Focus = []                 
            
    for idx, elem in enumerate(Cons):
        for idx2, elem2 in enumerate(elem):    
            try:
                hdulist = fits.open(filenames[CollFoc[Cons[idx][idx2][0]][2]], mode='update' ,ignore_missing_end=True)
            except IOError:
                logging.error('cannot open file %s' % filename)
                print('ERROR: cannot open file %s' % filename)
                filenames.pop(idx)
                continue
            if Series[idx] == 'ARCS':
                print('%s changed to %s' % (hdulist[0].header[obsparam['obstype']], 'ARCS'))
                hdulist[0].header[obsparam['obstype']] = 'ARCS'
            else:
                print('%s changed to %s' % (hdulist[0].header[obsparam['obstype']], 'FOCUS'))
                hdulist[0].header[obsparam['obstype']] = 'FOCUS'
        
            hdulist.close()
                  
        
       
        
        
        
                    
                    
                    
        
    
    