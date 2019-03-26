#!/usr/bin/env python

""" SP_CheckGrating - Script that check the grating for each file and correct it if necessary
    v1.0: 2018-03-15, michael.mommert@nau.edu
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

from astropy.io import fits

from SP_CheckInstrument import CheckInstrument



# setup logging
logging.basicConfig(filename = _SP_conf.log_filename,
                    level    = _SP_conf.log_level,
                    format   = _SP_conf.log_formatline,
                    datefmt  = _SP_conf.log_datefmt)


def CheckGrating(filenames):
    
    
    telescope, obsparam = CheckInstrument(filenames)

    ### read grating information from fits headers
    # check that they are the same for all images
    logging.info(('check for same grating for %d ' + \
                  'frames') % len(filenames))
    grating = []
    for idx, filename in enumerate(filenames):
        try:
            hdulist = fits.open(filename, ignore_missing_end=True)
        except IOError:
            logging.error('cannot open file %s' % filename)
            print('ERROR: cannot open file %s' % filename)
            filenames.pop(idx)
            continue

        header = hdulist[0].header
        grating.append(header[obsparam['grating']])
        print('%s \t %s' % (filename,header[obsparam['grating']]))

    if len(grating) == 0:
        raise KeyError('cannot identify filter; please update' + \
                       'setup/telescopes.py accordingly')

    if len(set(grating)) == 1:
        print('All the files possess the same grating')
    else:
        out = 0
        while out == 0:
            print('ERROR: multiple gratings used in dataset: %s' % str(set(grating)))
            text = raw_input("Do you want to modify the gratings of all the files ? : (['y'],'n']) ")
            if text == '' or text == 'y':
                grating = raw_input("What grating do you want to put in the header ? : ")
                out = 1
            else:
                if text == 'n':
                    print('The headers has not been changed.')
                    print('Exiting...')
                    sys.exit()
                else:
                    print("SP did not recognized you answer, please write 'y' or 'n' ")
                    text = raw_input("Do you want to modify the gratings of all the files ? : (['y'],'n']) ")
    
    if text == 'y' or '':
        for idx, filename in enumerate(filenames):
            try:
                hdulist = fits.open(filename, mode='update')
            except IOError:
                logging.error('cannot open file %s' % filename)
                print('ERROR: cannot open file %s' % filename)
                filenames.pop(idx)
                continue
            hdulist[0].header[obsparam['grating']] = grating
            hdulist.close()
        
        

    
    
    