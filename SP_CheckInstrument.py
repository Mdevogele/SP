#!/usr/bin/env

""" SP_CheckInstrument - Script that check the instrument for each file
    v1.0: 2018-03-15, mdevogele@lowell.edu
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

import SP_Toolbox as SP


def CheckInstrument(filenames):
    instruments = []
    for idx, filename in enumerate(filenames):
        try:
            hdulist = fits.open(filename, ignore_missing_end=True)
        except IOError:
            logging.error('cannot open file %s' % filename)
            print('ERROR: cannot open file %s' % filename)
            continue

        header = hdulist[0].header
        for key in _SP_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])
                break

    if len(instruments) == 0:
        raise KeyError('cannot identify telescope/instrument; please update' + \
                       '_SP_conf.instrument_keys accordingly')


    # check if there is only one unique instrument
    if len(set(instruments)) > 1:
        print('ERROR: multiple instruments used in dataset: %s' % \
            str(set(instruments)))
        logging.error('multiple instruments used in dataset: %s' %
                      str(set(instruments)))
        for i in range(len(filenames)):
            logging.error('%s %s' % (filenames[i], instruments[i]))
        sys.exit()

    telescope = _SP_conf.instrument_identifiers[instruments[0]]
    obsparam = _SP_conf.telescope_parameters[telescope]    
    
    
    # for the different detectors of GMOS North and Souths
    # Definition of the location of the CHIP GAPS and detector pixscale
    if telescope == 'GMOSN':
        if hdulist[0].header[obsparam['detector']] == 'GMOS + e2v DD CCD42-90':
            obsparam['pixscale'] = 0.07288
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 2:
                obsparam['ship_gap'] = ([1023,1046],[2068,2091])
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 4:
                obsparam['ship_gap'] = ([511,523],[1033,1045])
        if hdulist[0].header[obsparam['detector']] == 'GMOS-N + Hamamatsu':
            obsparam['pixscale'] = 0.0807
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 2:    #to be defined
                obsparam['ship_gap'] = ([0,0],[0,0])
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 4:    #to be defined
                obsparam['ship_gap'] = ([0,0],[0,0])
        if obsparam['pixscale'] == 0:
                    print('The detector: ' + hdulist[0].header[obsparam['detector']] + 'is not found it the database, use of default 0.08 arcsec/pixel default prixel scale')

    if telescope == 'GMOSS':
        if hdulist[0].header[obsparam['detector']] == 'GMOS + Blue1 + new CCD1':
            obsparam['pixscale'] = 0.073
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 2:    #to be defined
                obsparam['ship_gap'] = ([0,0],[0,0])
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 4:    #to be defined
                obsparam['ship_gap'] = ([0,0],[0,0])
        if hdulist[0].header[obsparam['detector']] == 'GMOS + Hamamatsu_new':  
            obsparam['pixscale'] = 0.08
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 2:    #to be defined
                obsparam['ship_gap'] = ([0,0],[0,0])
            if SP.get_binning(hdulist[1].header,obsparam)[0] == 4:    #to be defined
                obsparam['ship_gap'] = ([0,0],[0,0])
        
        if obsparam['pixscale'] == 0:
                    print('The detector: ' + hdulist[0].header[obsparam['detector']] + 'is not found it the database, use of default 0.08 arcsec/pixel default prixel scale')
        
    
    
    logging.info('%d %s frames identified' % (len(filenames), telescope))
    
    return telescope, obsparam