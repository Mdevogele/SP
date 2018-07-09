#!/usr/bin/env

""" SP_CheckInstrument - Script that check the instrument for each file
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



# setup logging
logging.basicConfig(filename = _SP_conf.log_filename,
                    level    = _SP_conf.log_level,
                    format   = _SP_conf.log_formatline,
                    datefmt  = _SP_conf.log_datefmt)


def CheckInstrument(filenames):
    instruments = []
    for idx, filename in enumerate(filenames):
        try:
            hdulist = fits.open(filename, ignore_missing_end=True)
        except IOError:
            logging.error('cannot open file %s' % filename)
            print('ERROR: cannot open file %s' % filename)
            #filenames.pop(idx)
            continue

        header = hdulist[0].header
        for key in _SP_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])
                break

    if len(filenames) == 0:
        raise IOError('cannot find any data...')

    if len(instruments) == 0:
        raise KeyError('cannot identify telescope/instrument; please update' + \
                       '_SP_conf.instrument_keys accordingly')


    # check if there is only one unique instrument
    if len(set(instruments)) > 1:
        print('ERROR: multiple instruments used in dataset: %s' % \
            str(set(instruemnts)))
        logging.error('multiple instruments used in dataset: %s' %
                      str(set(instruments)))
        for i in range(len(filenames)):
            logging.error('%s %s' % (filenames[i], instruments[i]))
        sys.exit()

    telescope = _SP_conf.instrument_identifiers[instruments[0]]
    obsparam = _SP_conf.telescope_parameters[telescope]
    logging.info('%d %s frames identified' % (len(filenames), telescope))
    
    return telescope, obsparam