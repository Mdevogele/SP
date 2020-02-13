#!/usr/bin/env python

""" SP_run - wrapper for automated data analysis
    v1.0: 2018-03-30, mdevogele@lowell.edu
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


import argparse, shlex

from astropy.io import fits

from SP_Prepare import Prepare as Prepare
from SP_CheckInstrument import CheckInstrument
import SP_Deveny_Toolbox as SP_Dev_Tool
import SP_Toolbox as Toolbox
from SP_Bias import Create_Bias
from SP_Flat import Create_Flat


import SP_diagnostics as diag


# setup logging
logging.basicConfig(filename = _SP_conf.log_filename,
                    level    = _SP_conf.log_level,
                    format   = _SP_conf.log_formatline,
                    datefmt  = _SP_conf.log_datefmt)


def run_the_pipel(filenames, Verbose):
    
    """
    wrapper to run the photometry pipeline
    """

    # reset diagnostics for this data set
    _SP_conf.dataroot, _SP_conf.diagroot, \
    _SP_conf.index_filename, _SP_conf.reg_filename, _SP_conf.cal_filename, \
    _SP_conf.res_filename = _SP_conf.setup_diagnostics()
    

    ## diag.create_index(filenames,_SP_conf.dataroot) # commented for development uncomment for dealing with new unprepared data 

    ### read telescope information from fits headers
    # check that they are the same for all images
    logging.info('##### new spectroscopic process in %s #####' % _SP_conf.dataroot)
    logging.info(('check for same telescope/instrument for %d ' + \
                  'frames') % len(filenames))

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

    if len(grating) == 0:
        raise KeyError('cannot identify filter; please update' + \
                       'setup/telescopes.py accordingly')

    if len(set(grating)) > 1:
        print('ERROR: multiple gratings used in dataset: %s' % str(set(grating)))
        print('Check the grating of each files by running SP_CheckGrating')
#        logging.error('multiple grating used in dataset: %s' %
#                      str(set(grating)))
#        for i in range(len(filenames)):
#            logging.error('%s %s' % (filenames[i], grating[i]))
        sys.exit()    
    
    """ run SP_prepare """
    
    ### Prepare(filenames) # commented for development uncomment for dealing with new unprepared data 
    
    _SP_conf.filenames = filenames # for development purpose, delete if Prepare(filenames) uncommented 

###    Toolbox.Get_BiasList()
     
    """ Add the bias file list to the diagnostic.html file """
    
###    diag.add_BiasList()
    
    """ Creation of the Master Bias, SP_run assumed that there is only one Master Bias to be created. One can use SP_Bias manually otherwise """
    
###    MasterBiasName = 'MasterBias.fits'
###    Create_Bias(_SP_conf.Bias_filenames,MasterBiasName,Verbose)
    
###    diag.add_BiasSummary(MasterBiasName)
    
###    Toolbox.Get_FlatList()
        
    """ Creation of the Master Flat """

###    diag.add_FlatList()
    
###    MasterFlatName = 'MasterFlat.fits'
###    Create_Flat(_SP_conf.Flat_filenames,MasterFlatName,Verbose,MasterBiasName,'index',True)
    
    
    """ processing of the target spectra """
    
    Toolbox.Get_AsteroidList()
    
    print(_SP_conf.Asteroid_filenames)
    
    
    
    
    
#    print('run spectroscopic pipeline on %d %s frames' % \
#    (len(filenames), telescope))
    
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='automated WCS registration')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
#    parser.add_argument('-target', help='primary targetname override',
#                        default=None)
#    parser.add_argument('-filter', help='filter name override',
#                        default=None)
#    parser.add_argument('-fixed_aprad', help='fixed aperture radius (px)',
#                        default=0)
#    parser.add_argument('-source_tolerance',
#                        help='tolerance on source properties for registration',
#                        choices=['none', 'low', 'medium', 'high'],
#                        default='high')
#    parser.add_argument('-solar',
#                        help='restrict to solar-color stars',
#                        action="store_true", default=False)
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")       
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
#    man_targetname = args.target
#    man_filtername = args.filter
#    fixed_aprad = float(args.fixed_aprad)
#    source_tolerance = args.source_tolerance
#    solar = args.solar
    Verbose = args.v
    filenames = args.images    
    
    # call run_the_pipeline only on filenames
    run_the_pipel(filenames,Verbose)
    pass
    
    
    