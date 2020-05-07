"""
Spectroscopic Pipeline Configuation File
2018-03-15, mdevogele@lowell.edu
"""

# Spectroscopic Pipiline
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

##### telescope/instrument configurations


# NOT La Palma
NOT_param = {
    'telescope_instrument' : 'ALFOSC_FASU', # telescope/instrument name
    'telescope_keyword'    : 'NOT',      # telescope/instrument keyword
    'secpix'               : 0.34, # pixel size (arcsec)
                                            # before binning

    # image orientation preferences
    'flipx'                : True,      # Is the wavelength increase with increasing X values ?

    # instrument-specific FITS header keywords
    'binning'              : ('CCDSUM#blank0', 'CCDSUM#blank1'),
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword
    'filter'               : 'FILTREAR',  # filter keyword
    'filter_translations'  : {},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword
    
    'grating'              : 'ALFLTNM', # grating used
    'obstype'              : 'OBSTYPE', # get the type of observation (bias, flat, arcs, focus)
    
    'focus'                : 'COLL_FOC', # Collimator focus header
    
    'crop'                 : True,
    'Y_crop'               : [800,-100],
    'X_crop'               : False,
    
    'transpose'            : True, # Is the spectrum vertical or horizontal? True if vertical 
    
    'index_data'           : 1,
    
    # Limits for BIAS, FLAT, ARCS, OBJECT auto detection
    
    'flat_median'          : 10000,   # if the median of all pixels is > 10000 files considered as a flat
    'bias_std'             : 10,      # if the sum of the std of the median of each axis > 10 files considered as Bias
    'arc_object_std'       : 100      # if the ratio of std(median(axis=0)) over std(median(axis=1)) > 100 files considered as Arcs, Object otherwise     

}



# GOODMAN@SOAR Chile
Soar_param = {
    'telescope_instrument' : 'Goodman Spectro', # telescope/instrument name
    'telescope_keyword'    : 'SOAR',      # telescope/instrument keyword
    'secpix'               : 0.34, # pixel size (arcsec)
                                            # before binning

    # image orientation preferences
    'flipx'                : True,      # Is the wavelength increase with increasing X values ?

    # instrument-specific FITS header keywords
    'binning'              : ('CCDSUM#blank0', 'CCDSUM#blank1'),
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'OBSRA',  # telescope pointing, RA
    'dec'                  : 'OBSDEC', # telescope pointin, Dec
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword
    'filter'               : 'FILTREAR',  # filter keyword
    'filter_translations'  : {},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword
    
    'posangle'             : 'POSANGLE', 
    'grating'              : 'GRATING', # grating used
    'grat_ang'             : 'GRT_ANG', # grating angle
    'obstype'              : 'OBSTYPE', # get the type of observation (bias, flat, arcs, focus)
    
    'focus'                : 'COLL_FOC',
    
    'crop'                 : True,
    'Y_crop'               : [20,-20],
    'X_crop'               : False,
    
    'lamp1'                : 'LAMP_HGA',
    'lamp1_name'           : 'HgAr',
    'lamp2'                : 'LAMP_NE',
    'lamp2_name'           : 'Ne',
    'lamp3'                : 'LAMP_AR',
    'lamp3_name'           : 'Ar',
    'lamp4'                : 'LAMP_FE',
    'lamp4_name'           : 'Fe',    
    
    'transpose'            : False, # Is the spectrum vertical or horizontal? True if vertical

    'index_data'           : 0,    
    
    # Limits for BIAS, FLAT, ARCS, OBJECT auto detection
    
    'flat_median'          : 10000,   # if the median of all pixels is > 10000 files considered as a flat
    'bias_std'             : 10,      # if the sum of the std of the median of each axis > 10 files considered as Bias
    'arc_object_std'       : 100      # if the ratio of std(median(axis=0)) over std(median(axis=1)) > 100 files considered as Arcs, Object otherwise     

}




# Deveny@DCT Lowell Observatory
Deveny_param = {
    'telescope_instrument' : 'Deveny/DCT', # telescope/instrument name
    'telescope_keyword'    : 'DCT',      # telescope/instrument keyword
    'secpix'               : 0.34, # pixel size (arcsec)
                                            # before binning

    # image orientation preferences
    'flipx'                : True,      # Is the wavelength increase with increasing X values ?

    # instrument-specific FITS header keywords
    'binning'              : ('DETXBIN', 'DETYBIN'),
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec
    'radec_separator'      : 'XXX',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'SCITARG',  # object name keyword
    'filter'               : 'ALFLTNM',  # filter keyword
    'filter_translations'  : {'Open': None},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword
    
    'posangle'             : 'IPA', 
    'grat_ang'             : 'GRANGLE', # grating angle 
    'grating'              : 'GRATING', # grating used
    'obstype'              : 'IMAGETYP', # get the type of observation (bias, flat, arcs, focus)
    
    'focus'                : 'COLLFOC',

    'crop'                 : False,
    'X_crop'               : False,
    'Y_crop'               : False,
    
    'lamp1'                : 'LAMP1ON',
    'lamp1_name'           : 'Cd',
    'lamp2'                : 'LAMP2ON',
    'lamp2_name'           : 'Ar',
    'lamp3'                : 'LAMP3ON',
    'lamp3_name'           : 'Ne',
    'lamp4'                : 'LAMP4ON',
    'lamp4_name'           : 'Hg',    
    
    
    'transpose'            : False, # Is the spectrum vertical or horizontal? True if vertical
    
    'data_index'           : 0, # index of the hdu cart containing the actual data

    
    # Limits for BIAS, FLAT, ARCS, OBJECT auto detection
    
    'flat_median'          : 10000,   # if the median of all pixels is > 10000 files considered as a flat
    'bias_std'             : 10,      # if the sum of the std of the median of each axis > 10 files considered as Bias
    'arc_object_std'       : 100      # if the ratio of std(median(axis=0)) over std(median(axis=1)) > 100 files considered as Arcs, Object otherwise     

}


##### THE GMOS STILL NEED TO BE MODIFIED

# GMOS-N@Gemini North 
GMOSN_param = {
    'telescope_instrument' : 'GMOS-N', # telescope/instrument name
    'telescope_keyword'    : 'GMOS-N',  # telescope/instrument keyword
    'observatory_code'     : '568',         # MPC observatory code
    'secpix'               : 0.080, # pixel size (arcsec)
                                            # before binning
    # image orientation preferences
    'flipx'                : True,      # Is the wavelength increase with increasing X values ?

    # instrument-specific FITS header keywords
    'binning'              : ('CCDSUM#blank0', 'CCDSUM#blank1'),
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec
    'radec_separator'      : 'XXX',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'OBSEPOCH', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword
    'filter'               : 'FILTREAR',  # filter keyword
    'filter_translations'  : {},
                             # filtername translation dictionary
    'exptime'              : 'DARKTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword
    
    'grating'              : 'GRATING', # grating used
    'obstype'              : 'OBSTYPE', # get the type of observation (bias, flat, arcs, focus)
    'obsclass'             : 'OBSCLASS', # get the type of observation (bias, flat, arcs, focus)
    
    'focus'                : 'DTAZME',
    
    'tilt'                 : 'WAVELENG',
    
    'pixscale'             : 0, # Actualy defined in Check_Instrument based on the detector 
    'detector'             : 'DETECTOR',
    'ship_gap'             : ([0,0],[0,0]), # Actualy defined in Check_Instrument based on the detector


    'crop'                 : False,
    'X_crop'               : False,
    'Y_crop'               : False,
    
    'transpose'            : False, # Is the spectrum vertical or horizontal? True if vertical

    'index_data'           : 0, # After specific opening from 
    
    # Limits for BIAS, FLAT, ARCS, OBJECT auto detection
    
    'flat_median'          : 10000,   # if the median of all pixels is > 10000 files considered as a flat
    'bias_std'             : 10,      # if the sum of the std of the median of each axis > 10 files considered as Bias
    'arc_object_std'       : 100      # if the ratio of std(median(axis=0)) over std(median(axis=1)) > 100 files considered as Arcs, Object otherwise     
}





# GMOS-S@Gemini South 
GMOSS_param = {
    'telescope_instrument' : 'GMOS-S', # telescope/instrument name
    'telescope_keyword'    : 'GMOS-S',  # telescope/instrument keyword
    'observatory_code'     : 'I11',         # MPC observatory code
    'secpix'               : 0.080, # pixel size (arcsec)
                                            # before binning
    # image orientation preferences
    'flipx'                : True,      # Is the wavelength increase with increasing X values ?

    # instrument-specific FITS header keywords
    'binning'              : ('CCDSUM#blank0', 'CCDSUM#blank1'),
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec
    'radec_separator'      : 'XXX',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'OBSEPOCH', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword
    'filter'               : 'FILTREAR',  # filter keyword
    'filter_translations'  : {},
                             # filtername translation dictionary
    'exptime'              : 'DARKTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword
    
    'grating'              : 'GRATING', # grating used
    'obstype'              : 'OBSTYPE', # get the type of observation (bias, flat, arcs, focus)
    'obsclass'             : 'OBSCLASS', # get the type of observation (bias, flat, arcs, focus)
    
    'focus'                : 'DTAZME',
    
    'tilt'                 : 'WAVELENG',
    
    'pixscale'             : 0, # Actualy defined in Check_Instrument based on the detector 
    'detector'             : 'DETECTOR',
    
    'ship_gap'             : ([0,0],[0,0]), # Actualy defined in Check_Instrument based on the detector

    'crop'                 : False,
    'X_crop'               : False,
    'Y_crop'               : False,
    
    'transpose'            : False, # Is the spectrum vertical or horizontal? True if vertical

    'index_data'           : False, # Special treatment for GMOS

    
    # Limits for BIAS, FLAT, ARCS, OBJECT auto detection
    
    'flat_median'          : 10000,   # if the median of all pixels is > 10000 files considered as a flat
    'bias_std'             : 10,      # if the sum of the std of the median of each axis > 10 files considered as Bias
    'arc_object_std'       : 100      # if the ratio of std(median(axis=0)) over std(median(axis=1)) > 100 files considered as Arcs, Object otherwise     

}


##### access functions for telescope configurations


implemented_telescopes = ['SOAR','DEVENY', 'GMOS-S','GMOS-N','NOT']

# translate INSTRUME (or others, see _pp_conf.py) header keyword into
# PP telescope keyword
instrument_identifiers = {'ALFOSC_FASU':        'NOT',
                          'Goodman Spectro':        'SOAR',
                          'Deveny':        'DEVENY',
                          'GMOS-S':               'GMOSS',
                          'GMOS-N':              'GMOSN'
}

# translate telescope keyword into parameter set defined here
telescope_parameters = {'NOT' :       NOT_param,
                         'SOAR' :       Soar_param,
                        'DEVENY' :       Deveny_param,
                        'GMOSS':        GMOSS_param,
                        'GMOSN':        GMOSN_param
}


#### append mytelescopes.py, if available
#
# mytelescopes.py allows you to setup your own telescope; that file is
# not part of the github repository, hence it will not be affected by
# pulls and pushes
#
# an example mytelescopes.py file is available here:
# http://134.114.60.45/photometrypipeline/mytelescopes.py
# more information are available on the PP documentation website:
# http://mommermi.github.io/pp/install.html

try:
    execfile(rootpath+'/setup/mytelescopes.py')
except IOError:
    pass
