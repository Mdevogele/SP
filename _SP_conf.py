#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
configuration file for spectroscopy pipeline

2018-03-13, mdevogele@lowell.edu
"""
from __future__ import print_function

import os
import sys
import logging
import warnings

try:
#    from astropy import wcs
    from astropy.io import fits
except ImportError:
    print('Module astropy not found. Please install with: pip install astropy')
    sys.exit()

try:
    import numpy as np
except ImportError:
    print('Module numpy not found. Please install with: pip install numpy')
    sys.exit()
    

# potential FITS header keywords for looking up the instrument
# any unique header keyword works as a potential identifier
instrument_keys = ['INSTRUME', 'LCAMMOD']

# suppress runtime and astropy warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# only import if Python3 is used
if sys.version_info > (3, 0):
    from past.builtins import execfile


# read spectroscopic pipeline root path from environment variable
rootpath = os.environ.get('SPECPIPEDIR')
if rootpath is None:
    print('ERROR: SPECPIPEDIR variable has not been set')
    sys.exit(0)


# potential FITS header keywords for looking up the instrument
# any unique header keyword works as a potential identifier
instrument_keys = ['INSTRUME', 'LCAMMOD']

execfile(rootpath+'/setup/telescopes.py')
