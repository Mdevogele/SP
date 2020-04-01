#!/usr/bin/env python

""" SP_CosmCorr - Correct for Cosmic rays
    v1.0: 2018-06-15, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import numpy as np

from SP_CheckInstrument import CheckInstrument

import SP_diagnostics as diag


from astropy.io import fits

import cosmic

from astroscrappy import detect_cosmics



def Cosmic(filename,Diagnostic):
    
    LACOSMIC = False
    
    Out_Name = []
    for idx, elem in enumerate(filename): 
        
        # load the fits files
        print(elem)
        hdulist = fits.open(elem)
        # Extract data from fits files
        data=hdulist[0].data
        
        #Start cosmic correction
        if LACOSMIC:
            # Does not work for python 3 as the cosmic code has been developped for python 2 only
            # However, for some reason the detect_cosmic code happen to not work for some images without throwing any errors ?? 
            # Keep the cosmic option here for such rare cases 
            import cosmic
            c = cosmic.cosmicsimage(data,satlevel=65000)
            c.run(maxiter = 1)
            hdulist[0].data = c.cleanarray
        else:
            crmask, c = detect_cosmics(data, inmask=None, satlevel=65000)
            hdulist[0].data = c

        # Create the output name. Remove _Procc and add _CosmCorr to the name
        Out_Name.append(elem.replace('.fits','').replace('_Procc','') + '_' + '_CosmCorr' + '.fits')
        print(Out_Name[idx])
        hdulist.writeto(Out_Name[idx])
        
    if Diagnostic: 
        diag.create_website('Cosmic-Correction_Log.html')
        diag.add_CosmCorr(Out_Name,'Cosmic-Correction_Log.html')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Create a log file')

    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  

    args = parser.parse_args()
    filenames = args.images    
    Diagnostic = args.d
    
    # call run_the_pipeline only on filenames
    Cosmic(filenames,Diagnostic)
    pass
