#!/usr/bin/env python

""" SP_Extract - Extract spectrum
    v1.0: 2018-04-19, mdevogele@lowell.edu
  
"""

import argparse, shlex

import SP_Toolbox as SP
import numpy as np

import simplejson

from astropy.io import fits
import operator

import SP_diagnostics as diag


from SP_CheckInstrument import CheckInstrument


def Extract_Spectrum(filename,Verbose,Spec_loc,Diagnostic,Spec_FWHM): 
    
    
    Out_Files = []
    Out_Spec = []
    telescope, obsparam = CheckInstrument([filename[0]])    
    DetecFlags = {'Instrument':telescope}
    if telescope == 'DEVENY' or telescope == 'SOAR':
        for elem in filename:
            hdulist = fits.open(elem)
            data = hdulist[0].data
            if telescope == 'DEVENY':
                DetecFlags = {'Instrument':'Deveny'}  
            if telescope == 'SOAR':
                DetecFlags = {'Instrument':'Soar'}  

            Center = SP.Detect_Spectra(data,Bin=2,**DetecFlags)
            Start = (1415,Center)
            
            Trace, bkg, MASK1 = SP.Fit_Trace(data,Start,Range = 15, SClip = True,**DetecFlags)

            TR=[]
            for idx,elem2 in enumerate(Trace):
                if elem2 > 50 and elem2 < np.size(data,0)-50:
                    TR.append(data[int(elem2)-45:int(elem2)+45,idx])
                else: 
                    TR.append(data[0:90,idx])
            TR = np.array(TR)
            hdulist[0].data = np.transpose(TR)
            hdulist.writeto(elem.replace('.fits','') + 'Trace.fits' )
            Out_Files.append(elem.replace('.fits','') + 'Trace.fits')

            SSpec = []
            for FW in range(20):
                Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=FW,Mask = MASK1,**DetecFlags)
                SSpec.append(Spec1)
                
#            max_index, max_value = max(enumerate(np.nanmedian(SSpec,axis=1)/np.nanstd(SSpec,axis=1)), key=operator.itemgetter(1))
            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=Spec_FWHM,Mask = MASK1,**DetecFlags)
            
#            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=50,Mask = MASK1,**DetecFlags)
            Spec1N = Spec1/np.abs(np.nanmedian(Spec1[1500:1600]))
            f = open(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt','w')
            Out_Spec.append(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt')
            for item in Spec1N:
                f.write("%s\n" % item)
            
            f.close()
    if telescope == 'GMOSS' or telescope == 'GMOSN':
        OffFile = Spec_loc + '_Offset.txt'
        with open(OffFile,'r') as f:
            Offset = simplejson.load(f)
        SpecFile = Spec_loc + '_Spec_loc.txt'    
        with open(SpecFile,'r') as f:
            Spec_Loc = simplejson.load(f)
        
        Offset = np.sort(np.array(Offset).astype(int))
        Spec_Loc = np.sort(Spec_Loc)
        for elem in filename:
            hdulist = fits.open(elem)
            Gratting = hdulist[0].header[obsparam['grating']]
            Binning = hdulist[1].header['CCDSUM'][0]
            Detector = hdulist[0].header['DETECTOR']
            
            DetecFlags = {'Instrument':telescope,'Gratting':Gratting,'Binning':Binning,'Detector':Detector}
            
            Off = int(hdulist[0].header['YOFFSET'])
            idx = np.where(Offset == Off)[0][0]
            if telescope == 'GMOSS':
                Spec = Spec_Loc[2-idx]
            if telescope == 'GMOSN':
                Spec = Spec_Loc[idx]
            if Binning == '2':
                Start = (1415,Spec)
            if Binning == '4':
                Start = (800,Spec)
            print(Start)
            data = hdulist[0].data
            data[data==0] = np.nan
            Trace, bkg, MASK1 = SP.Fit_Trace(data,Start,Range = 15, SClip = True,**DetecFlags)
            
            
            # write a fits file with the extracted trace
            TR=[]
            for idx,elem2 in enumerate(Trace):
                if elem2 > 15:
                    TR.append(data[int(elem2)-15:int(elem2)+15,idx])
                else: 
                    TR.append(data[0:30,idx])
            TR = np.array(TR)
            hdulist[0].data = np.transpose(TR)
            hdulist.writeto(elem.replace('.fits','') + 'Trace.fits' )
            Out_Files.append(elem.replace('.fits','') + 'Trace.fits')

            
#            SSpec = []
#            for FW in range(20):
#                Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=FW,Mask = MASK1,**DetecFlags)
#                SSpec.append(Spec1)
                
#            max_index, max_value = max(enumerate(np.nanmedian(SSpec,axis=1)/np.nanstd(SSpec,axis=1)), key=operator.itemgetter(1))
#            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=max_index+1,Mask = MASK1,**DetecFlags)
            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=Spec_FWHM,Mask = MASK1,**DetecFlags)
    
            
            
            if Binning == '2':
                Spec1N = Spec1/np.abs(np.nanmedian(Spec1[1500:1600]))
            if Binning == '4':
                Spec1N = Spec1/np.abs(np.nanmedian(Spec1[800:850]))

            f = open(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt','w')
            Out_Spec.append(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt')

            for item in Spec1N:
                f.write("%s\n" % item)
            
            f.close()

    if Diagnostic: 
        diag.create_website('Extract_Log.html')
        diag.add_Extract(Out_Files,Out_Spec,'Extract_Log.html')
            
        




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract spectrum')

    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    

    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    parser.add_argument('-g',
                        help='Generic name of the offset and spectra location')
    
    parser.add_argument('-fwhm',
                        default = 5,
                        help='FWHM of the spectrum trace')

    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  

    args = parser.parse_args()

    Verbose = args.v
    Spec_loc = args.g
    filenames = args.images  
    fwhm = int(args.fwhm)
    Diagnostic = args.d
    
    print(filenames)

    
    Extract_Spectrum(filenames,Verbose,Spec_loc,Diagnostic,fwhm)
    pass

