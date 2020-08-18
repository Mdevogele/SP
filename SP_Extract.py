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


def Extract_Spectrum(filename,Verbose,Spec_loc,Diagnostic,Spec_FWHM,Yloc,Live = 0,Live2 = False, Bckg = 'Auto'): 
    
    
    Out_Files = []
    Out_Spec = []
    telescope, obsparam = CheckInstrument([filename[0]])    
    DetecFlags = {'Instrument':telescope}
    if telescope == 'DEVENY' or telescope == 'SOAR' or telescope == 'SOAR 4.1m' or telescope =='NOT':
        for elem in filename:
            # Read the fits file 
            hdulist = fits.open(elem)
            # Extract data from the fits file
            data = hdulist[0].data
            # Output filename for the extracted background for errorbar computation
            if Bckg == 'Auto':
                Backg = np.loadtxt(elem.replace('.fits','.txt'))
                print(Backg)
            
            # Check the telescope which is used
            if telescope == 'DEVENY':
                DetecFlags = {'Instrument':'Deveny'}  
            if telescope == 'SOAR' or telescope == 'SOAR 4.1m':
                DetecFlags = {'Instrument':'Soar'}  
            if telescope == 'NOT':
                DetecFlags = {'Instrument':'ALFOSC_FASU'} 
            # Try to detect the spectrum to extract 
            
            if Yloc == False:
                Center = SP.Detect_Spectra(data,Bin=2,**DetecFlags)
            else:
                Center = Yloc


            if telescope == 'DEVENY':
                Start = (1415,Center)
            if telescope == 'SOAR' or telescope == 'SOAR 4.1m':
                Start = (1415,Center)
            if telescope == 'NOT':
                Start = (1000,Center)                
#            Start = (1415,Center)
            
            Trace, bkg, MASK1 = SP.Fit_Trace(hdulist,Start,Range = 15, SClip = True,Live = Live,Live2 = Live2, **DetecFlags)

            TR=[]
            for idx,elem2 in enumerate(Trace):
                if elem2 > 50 and elem2 < np.size(data,0)-50:
                    TR.append(data[int(elem2)-45:int(elem2)+45,idx])
                else: 
                    TR.append(data[0:90,idx])
            TR = np.array(TR)
            hdulist[0].data = np.transpose(TR)
            hdulist.writeto(elem.replace('.fits','') + 'Trace.fits',overwrite = True,output_verify="ignore")
            Out_Files.append(elem.replace('.fits','') + 'Trace.fits')

#            SSpec = []
#            for FW in range(20):
#                Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=FW,Mask = MASK1,**DetecFlags)
#                SSpec.append(Spec1)
                
#            max_index, max_value = max(enumerate(np.nanmedian(SSpec,axis=1)/np.nanstd(SSpec,axis=1)), key=operator.itemgetter(1))
            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=Spec_FWHM,Mask = MASK1,**DetecFlags)
            
#            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=50,Mask = MASK1,**DetecFlags)
            
            # Error bar computation
            print(Bckg)
            Err = np.sqrt(np.array(Spec1).astype(float)*2.1) + np.sqrt(Spec_FWHM*Backg*2.1) + np.sqrt(4)
            Err = Err/2.1
            # Normalization of the spectrum 
            print(Spec1 )
            if telescope == 'NOT':
                Median = np.abs(np.nanmedian(Spec1[200:250]))                
            else:
                Median = np.abs(np.nanmedian(Spec1[1500:1600]))
            Spec1N = Spec1/Median
            Err= Err/Median
            fname = elem.replace('_CosmCorr','').replace('Bckg','Extracted').replace('.fits','.txt')
            np.savetxt(fname, np.array([Spec1N,Err]).transpose())
#            f = open(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt','w')
#            Out_Spec.append(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt')
#            for item in Spec1N:
#                f.write("%s\n" % item)
#            f.close()
    if telescope == 'GMOSS' or telescope == 'GMOSN':
        
        # Chip gaps for
        #
        # Hamamatsu:
        #
        # GMOS S: 61 pixels
        # GMOS N: 67 pixels
        #
        # Other chips:
        #
        # 37 pixels
        
        OffFile = Spec_loc + '_Offset.txt'
        with open(OffFile,'r') as f:
            Offset = simplejson.load(f)
        SpecFile = Spec_loc + '_Spec_loc.txt'    
        with open(SpecFile,'r') as f:
            Spec_Loc = simplejson.load(f)
        
        Offset = np.sort(np.array(Offset).astype(int))
        Spec_Loc = np.sort(Spec_Loc)
        for elem in filename:
            
            # Read the fits file 
            hdulist = fits.open(elem)
            
            # Extract data from the fits file
            data = hdulist[0].data
            data[data==0] = np.nan # Convert 0 to np.nan
            
            # Output filename for the extracted background for errorbar computation
            if Bckg == 'Auto':
                Backg = np.loadtxt(elem.replace('.fits','.txt'))
            else:
                Backg = np.loadtxt(Bckg)
                
            
            # Get grating information    # USED ? 
            Gratting = hdulist[0].header[obsparam['grating']]
            
            # Get binning information    # USED ? 
            Binning = hdulist[1].header['CCDSUM'][0]
            
            # Get detector information   # USED ?  
            Detector = hdulist[0].header['DETECTOR']
            
            # Get the detectors gain and read-out-noise
            Gain = []
            RON = []
            for HDU in hdulist: # loop over all header in the fits file
                try:
                    Gain.append(HDU.header['GAIN'])
                    RON.append(HDU.header['RDNOISE'])
                except:
                    pass
            # Use the mean as an approximation of the gain for the whole detector (only used for errorbar computation)
            Gain = np.mean(Gain)
            RON = np.mean(RON)
                
            
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
            Trace, bkg, MASK1 = SP.Fit_Trace(hdulist,Start,Range = 15, SClip = True,**DetecFlags)
                        
#            # write a fits file with the extracted trace
#            TR=[]
#            for idx,elem2 in enumerate(Trace):
#                if elem2 > 15:
#                    TR.append(data[int(elem2)-15:int(elem2)+15,idx])
#                else: 
#                    TR.append(data[0:30,idx])
#            TR = np.array(TR)
#            hdulist[0].data = np.transpose(TR)
#            hdulist.writeto(elem.replace('.fits','') + 'Trace.fits' )
#            Out_Files.append(elem.replace('.fits','') + 'Trace.fits')

            
#            SSpec = []
#            for FW in range(20):
#                Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=FW,Mask = MASK1,**DetecFlags)
#                SSpec.append(Spec1)
                
#            max_index, max_value = max(enumerate(np.nanmedian(SSpec,axis=1)/np.nanstd(SSpec,axis=1)), key=operator.itemgetter(1))
#            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=max_index+1,Mask = MASK1,**DetecFlags)
            Spec1, Backgroung = SP.Extract_Spectrum(data,Trace,Backg,FWHM=Spec_FWHM,Mask = MASK1,**DetecFlags)
    

            # Error bar computation
            
            Err = np.sqrt(np.array(Spec1).astype(float)*Gain) + np.sqrt(Spec_FWHM*Gain*np.array(Backgroung).astype(float)) + np.sqrt(RON)
            Err = Err/Gain            
            
            
            if Binning == '2':
                Spec1N = Spec1/np.abs(np.nanmedian(Spec1[1500:1600]))
                Err= Err/np.abs(np.nanmedian(Spec1[1500:1600]))
            if Binning == '4':
                Spec1N = Spec1/np.abs(np.nanmedian(Spec1[800:850]))
                Err= Err/np.abs(np.nanmedian(Spec1[800:850]))
#            f = open(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt','w')
#            Out_Spec.append(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt')
            
            
            
            
            # Error bar computation
            
#            Err = np.sqrt(np.array(Spec1).astype(float)*Gain) + np.sqrt(Spec_FWHM*Gain*np.array(Backgroung).astype(float)) + np.sqrt(RON)
#            Err = Err/Gain
            # Normalization of the spectrum 
#            Median = np.abs(np.nanmedian(Spec1[1500:1600]))
#            Spec1N = Spec1/Median
#            Err= Err/Median
            
            
            
            fname = elem.replace('_CosmCorr','').replace('Bckg','Extracted').replace('.fits','.txt')
            np.savetxt(fname, np.array([Spec1N,Err]).transpose())
            
            
            
            

#            for item in Spec1N:
#                f.write("%s\n" % item)
            
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
    
    parser.add_argument('-Yloc',
                        default = False,
                        help='Y location of the spectrum')
    

    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  

    args = parser.parse_args()

    Verbose = args.v
    Spec_loc = args.g
    Yloc = args.Yloc
    filenames = args.images  
    fwhm = int(args.fwhm)
    Diagnostic = args.d
    
    print(filenames)

    
    Extract_Spectrum(filenames,Verbose,Spec_loc,Diagnostic,fwhm,Yloc,Live = 0,Live2 = False, Bckg = 'Auto')
    pass

