#!/usr/bin/env python

""" SP_Extract - Extract spectrum
    v1.0: 2018-04-19, mdevogele@lowell.edu
  
"""

import argparse, shlex

import SP_Toolbox as SP
import SP_Deveny_Toolbox as SP_Dev_Tool
import numpy as np
import SP_diagnostics as diag

import simplejson

from astropy.io import fits
import operator

from SP_CheckInstrument import CheckInstrument


def Extract_Spectrum(filename,Verbose,Spec_loc): 
    
    telescope, obsparam = CheckInstrument([filename[0]])    
    DetecFlags = {'Instrument':telescope}
    if telescope == 'DEVENY':
        for elem in filename:
            hdulist = fits.open(elem)
            data = hdulist[0].data
            DetecFlags = {'Instrument':'Deveny'}  
            Center = SP.Detect_Spectra(data,Bin=2,**DetecFlags)
            Start = (1415,Center)
            
            Trace, bkg, MASK1 = SP.Fit_Trace(data,Start,Range = 15, SClip = True,**DetecFlags)

            TR=[]
            for idx,elem2 in enumerate(Trace):
                if elem2 > 15:
                    TR.append(data[int(elem2)-15:int(elem2)+15,idx])
                else: 
                    TR.append(data[0:30,idx])
            TR = np.array(TR)
            hdulist[0].data = np.transpose(TR)
            hdulist.writeto(elem.replace('.fits','') + 'Trace.fits' )

            SSpec = []
            for FW in range(20):
                Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=FW,Mask = MASK1,**DetecFlags)
                SSpec.append(Spec1)
                
            max_index, max_value = max(enumerate(np.nanmedian(SSpec,axis=1)/np.nanstd(SSpec,axis=1)), key=operator.itemgetter(1))
            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=max_index+1,Mask = MASK1,**DetecFlags)
            
#            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=50,Mask = MASK1,**DetecFlags)
            Spec1N = Spec1/np.abs(np.nanmedian(Spec1[1500:1600]))
            f = open(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt','w')
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
            
            SSpec = []
            for FW in range(20):
                Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=FW,Mask = MASK1,**DetecFlags)
                SSpec.append(Spec1)
                
            max_index, max_value = max(enumerate(np.nanmedian(SSpec,axis=1)/np.nanstd(SSpec,axis=1)), key=operator.itemgetter(1))
            Spec1 = SP.Extract_Spectrum(data,Trace,bkg,FWHM=max_index+1,Mask = MASK1,**DetecFlags)
    
            
            
            if Binning == '2':
                Spec1N = Spec1/np.abs(np.nanmedian(Spec1[1500:1600]))
            if Binning == '4':
                Spec1N = Spec1/np.abs(np.nanmedian(Spec1[800:850]))

            f = open(elem.replace('.fits','').replace('_bkgSub','').replace('_Procc','') + '.txt','w')
            for item in Spec1N:
                f.write("%s\n" % item)
            
            f.close()

        
            
        




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract spectrum')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
#    parser.add_argument('-target', help='primary targetname override',
#                        default=None)
#    parser.add_argument('-m', help='add flats information to diagnostic.htlm file',
#                        choices=['auto','range'])
#    parser.add_argument('-s', help='If there is several series of flat \n || Options: none: Only one serie \n || index: split according to the index of the files \n || target: split according to the target name in fits headers \n || pointing: split according to telescope pointing according to fits headers',
#                        choices=['none','index','target','pointing'],
#                        default = 'None')
#    parser.add_argument('-b',
#                        help='Name of the master bias to use \n || Can use None if no bias are available',
#                        default='MasterBias.fits')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
#    parser.add_argument('-r',
#                        help='Range of pixels to use for background subtraction',
#                        nargs=2)    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    parser.add_argument('-g',
                        help='Generic name of the offset and spectra location')

    args = parser.parse_args()
#    prefix = args.prefix
#    man_targetname = args.target
#    Method = args.m
#    Series = args.s
#    MasterBias = args.b
    Verbose = args.v
    Spec_loc = args.g
#    Range = args.r
    filenames = args.images  
    
    print(filenames)

    
    
    # call run_the_pipeline only on filenames
    Extract_Spectrum(filenames,Verbose,Spec_loc)
    pass

