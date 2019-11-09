#!/usr/bin/env python

""" SP_BckgSub - Substract background to science data  
    v1.0: 2018-04-19, mdevogele@lowell.edu
  
"""
import argparse, shlex

from astropy.io import fits
import SP_Toolbox as SP

import simplejson

import SP_diagnostics as diag

from SP_CheckInstrument import CheckInstrument

import numpy as np


from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askopenfilenames
from tkinter.messagebox import showerror

from PIL import Image
from PIL import ImageTk
from PIL import ImageDraw



from past.utils import old_div
from scipy.ndimage import interpolation as interp


def BckgSub(FileName, Verbose, Method ,Suffix,Spec_loc,Diagnostic,Area = [250,350],test = 0, Live = False):
   
    telescope, obsparam = CheckInstrument([FileName[0]])
    Out_File = []
    print('telescope')
    for elem in FileName:
        
        hdulist = fits.open(elem)
        image = hdulist[0].data
        Mask = np.zeros_like(image)
        Mask[:,:] =  0
        
        
        if telescope == 'Deveny' or telescope == 'DEVENY' or telescope == 'SOAR':
            DetecFlags = {'Instrument':'Deveny'}
            if Method == 'auto':               
                Center = SP.Detect_Spectra(image,Bin =2,**DetecFlags)
                (Center)
                Area[0] = int(Center-70)
                Area[1] = int(Center+70)
            if Method == 'range':
                Center = (Area[1]-Area[0])/2                
            
            print(Area)
            
            Mask[0:Area[0],:] = 1
            Mask[Area[1]:,:] = 1
            Mask = Mask.astype(bool)
            
            image2 = np.ma.masked_array(image,mask = Mask)
            
            X = range(image.shape[0])
            X = np.array(X)
            
            Fdc = []
            
            if Live:
                median = np.median(image[int(image.shape[1]*0.25):
                                             int(image.shape[1]*0.75),
                                             int(image.shape[0]*0.25):
                                             int(image.shape[0]*0.75)])
                std    = np.std(image[int(image.shape[1]*0.25):
                                          int(image.shape[1]*0.75),
                                          int(image.shape[0]*0.25):
                                          int(image.shape[0]*0.75)])
    
                imgdat= image
                imgdat = old_div(np.clip(imgdat, median-0.5*std,
                                   median+0.5*std),(old_div(std,256)))
                imgdat = imgdat - np.min(imgdat)
    
                imgdat = interp.zoom(imgdat, test.zoom)            
            
            for i in range(image.shape[1]):
                index = np.argwhere(np.isnan(image2[:,i]))
                image2[index,i]=0
                image3 = SP.Sigma_Clip(image2[:,i])
                XX = np.ma.masked_array(X,image3.mask)
                z = np.ma.polyfit(XX,image3 , 1)
                p = np.poly1d(z)
                Fdc.append(p(Center))
                image[:,i] = image[:,i] - p(X)
                
                if Live: 
                    imgdat[:,int(i/2)] = old_div(np.clip(interp.zoom(image[:,i], test.zoom),0,100),(old_div(100,256)))
                    
                    test.images[test.index] = Image.fromarray(imgdat)
                    test.canvas.update()
                    test.canvas.delete("all")
                    test.canvas.forget()
            
            
                    test.tkimage = ImageTk.PhotoImage(test.images[0], palette=256)
                    test.canvas = Canvas(test.frame_fits, height=test.tkimage.height(), width=
                                 test.tkimage.width())
                    test.canvas.pack()
                    test.image = test.canvas.create_image(0, 0, anchor='nw',
                                                  image=test.tkimage)  
            
                  # select first image
            
                    test.index = 0
                    im = test.images[test.index]
                    test.tkimage.paste(im)                
                    test.canvas.update()
                
                   
            Out_File.append(elem.replace('.fits','').replace('_Procc','') + '_' + Suffix + '.fits')
            hdulist.writeto(elem.replace('.fits','').replace('_Procc','') + '_' + Suffix + '.fits',overwrite = True )
            np.savetxt(elem.replace('.fits','').replace('_Procc','') + '_' + Suffix + '.txt',Fdc)
        if telescope == 'GMOSS' or telescope == 'GMOSN':
            DetecFlags = {'Instrument':'GMOS'}
            OffFile = Spec_loc + '_Offset.txt'
            with open(OffFile,'r') as f:
                Offset = simplejson.load(f)
            SpecFile = Spec_loc + '_Spec_loc.txt'    
            with open(SpecFile,'r') as f:
                Spec_Loc = simplejson.load(f)
                
            Offset = np.sort(np.array(Offset).astype(int))
            Spec_Loc = np.sort(Spec_Loc)
            
            Off = int(hdulist[0].header['YOFFSET'])
            print(Off)
            idx = np.where(Offset == Off)[0][0]
    
            if telescope == 'GMOSS':
                Center = Spec_Loc[2-idx]
            if telescope == 'GMOSN':
                Center = Spec_Loc[idx]
                print(Center)
                
            Area[0] = int(Center-200)
            Area[1] = int(Center+200) 

            print(Area)
            
            Mask[0:Area[0],:] = 1
            Mask[Area[1]:,:] = 1
            Mask = Mask.astype(bool)
            
            image2 = np.ma.masked_array(image,mask = Mask)
            
            X = range(image.shape[0])
            X = np.array(X)
            
            
            for i in range(image.shape[1]):
                index = np.argwhere(np.isnan(image2[:,i]))
                image2[index,i]=0
                image3 = SP.Sigma_Clip(image2[:,i],sig = 3)
                XX = np.ma.masked_array(X,image3.mask)
                z = np.ma.polyfit(XX,image3 , 1)
                p = np.poly1d(z)
                image[:,i] = image[:,i] - p(X)
            
            hdulist.writeto(elem.replace('.fits','').replace('_Procc','').replace('Sub','') + '_' + Suffix + '.fits', v)
            Out_File.append(elem.replace('.fits','').replace('_Procc','') + '_' + Suffix + '.fits')
            
    if Diagnostic: 
        diag.create_website('Background-Sub_Log.html')
        diag.add_CosmCorr(Out_File,'Background-Sub_Log.html')






if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master flats')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
    parser.add_argument('-s',
                        help='Suffix to add to processed files',
                        default='Bckg')
    parser.add_argument('-m', help='Method to use for the selection of region to be considered for the evaluation of the background',
                        choices=['auto','range'],
                        default = 'auto')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-r',
                        help='Range of pixels to use for background subtraction. 2 arguments, the pixel center and the number of pixels to consider',
                        default = '300 200',
                        nargs=2)  
    
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  

    
    
    parser.add_argument('-g',
                        help='Generic name of the offset and spectra location')
    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
#    prefix = args.prefix
    Suffix = args.s
    Method = args.m
    Verbose = args.v
    Range = args.r
    filenames = args.images    
    Spec_loc = args.g
    Diagnostic = args.d

    
    Ran = []
    Ran.append(int(Range[0])-int(Range[1])/2)
    Ran.append(int(Range[0])+int(Range[1])/2)
    
    
    # call run_the_pipeline only on filenames
    BckgSub(filenames,Verbose,Method,Suffix,Spec_loc,Diagnostic,Area = Ran, test = 0, Live = False)
    pass