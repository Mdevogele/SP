#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 14:37:19 2017

@author: maximedevogele

Spectroscopy pipeline: 
    

"""

#! /usr/bin/env python

#import sys

import fileSelect as fs
import os, sys
import copy

import numpy as np
import operator
from astropy.io import fits
import sqlite3
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
try:
    import imreg_dft as ird
except ImportError:
    print('Module imreg_dft not found. Please install with: pip install imreg_dft')
    sys.exit()
from scipy.optimize import curve_fit
import itertools
import datetime
import numpy.ma as ma
import prompt_hack
import csv
from matplotlib.cbook import get_sample_data

from scipy.interpolate import UnivariateSpline
from astropy.convolution import convolve, Box1DKernel


from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astroquery.simbad import Simbad
import astropy.units as u
import astropy.coordinates as coord
from SP_CheckInstrument import CheckInstrument


import _SP_conf



CHIP_GAP_B4 = ([500,528],[1029,1058])
CHIP_GAP_B2 = ([1007,1063],[2062,2114])

CHIP_GAP2_B4 = ([500,528],[1029,1058])
CHIP_GAP2_B2 = ([732,752],[1426,1446])

Pipe_Path = _SP_conf.rootpath

def Check(filename):
    
    instruments =[] 
    for elem in filename:
        if not os.path.isfile(elem):
            raise IOError('No such file: ' + str(elem))
            
        try:
            hdulist = fits.open(elem, ignore_missing_end=True)
        except IOError:
            print('ERROR: cannot open file %s' % elem)
            continue

        header = hdulist[0].header
        for key in _SP_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])
                break

    if len(np.unique(instruments)) == 0:
        raise KeyError('cannot identify telescope/instrument; please update' + \
                       '_SP_conf.instrument_keys accordingly')

    if len(np.unique(instruments)) > 1:
        raise Warning('More than one instrument identified')    
    
    
    if len(filename) == 0:
        raise IOError('cannot find any data...')

def Auto_Detect_Lines(Arcs, Tresh_Det = 1.5, Tresh_Arcs = [8, 20]):


    Arcs[np.isnan(Arcs)] = 0
    
    Arcs_loc = []
    max_index, max_value = max(enumerate(Arcs), key=operator.itemgetter(1))
    plt.plot(Arcs)
    while np.median(Arcs) < max_value-5*np.std(Arcs):
        max_index, max_value = max(enumerate(Arcs), key=operator.itemgetter(1))
        # Search for the size of the arc lines 
        
        value = 999999
        min_ind = max_index
        while value> max_value/Tresh_Det:
            min_ind -= 1
            value = Arcs[min_ind]
    
        value = 999999
        max_ind = max_index
        while value> max_value/Tresh_Det:
            max_ind += 1
            value = Arcs[max_ind]        
            
        Arcs_Size = (max_ind - min_ind)
        if Arcs_Size < Tresh_Arcs[1] and Arcs_Size > Tresh_Arcs[0] :
            print(Arcs_Size)
            Arcs_loc.append(int(max_ind+min_ind)/2)
        
        if min_ind>4:
            Arcs[min_ind-4:max_ind+4] = np.nanmedian(Arcs)
        else:
            Arcs[:max_ind+4] = np.nanmedian(Arcs)
    
    for elem in Arcs_loc:
        plt.plot([elem,elem],[0,60000])
    
    return Arcs_loc

    

def Get_Consecutive(seq):
    
    """ From a list on numbers, provide the sublists of consecutive numbers 
     Ex: MyList = [0,1,2,5,6,7,9,10,11]
         OutLists = [[0,1,2],[5,6,7],[9,10,11]] """
    
    Series = []
    subseries=[]
    for idx,elem in enumerate(seq):
        if idx == 0:
            subseries.append((0,elem))
        elif seq[idx-1] +1  == seq[idx]:
            subseries.append((idx ,elem))
        else:
            Series.append((subseries))
            subseries = []
            subseries.append((idx,elem))
    Series.append((subseries))
    
    return Series




## Which one is actually needed? 
def Get_Binning(files):
    
    """ derive binning from image header"""
    
    hdulist = fits.open(files)
    
    Binning = hdulist[1].header['CCDSUM']
    Bin = int(Binning[0])
    
    return Bin

def get_binning(header, obsparam):
    
    """ derive binning from image header
        use obsparam['binning'] keywords, unless both keywords are set to 1
        return: tuple (binning_x, binning_y)"""

    if obsparam['binning'][0] == 1 and obsparam['binning'][1] == 1:
        binning_x = 1
        binning_y = 1
    elif '#' in obsparam['binning'][0]:
        if '#blank' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].\
                                     split('#')[0]].split()[0])
            binning_y = float(header[obsparam['binning'][1].\
                                     split('#')[0]].split()[1])
        elif '#x' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].\
                                     split('#')[0]].split('x')[0])
            binning_y = float(header[obsparam['binning'][1].\
                                     split('#')[0]].split('x')[1])
        elif '#CH#' in obsparam['binning'][0]:
            # only for RATIR
            channel = header['INSTRUME'].strip()[1]
            binning_x = float(header[obsparam['binning'][0].
                                     replace('#CH#', channel)])
            binning_y = float(header[obsparam['binning'][1].
                                     replace('#CH#', channel)])
    else:
        binning_x = header[obsparam['binning'][0]]
        binning_y = header[obsparam['binning'][1]]

    return (binning_x, binning_y)




def Plot_Taxonomy(Wav,Spec,SpecName,Date,Facility):
    
    """ Generate a figure containing the asteroid spectrum with overlay the best taxonomic type """
    
    # Variables definition
    Wv = []
    Data_Tax = []
    Data_Err = []
    
    # Read the busdemeo-meanspectra.csv containing the templates for each taxonpmic type
    
    f = open(Pipe_Path +'/busdemeo-meanspectra.csv', 'rU')
    Tax = csv.reader(f)
    T1 = next(Tax)
    T1 = next(Tax)
    Head = T1
    Head_Tax = np.array(Head[1::2])
    Head_Err = np.array(Head[2::2])

    for row in Tax:
        Wv.append(row[0])
        Data_Tax.append(row[1::2])
        Data_Err.append(row[2::2])
    
    Wv = np.array(Wv).astype(float)
    
    Data_Tax = np.array(Data_Tax).astype(float)
    Data_Err = np.array(Data_Err).astype(float)

    Data_Tax = Data_Tax.T

    Data_Err = Data_Err.T

    # Some types have undefined enveloppes
    Data_Err[4,:] = 0.03
    Data_Err[10,:] = 0.03
    Data_Err[12,:] = 0.03

    # Rebin the asteroid spectra to match the taxonomic templates
    Wave = np.linspace(0.425,0.975,12)
    Dw = (0.975-0.425)/11
    
    Wave_out = []
    Spec_out = []
    Error = []
    
    for k in Wave:
        Wave_out.append(k+Dw/2)
        Spec_Tamp = []
    
        condition1 = Wav > k
        condition2 = Wav < k + Dw
        condition = condition1*condition2
        Spec_Tamp.append(np.nanmedian(np.extract(condition, Spec)))
        Spec_out.append(Spec_Tamp)
        Error.append(np.nanstd(np.extract(condition, Spec)*~np.isinf(np.extract(condition, Spec)))/np.sqrt(len(np.extract(condition, Spec))-1))

    Wave_out = np.array(Wave_out).astype(float)
    Spec_out = np.array(Spec_out).astype(float).flatten()
    Spec_out = Spec_out/Spec_out[2]
    Error = np.array(Error).astype(float).flatten()
    
    Error[Error == 0] = 100

    
    ChiSq = [] 
    for Taxon,Err in zip(Data_Tax,Data_Err):
        Taxon = np.array(Taxon).astype(float)
        Err = np.array(Err).astype(float)
        Err[2]=1
        ChiSq.append(np.sqrt(np.nansum((Spec_out-Taxon[0:12])**2/Error**2)))

    Index = []
    Chi = []
    min_index, min_value = min(enumerate(ChiSq), key=operator.itemgetter(1))

    # Rebin the asteroid spectra for plotting 
    Wave = np.linspace(0.39,0.99,31)
    Dw = (0.98-0.38)/30
    
    Wave_out = []
    Spec_out = []    
    Error = []
    
    for k in Wave:
        Wave_out.append(k+Dw/2)
        Spec_Tamp = []
    
        condition1 = Wav > k
        condition2 = Wav < k + Dw
        condition = condition1*condition2
        Spec_Tamp.append(np.nanmedian(np.extract(condition, Spec)))
        Spec_out.append(Spec_Tamp)
        Error.append(np.nanstd(np.extract(condition, Spec)*~np.isinf(np.extract(condition, Spec)))/np.sqrt(len(np.extract(condition, Spec))-1))

    Wave_out = np.array(Wave_out).astype(float)
    Spec_out = np.array(Spec_out).astype(float).flatten()
    Spec_out = Spec_out#/Spec_out[5] 
    Error = np.array(Error).astype(float).flatten()
    Taxonomy = Head_Tax[min_index].split('_')

    Num_Tax = 12
    
    
    # Generate the plot
    
    plt.plot(Wav,Spec,zorder=1)
    plt.errorbar(Wave_out,Spec_out,yerr = Error,fmt='.k',markersize=12,zorder=3, label= SpecName + ' (Binned)')
    plt.fill_between(Wv[0:Num_Tax], Data_Tax[min_index,0:Num_Tax] - Data_Err[min_index,0:Num_Tax], Data_Tax[min_index,0:Num_Tax] + Data_Err[min_index,0:Num_Tax],alpha=0.8,facecolor='grey',zorder=2, label= 'Taxonomy: ' + str(Taxonomy[0]))    

    plt.legend()

    axes = plt.gca()
    axes.set_ylim(0.55,1.6)

    fig = plt.gcf()
    
#    plt.title(Facility + ', ' + SpecName )
    plt.title("Lowell's 4.3m DCT; DeVeny; 2019-05-28")
    plt.xlabel(r'Wavelength [micron]', fontsize=14)
    plt.ylabel(r'Normalized reflectance', fontsize=14)
    plt.text(0.35,1.52, "Observers: B. Skiff, N. Moskovitz")
    plt.text(0.35,1.47, "Reduction: M. Devogele")


    im = plt.imread(get_sample_data(Pipe_Path +'/manos_splash.eps'))
    newax = fig.add_axes([0.82, 0.82, 0.12, 0.12], anchor='NE', zorder=1)
    newax.imshow(im)
    newax.axis('off')


    axes.set_ylim(0.5,1.6)
    
    f = open(SpecName + '_' + Date + '_' + Facility + '_spec.dat','w')
    
    for Wavel, Int, Error in zip(Wave_out,Spec_out,Error):
        f.write(str(Wavel) + ',' + format(Int, '.4f') + ',' + format(Error, '.4f') + '\n')
        
    f.close()

def Get_Taxonomy(Wav,Spec):
    Wv = []
    f = open(Pipe_Path +'/busdemeo-meanspectra.csv', 'rU')
    Tax = csv.reader(f)
    T1 = next(Tax)
    T1 = next(Tax)
    Head = T1
    Head_Tax = np.array(Head[1::2])
    Head_Err = np.array(Head[2::2])
    Data_Tax = []
    Data_Err = []
    for row in Tax:
        Wv.append(row[0])
        Data_Tax.append(row[1::2])
        Data_Err.append(row[2::2])
    
    Data_Tax = np.array(Data_Tax).astype(float)
    Data_Err = np.array(Data_Err).astype(float)

    Data_Tax = Data_Tax.T
    Data_Err = Data_Err.T

    Data_Err[4,:] = 0.03
    Data_Err[10,:] = 0.03
    Data_Err[12,:] = 0.03

    Wave = np.linspace(0.425,0.975,12)
    Dw = (0.975-0.425)/11
    
    Wave_out = []
    Spec_out = []
    Error = []
    
    for k in Wave:
        Wave_out.append(k+Dw/2)
        Spec_Tamp = []
    
        condition1 = Wav > k
        condition2 = Wav < k + Dw
        condition = condition1*condition2
        Spec_Tamp.append(np.nanmedian(np.extract(condition, Spec)))
        Spec_out.append(Spec_Tamp)
        Error.append(np.nanstd(np.extract(condition, Spec)*~np.isinf(np.extract(condition, Spec)))/np.sqrt(len(np.extract(condition, Spec))-1))

    Wave_out = np.array(Wave_out).astype(float)
    Spec_out = np.array(Spec_out).astype(float).flatten()
    Spec_out = Spec_out/Spec_out[2]
    print(Spec_out)
    Error = np.array(Error).astype(float).flatten()

    Error[Error == 0] = 100
    
    ChiSq = [] 
    for Taxon,Err in zip(Data_Tax,Data_Err):
        Taxon = np.array(Taxon).astype(float)
        Err = np.array(Err).astype(float)
        Err[2]=1
        ChiSq.append(np.sqrt(np.nansum((Spec_out-Taxon[0:12])**2/Error**2)))

    Index = []
    Chi = []
    min_index, min_value = min(enumerate(ChiSq), key=operator.itemgetter(1))
    ChiSq[min_index] = 999
    Index.append(min_index)
    Chi.append(min_value)
#    plt.plot(Wv[0:12],Data_Tax[min_index,0:12],'.')
    min_index, min_value = min(enumerate(ChiSq), key=operator.itemgetter(1))
    ChiSq[min_index] = 999
    Index.append(min_index)
    Chi.append(min_value)
#    plt.plot(Wv[0:12],Data_Tax[min_index,0:12],'.')
    min_index, min_value = min(enumerate(ChiSq), key=operator.itemgetter(1))
    ChiSq[min_index] = 999
    Index.append(min_index)
    Chi.append(min_value)
    min_index, min_value = min(enumerate(ChiSq), key=operator.itemgetter(1))
    ChiSq[min_index] = 999
    Index.append(min_index)
    Chi.append(min_value)
    min_index, min_value = min(enumerate(ChiSq), key=operator.itemgetter(1))
    ChiSq[min_index] = 999
    Index.append(min_index)
    Chi.append(min_value)
#    plt.plot(Wv[0:12],Data_Tax[min_index,0:12],'.')
#    plt.plot(Wv[0:12],Spec_out)
    
    for head,chi in zip(Head_Tax[Index],Chi):
        print(str(head) + '\t' +  str(chi))
        
    return zip(Head_Tax[Index],Chi)
    

def Cut_Image(data,Bin = 4):
    
    if Bin == 4:
        In_y = 400
        Fin_y = 700
        In_x = 670
        Fin_x = 700
    if Bin == 2:
        In_y = 800
        Fin_y = 1400
        In_x = 1530
        Fin_x = 1600
        
    xs = [In_x,Fin_x]
    ys = [In_y,Fin_y]
        
    SS = data[In_y:Fin_y,In_x:Fin_x]
        
    Trace = np.median(SS,axis=1)
    
    return Trace, xs, ys

def Auto_Detect_Spectra(files, Sig, Auto = True):

    hdulist = fits.open(files)
    data = hdulist[0].data
    
    data = np.nan_to_num(data)
    
    Binning = hdulist[1].header['CCDSUM']
    Bin = int(Binning[0])
    
    image,XS,YS = Cut_Image(data, Bin = Bin)
    
    Range = 30    
    image -= np.nanmedian(image)
#    image = np.abs(image)

    image[0:30] = 0
    image[-30:] = 0
    
    std = np.nanstd(image)
    med = np.nanmedian(image)
    max_index, max_value = max(enumerate(image), key=operator.itemgetter(1))
    
    print(max_value)
    print(std)
    print(med)
    
    Spec = []
    print(Sig)
    
    if Auto == True:
        counter = 0
        while max_value > med + Sig*std and counter < 10:
            counter +=1
            xs = range(max_index-Range,max_index+Range)
            ys = image[xs]
            p0 = [0,max_value,2.24,1.41,max_index] 
            Res = []
            for Conv in range(10):
                aa = np.zeros(31)
                aa[15-Conv:15+Conv+1] = 1./(Conv*2+1)
                try:
                    coeff, var_matrix = curve_fit(lambda x, S0,S1,a, b,x0: MOFFAT_NEW(x, S0,S1,a, b,x0 , Conv),xs,ys,p0,maxfev = 100000)
                    fit = MOFFAT(xs, *coeff)
                    Res.append(np.sum((image[xs] - np.convolve(fit,aa,'same'))**2))
                except RuntimeError:
                    Res.append(99999999999)
            
            min_index, min_value = min(enumerate(Res), key=operator.itemgetter(1))
            
            Conv = min_index
            
            aa = np.zeros(31)
            aa[15-Conv:15+Conv+1] = 1./(Conv*2+1)
            coeff, var_matrix = curve_fit(lambda x, S0,S1,a, b,x0: MOFFAT_NEW(x, S0,S1,a, b,x0 , Conv),xs,ys,p0,maxfev = 100000)
            fit = MOFFAT(xs, *coeff)
    
            image[xs] = image[xs] - np.convolve(fit,aa,'same')
        
            std = np.nanstd(image)
            med = np.nanmedian(image)
        
            Spec.append(max_index)
            max_index, max_value = max(enumerate(image), key=operator.itemgetter(1))
            print(max_value)
            print(std)
            print(med)
    else:
        print('Manual detection of the spectra')
        plt.figure()
        plt.plot(image)
        prompt_hack.start()
        Spec = prompt_hack.input('Enter the coordinates of the spectra: ')
        prompt_hack.finish()
        plt.close()
        Spec = np.array(Spec)
    
    Spec = np.array(Spec)
    Spec = np.unique(Spec)
    return Spec


def Get_Spectra(Spec,Offset):
    Perm1 = np.array(list(itertools.permutations(Spec)))
    Perm2 = np.unique(Perm1[:,0:3],axis = 0)
    if len(Spec) >= len(Offset):
        Correl = []
        for Ra in Perm2:
            Correl.append(np.correlate((Ra[0:len(Offset)] - np.mean(Ra[0:len(Offset)])) / (np.std(Ra[0:len(Offset)]) * len(Ra[0:len(Offset)])),(Offset - np.mean(Offset))/np.std(Offset)))
    max_index, max_value = max(enumerate(Correl), key=operator.itemgetter(1))
    Pos_Spec = Perm2[max_index][0:3]
    
    counter = 0
    while np.std((Pos_Spec-Offset)/np.mean(Pos_Spec-Offset)) >0.1 and counter < 100:
        counter += 1
        print(counter)
        Correl[max_index] = 0
        max_index, max_value = max(enumerate(Correl), key=operator.itemgetter(1))
        Pos_Spec = Perm2[max_index][0:3]

    if counter == 100:
        print('Spectra not detected \n')
        Pos_Spec = False
    
    return Pos_Spec
    
    

def Get_Offset(image):
    hdulist = fits.open(image)
    header = hdulist[0].header
    Offset = header['YOFFSET']
    
    return Offset

def Image_Subtract(im1,im2,out,data = 0):
    
    hdulist1 = fits.open(im1 )
    hdulist2 = fits.open(im2)
    
    hdulist1[data].data = hdulist1[data].data - hdulist2[data].data
    
    hdulist1.writeto(out)
    
def Image_Add(im1,im2,out,data = 0):
    
    hdulist1 = fits.open(im1 )
    hdulist2 = fits.open(im2)
    
    hdulist1[data].data = hdulist1[data].data + hdulist2[data].data
    
    hdulist1.writeto(out)    


## is GMOS_open not used anymore? 
def GMOS_open2(image):
    
    telescope, obsparam = CheckInstrument(image)
    hdulist = fits.open(image[0])
    
    Bin = get_binning(hdulist[1].header,obsparam)
    
    Data =[]
    DATASEC = []
    
    
    if hdulist[0].header['DETECTOR'] == 'GMOS + e2v DD CCD42-90':
        print('bla')
        Order = []
        for idx, elem in enumerate(hdulist):
            try: 
                Order.append(int(elem.header['FRAMEID']))
            except KeyError:
                pass

        Ord = np.argsort(Order) + 1             
            
        for idx, elem in enumerate(hdulist):
            if idx < 6:
                try: 
                    NCCD = hdulist[Ord[idx]].header['NCCDS']
                except KeyError:
                    pass
                try:
                    print(hdulist[Ord[idx]].header['FRAMEID'])
                    DATASEC = hdulist[Ord[idx]].header['DATASEC']
                    print(DATASEC)
                    data = hdulist[Ord[idx]].data
                    Xaxis = hdulist[Ord[idx]].header['NAXIS1']
                    Yaxis = hdulist[Ord[idx]].header['NAXIS2']
                    Shape = hdulist[Ord[idx]].data.shape
                    if Bin[0] == 4:
                        GAP = np.zeros([Shape[0]-1,11])
                    if Bin[0] == 2:   
                        GAP = np.zeros([Shape[0]-1,22])
                    if Bin[0] == 1:   
                        GAP = np.zeros([Shape[0]-1,44])    
                    GAP[:,:] = np.nan
                    if int(hdulist[Ord[idx]].header['FRAMEID']) == 2 or int(hdulist[Ord[idx]].header['FRAMEID']) == 4:
                        Data.append(GAP)
                    Data.append(data[int(DATASEC.replace('[','').replace(']','').split(',')[1].split(':')[0])-1 : int(DATASEC.replace('[','').replace(']','').split(',')[1].split(':')[1])-1, int(DATASEC.replace('[','').replace(']','').split(',')[0].split(':')[0]):int(DATASEC.replace('[','').replace(']','').split(',')[0].split(':')[1])+1])
                except KeyError:
                    pass
    
        hdulist[0].data = np.concatenate(Data[:],axis=1)
    else:    
     
        for idx, elem in enumerate(hdulist):
            try: 
                NCCD = elem.header['NCCDS']
            except KeyError:
                pass
            try:
                DATASEC = elem.header['DATASEC']
                Xaxis = elem.header['NAXIS1']
                Yaxis = elem.header['NAXIS2']
                Shape = elem.data.shape
                GAP = np.zeros([Shape[0]-1,int(76/Bin[0])])
                GAP[:,:] = np.nan
                if int(elem.header['FRAMEID']) == 4 or int(elem.header['FRAMEID']) == 8:
                    Data.append(GAP)
                Data.append(elem.data[int(DATASEC.replace('[','').replace(']','').split(',')[1].split(':')[0])-1 : int(DATASEC.replace('[','').replace(']','').split(',')[1].split(':')[1])-1, int(DATASEC.replace('[','').replace(']','').split(',')[0].split(':')[0])-1:int(DATASEC.replace('[','').replace(']','').split(',')[0].split(':')[1])-1])
            except KeyError:
                pass
    
        hdulist[0].data = np.concatenate(Data[:],axis=1)
    
    return hdulist


def GMOS_open(image):
    
    telescope, obsparam = CheckInstrument(image)
    hdulist = fits.open(image[0])
    
    
    X1 = hdulist[1].data
    X2 = hdulist[2].data
    X3 = hdulist[3].data
    X4 = hdulist[4].data
    
    if hdulist[0].header[obsparam['grating']] != 'MIRROR': 
        X5 = hdulist[5].data
        X6 = hdulist[6].data
        X7 = hdulist[7].data
        X8 = hdulist[8].data
        X9 = hdulist[9].data
        X10 = hdulist[10].data
        X11 = hdulist[11].data
        X12 = hdulist[12].data
    
        Chip_Gap = np.zeros((2088, 52))
    
        X = np.concatenate([X1[24:,1:256],X2[24:,32:288],X3[24:,0:256],X4[24:,32:273],Chip_Gap,X5[24:,6:256],X6[24:,32:288],X7[24:,0:256],X8[24:,32:272],Chip_Gap,X9[24:,6:256],X10[24:,32:288],X11[24:,0:256],X12[24:,32:288]],axis = 1)    
    else:
        X = np.concatenate([X1[24:,0:255],X2[24:,32:287],X3[24:,0:255],X4[24:,32:287]],axis = 1)    
    
    hdulist[0].data = X 
    
    return hdulist
    




def Create_Flat(image_list,**kw):
    ''' 
    Create Master_Flat
    
    !!!! the region to consider for the normalization depend on the disperser and the binning
    for now just work for R150 bin 2x2
    
    need to implement the interpolation in the chip gap for the flats
    
    '''


    ###### Parse the argument sent through the function #######

    # Do the program displays text during the processing ?
    # [False], True
    if 'Verbose' in kw:
        verbose = kw['Verbose']
    else:
        verbose = False
    
    # Do the program create a fits file with the MasterFlat ?
    # WriteFile = True : Create a fits file named MasterFlat.fits
    # WriteFile = False : Do not create any fits file
    # WriteFile = 'any string' : Create a fits file named the string given
    # [False]
    if 'WriteFile' in kw:
        WriteFile = kw['WriteFile']
        if WriteFile is True:
            print('No flat name was provided, MasterFlat.fits is used as default')
            NameFlat = 'MasterFlat.fits'
        NameFlat = kw['WriteFile']
        WriteFile = True
    else:
        WriteFile = False
    
    
    # Indicate where the raw data are located 
    # Take the current directory as default 
    if 'RawPath' in kw:
        RawPath = kw['RawPath']
    else:
        RawPath = './'  
        
        
    # Check if a bias file is provided
    # Stop the program if not
    if 'Bias' in kw:
        if isinstance(kw['Bias'], str):
            hdulist = fits.open(kw['Bias'])
            MasterBias = hdulist[0].data
        if isinstance(kw['Bias'], np.ndarray):
            MasterBias = kw['Bias']
    else:
        print('****************************')
        print('ERROR: No bias file provided')
        print('Execution stopped')
        return

    # Do you want to want to overwrite if the file already exist ? 
    # [False]
    if 'OverWrite' in kw:
        OverWrite = kw['OverWrite']
    else:
        OverWrite = False   


    if 'AddFits' in kw:
        AddFits = kw['AddFits']
    else:
        AddFits = True 
        
        
    if 'IsGMOS' in kw:
        IsGMOS = kw['IsGMOS']
    else:
        IsGMOS = True        
    ###### End of argument parsing #######



    ##### NEED TO IMPROVE THE FITTING OF THE FLAT RESPONSE #####

    Flat = []
    for image in image_list:
        if AddFits:
            toopen = RawPath + image + '.fits'
        else:
            toopen = RawPath + image
        print(toopen)
        if IsGMOS:
            hdulist = GMOS_open(toopen)
        else:
            hdulist = fits.open(toopen)
        DFlat = hdulist[0].data-MasterBias
        Bias = DFlat

        index = np.argwhere(np.isnan(DFlat[10,:]))
        mask = np.ones(len(DFlat[10,:]), dtype=bool)
        mask[index] = False        
        xdata = np.arange(DFlat.shape[0])
        flat_1d = np.median(DFlat[int(len(xdata)*0.2):-int(len(xdata)*0.2),:],axis=0)#convolve(np.median(DFlat,axis=0), Box1DKernel(5))
        #spl = UnivariateSpline(xdata[mask], flat_1d[mask], ext=0, k=2 ,s=1000)
        #aa = 10.0**spl(xdata)
        for i in range(len(DFlat)):
            Bias[i,:] = DFlat[i,:]/flat_1d
        
#        index = np.argwhere(np.isnan(DFlat[10,:]))
#        mask = np.ones(len(DFlat[10,:]), dtype=bool)
#        mask[index] = False
 #       for i in range(10,len(DFlat)):
 #           index = np.argwhere(np.isnan(DFlat[i,:]))
 #           mask = np.ones(len(DFlat[i,:]), dtype=bool)
 #           mask[index] = False
#            xdata = np.arange(DFlat.shape[1])
 #           spl = UnivariateSpline(xdata[mask], DFlat[i,mask], ext=0, k=1,s=9999999)
 #          flat_curve = spl(xdata)
 #           aa = savitzky_golay(DFlat[i,:],21,2)
 #           Bias[i,:] = DFlat[i,:]/aa #flat_curve
        Flat.append(Bias) #/np.median((hdulist[0].data[592:2227,:]-MasterBias[592:2227,:])))
    
    
    
    MasterFlat = np.nanmedian(Flat,axis = 0)
    MasterFlat[MasterFlat == np.inf] = 1
    MasterFlat[MasterFlat == -np.inf] = 1
    
    indices = np.argwhere(np.isnan(MasterFlat))
    
    for ind in indices:
        MasterFlat[ind[0],ind[1]] = 1
     
    indices = np.argwhere(MasterFlat == 0)    
    for ind in indices:
        MasterFlat[ind[0],ind[1]] = 1    
    
    # Create a fits file containing the master flat if WriteFile = True   
    hdulist[0].data = MasterFlat
    if WriteFile:
        hdulist.writeto(NameFlat, overwrite = OverWrite)
    
    hdulist.close()
    
    return MasterFlat

def Reduce(image_list,**kw):


    ###### Parse the argument sent through the function #######

    # Do the program displays text during the processing ?
    # [False], True
    if 'Verbose' in kw:
        verbose = kw['Verbose']
    else:
        verbose = False
        
    # Indicate where the raw data are located 
    # Take the current directory as default 
    if 'RawPath' in kw:
        RawPath = kw['RawPath']
    else:
        RawPath = './'   
        
        
    # Check if a bias file is provided
    # Stop the program if not
    if 'Bias' in kw:
        if isinstance(kw['Bias'], str):
            hdulist = fits.open(kw['Bias'])
            MasterBias = hdulist[0].data
        if isinstance(kw['Bias'], np.ndarray):
            MasterBias = kw['Bias']
    else:
        print('****************************')
        print('ERROR: No bias file provided')
        print('Execution stopped')
        return

    # Check if a flat file is provided
    # Stop the program if not
    if 'Flat' in kw:
        if isinstance(kw['Flat'], str):
            hdulist = fits.open(kw['Flat'])
            MasterFlat = hdulist[0].data
        if isinstance(kw['Flat'], np.ndarray):
            MasterFlat = kw['Flat']
    else:
        print('****************************')
        print('ERROR: No flat file provided')
        print('Execution stopped')
        return


    # Do you want to want to overwrite if the file already exist ? 
    # [False]
    if 'OverWrite' in kw:
        OverWrite = kw['OverWrite']
    else:
        OverWrite = False   


    # Do the user provided a list of outputfiles ? 
    if 'OutFilelist' in kw:
        OutList = True
        OutFileList = kw['OutFilelist']
        Out_Suffix = ''
    else:
        OutList = False
        Out_Suffix = '_proc'   

    if 'IsGMOS' in kw:
        IsGMOS = kw['IsGMOS']
    else:
        IsGMOS = True 
    ###### End of argument parsing #######
    
    
    if OutList:
        for ImageIn,ImageOut in zip(image_list,OutFileList):
            toopen = RawPath + ImageIn + '.fits'
            if IsGMOS:
                hdulist = GMOS_open(toopen)
            else:
                hdulist = fits.open(toopen)
            red = (hdulist[0].data-MasterBias)/MasterFlat
            hdulist[0].data = red
            hdulist.writeto(ImageOut + '.fits',overwrite = OverWrite)
    else:
        for ImageIn in image_list:
            toopen = RawPath + ImageIn + '.fits'
            if IsGMOS:
                hdulist = GMOS_open(toopen)
            else:
                hdulist = fits.open(toopen)
            red = (hdulist[0].data-MasterBias)/MasterFlat
            hdulist[0].data = red
            hdulist.writeto(ImageIn + Out_Suffix + '.fits',overwrite = OverWrite)
      


def Wav_Cal2(Arc ,**kw):
    
    ''' need to improve the function to work with any bin value and any R values ''' 
    
    if 'Instrument' in kw:
        Instrument = kw['Instrument']
    else:
        Instrument = 'GMOS'     

    if 'Binning' in kw:
        Binning = kw['Binning']
    else:
        Binning = '2'
        
    if 'Gratting' in kw:
        Gratting = kw['Gratting']
    else:
        Gratting = 'R400+_G5325'  
        
    if 'Detector' in kw:
        Detector = kw['Detector']
    else:
        Detector = 'Det'



    if Instrument == 'Soar':
        print('Instrument = Soar')

        Master = Pipe_Path +'/Soar_Arcs_Master'        

        with open(Master) as f:
            Master_Arcs = f.read().splitlines()
        Master_Arcs = np.array(Master_Arcs).astype(float)

        data = Arc
        
        data = data - np.median(data);
        data = data / np.mean(data);
    
        Lines = np.median(data[250:260,:],axis=0)
        LL = [Lines, Lines,Lines, Lines,Lines, Lines,Lines, Lines]
        LL2 = [Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs]
        
        result = ird.similarity(np.array(LL2), np.array(LL), numiter=3)
        
        x = np.array(range(2071))
        Wav = (x+result['tvec'][1])*1.9752+4956.4



    if Instrument == 'Deveny':
        print('Instrument = Deveny')
        Master = Pipe_Path + '/Deveny_Arcs_Master'
        with open(Master) as f:
            Master_Arcs = f.read().splitlines()
        Master_Arcs = np.array(Master_Arcs).astype(float)
        Master_Arcs = Master_Arcs/39
                
        data = Arc
        
        data = data - np.median(data);
        data = data / np.mean(data);
    
        Lines = np.median(data[250:260,:],axis=0)
        LL = [Lines, Lines,Lines, Lines,Lines, Lines,Lines, Lines]
        LL2 = [Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs]
        
        result = ird.similarity(np.array(LL2), np.array(LL), numiter=3)
        
        x = np.array(range(2148))
        Wav = (x+result['tvec'][1])*-4.277+11780
        
    if Instrument == 'GMOSS' or Instrument == 'GMOSN':
        print(Gratting)
        print(Detector)
        if 'R150+' in Gratting:
            if Binning == '2':
                if 'e2v DD CCD42-90' in Detector:
                    Master = Pipe_Path + '/GeminiN_B2_R150_e2vDD_Arcs_Master'
                    with open(Master) as f:
                        Master_Arcs = f.read().splitlines()
                    
                    Master_Arcs = np.array(Master_Arcs).astype(float)
                    
                    data = Arc
                    data = data - np.nanmedian(data);
                    data = data / np.nanmean(data);
            
                    Lines = np.nanmedian(data[1000:1100,:],axis=0)
            
                    Lines = Lines[1350:1700]
                    Master_Arcs = Master_Arcs[1350:1700]
                    
                else:
                    Master =  Pipe_Path + '/Gemini_B2_R150_Arcs_Master'
                    
                    with open(Master) as f:
                        Master_Arcs = f.read().splitlines()
                    
                    Master_Arcs = np.array(Master_Arcs).astype(float)
                    
                    data = Arc
                    data = data - np.nanmedian(data);
                    data = data / np.nanmean(data);
            
                    Lines = np.nanmedian(data[800:1100,:],axis=0)
            
                    Lines = Lines[500:770]
                    Master_Arcs = Master_Arcs[700:970]

                    Master_Arcs = Master_Arcs/np.nanmax(Master_Arcs)
                    Lines = Lines/np.nanmax(Lines)

                
            if Binning == '4':
                Master = '/Users/maximedevogele/Documents/Gemini/CuAr/Gemini_MasterArc_R150.fits'
                WavCalFile = '/Users/maximedevogele/Documents/Gemini/CuAr/Wav_Cal_MasterArc_R150'                
        if 'R400+' in Gratting:
            if Binning == '2':
                Master = '/Users/maximedevogele/Documents/Gemini/CuAr/Gemini_MasterArc_R400.fits'
#                WavCalFile = '/Users/maximedevogele/Documents/Gemini/CuAr/Wav_Cal_MasterArc_R400'
            if Binning == '4':
                print('Bla')
                if Instrument == 'GMOSS':
                    Master = Pipe_Path + '/Gemini_B4_R400_Arcs_Master'
                    print('master')
                if Instrument == 'GMOSN':
                    if 'e2v DD CCD42-90' in Detector:
                        Master = Pipe_Path + '/GeminiN_B4_R400_e2vDD_Arcs_Master'
                    else: 
                        Master = Pipe_Path + '/Gemini_B4_R400_Arcs_Master'
                        print('master')
                        
                
                with open(Master) as f:
                    Master_Arcs = f.read().splitlines()
                
                Master_Arcs = np.array(Master_Arcs).astype(float)
                
                data = Arc
                data = data - np.nanmedian(data);
                data = data / np.nanmean(data);
                
                Lines = np.nanmedian(data[500:550,:],axis=0)
        
                Lines = Lines[600:1000]
                Master_Arcs = Master_Arcs[600:1000]

        print(len(Lines))
        print(len(Master_Arcs))

        LL = [Lines, Lines,Lines, Lines,Lines, Lines,Lines, Lines]
        LL2 = [Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs]

        print(len(LL))
        print(len(LL2)) 
        
#        for elem in LL[1]:
#            print(elem)
        
        result = ird.similarity(np.array(LL2), np.array(LL), numiter=3)
  
        if Binning == '4':
            x = np.array(range(1562))
            if 'R400+' in Gratting:
                if 'e2v DD CCD42-90' in Detector:
                    Wav = (x+result['tvec'][1])*-2.7632+9152.5
                else:
                    Wav = (x+result['tvec'][1])*-2.9851+9334.3
        if Binning == '2':
            x = np.array(range(3136))
            if 'R150+' in Gratting:
                if 'e2v DD CCD42-90' in Detector:
                    Wav = (x-886+result['tvec'][1])*-3.5461+12505 # -886 because during the extraction we go untill 2250 only and not the full line
                else:
                    print(result['tvec'][1])
                    Wav = (x-886+200+result['tvec'][1])*-3.9237+13149 # -886 because during the extraction we go untill 2250 only and not the full line
                    print(Wav)
        
#        hdulist = fits.open(Master)
#        Master_Arcs = hdulist[0].data
#        hdulist.close()
        
#        WavCal = []
#        with open(WavCalFile) as f:
#             WavCal = [x.split() for x in f.readlines()]
#        WavCal = np.array(WavCal).astype(float)     
        
#        hdulist = fits.open(Arc)
#        data = hdulist[0].data
#        hdulist.close()
        
#        result = ird.similarity(Master_Arcs[990:1000,800:950], data[990:1000,800:950], numiter=3)
#        WavCalSci = copy.deepcopy(WavCal)
        
#        z = np.polyfit((WavCalSci[:,0]-np.mean(WavCalSci[:,0]))/np.std(WavCalSci[:,0]),WavCalSci[:,1], 4)
#        p = np.poly1d(z)
        
#        Xaxis = np.linspace(0,3132,3132)
#        Wav = p((Xaxis-np.mean(WavCalSci[:,0]))/np.std(WavCalSci[:,0]))

    return Wav

def Shift_Spec(Spectre,Err,Wavel,**kw):
    
    if 'Instrument'  in kw:
        Instrument = kw['Instrument']
    else:
        Instrument = 'GMOS'    
    
    
    Spectre = np.array(Spectre)
    Wavel = np.array(Wavel)
    
    Spectre = np.nan_to_num(Spectre)
    
    print(Instrument)
    
    if Instrument == 'Deveny':
        Inter = np.linspace(0.7100,0.7300,1000)
        
        f1 = interp1d(Wavel[0], Spectre[0])
        f2 = interp1d(Wavel[1], Spectre[1])
        
        diffX = []
        diffY = []
        sub = np.linspace(-0.004,0.004,50000)
        for i in sub:
            Inter2 = np.linspace(0.7100+i,0.7300+i,1000)
            New1 = f1(Inter2)
            New2 = f2(Inter)
            dev = np.nanstd(New1/New2)
            diffY.append(dev)
            diffX.append(i)
         
        diffY = np.nan_to_num(diffY)
        diffY[diffY==0] = 100
        plt.figure()
        plt.plot(diffY)
        min_index, min_value = min(enumerate(diffY), key=operator.itemgetter(1))
        Shift = diffX[min_index]
        print(Shift)
        SpecN = f2(Wavel[1][10:-10]-Shift)
        
        SpectreN = [Spectre[0][10:-10],SpecN]
        WavelN = [Wavel[0][10:-10],Wavel[0][10:-10]]
        ErrN = [Err[0][10:-10],Err[1][10:-10]]  
        
    if Instrument == 'GMOS':

        Inter = np.linspace(0.7500,0.7700,1000)
        
        f1 = interp1d(Wavel[0], Spectre[0])
        f2 = interp1d(Wavel[1], Spectre[1])
        
        diffX = []
        diffY = []
        sub = np.linspace(-0.004,0.004,50000)
        for i in sub:
            Inter2 = np.linspace(0.7500+i,0.7700+i,1000)
            New1 = f1(Inter2)
            New2 = f2(Inter)
            dev = np.nanstd(New1/New2)
            diffY.append(dev)
            diffX.append(i)
         
        diffY = np.nan_to_num(diffY)
        diffY[diffY==0] = 100
        plt.figure()
        plt.plot(diffY)
        min_index, min_value = min(enumerate(diffY), key=operator.itemgetter(1))
        Shift = diffX[min_index]
        print(Shift)
        SpecN = f2(Wavel[1][10:-10]-Shift)
        
        SpectreN = [Spectre[0][10:-10],SpecN]
        WavelN = [Wavel[0][10:-10],Wavel[0][10:-10]]   
        ErrN = [Err[0][10:-10],Err[1][10:-10]]  
    
    return SpectreN, ErrN, WavelN;
    
    
    
    


def Bin_Spec(Spectre,Wavel):
    
    
    
    Wave = np.linspace(4800,9100,50)
    Dw = (9100-4800)/50
    
    Wave_out = []
    Spec_out = []
    
    for k in Wave:
        Wave_out.append(k+Dw/2)
        Spec_Tamp = []

        condition1 = Wavel > k
        condition2 = Wavel < k + Dw
        condition = condition1*condition2
        Spec = Spectre
        Spec_Tamp.append(np.nanmedian(np.extract(condition, Spec)))
        Spec_out.append(Spec_Tamp)

    return Wave_out,Spec_out

    
def Divide_Spec(Spectra,Wavel):
    
    Wavel = np.array(Wavel)
    
    WL_max = max(Wavel.flatten())+0.001
    WL_min = min(Wavel.flatten())-0.001
    
    Wave = np.linspace(WL_min,WL_max,1400)
    Dw = (WL_max-WL_min)/1400
    
    Wave_out = []
    Spec_out = []
    
    for k in Wave:
        Wave_out.append(k+Dw/2)
        Spec_Tamp = []
        for i in range(len(Wavel)):
            condition1 = Wavel[i] > k
            condition2 = Wavel[i] < k + Dw
            condition = condition1*condition2
            Spec = Spectra[i]
            Spec_Tamp.append(np.nanmedian(np.extract(condition, Spec)))
        Spec_out.append(Spec_Tamp[0]/Spec_Tamp[1])

    return Wave_out,Spec_out
    


def Normalize_Spectrum(spec,wave,mask,wavelength=7000):
    
    condition1 = wave[mask] > wavelength
    condition2 = wave[mask] < wavelength+100
    condition = condition1*condition2
    Norm_Value = np.extract(condition, spec[mask])
    
    spec = spec/np.nanmedian(Norm_Value)
    
    return spec
    
    


def Average_Spec(Spec_list,Wave_list,Mask_list):

#    WL_max = max(WL.flatten())+0.001
#    WL_min = min(WL.flatten())-0.001
    
    WL_max = 10000
    WL_min = 5000
    
    Wave = np.linspace(WL_min,WL_max,1500)
    Dw = (WL_max-WL_min)/1500
    
    
    Wave_out = []
    Spec_out = []
    for k in Wave:
        Wave_out.append(k+Dw/2)
        Spec_Tamp = []
        for i in range(len(Wave_list)):
            condition1 = Wave_list[i][Mask_list[i]] > k
            condition2 = Wave_list[i][Mask_list[i]] < k + Dw
            condition = condition1*condition2
            Spec = Spec_list[i]
            Spec_Tamp.append(np.nanmedian(np.extract(condition, Spec[Mask_list[i]])))
        Spec_out.append(np.nanmedian(Spec_Tamp))

    return Wave_out, Spec_out


def Sig_Clip_Spec(Spec,n = 7, sig = 4):

    x_size = np.shape(np.array(Spec))[0]
    
    z = np.polyfit(range(x_size), Spec, n)

    p = np.poly1d(z)

    Fitted = p(range(x_size))
    
    SpecS = Sigma_Clip(Fitted-Spec,sig = sig)
    
    return SpecS

def Extract_Wave(name,Trace):
    Wave = Get_Wavtrans(name)
    
    Wavelength = []
    for i in range(1566):
        Wavelength.append(Wave[i,Trace.astype(int)[i]])
    
    return np.array(Wavelength)
    
    
def Detect_Spectra(data,Bin = 4,**kw):
    
    ###### Parse the argument sent through the function #######

    # Does the program displays text during the processing ?
    # [False], True
    if 'Verbose' in kw:
        verbose = kw['Verbose']
    else:
        verbose = False
        
    # Choose the method to detect the spectrum
    # Use the 'Maximum' method as default
    #
    # Methods describtion 
    #   Maximum : Fit a Moffat function around the pixels with highest count.
    #             The region is selected to be around where the spectra should be (center of the image)
    #
    #   Offset : Check the offset of the different acquisition (only if dithering is used). Use an image
    #            which is the sum of all the images and search for the different spectra using the offset
    #           as calibration. 
    if 'Method' in kw:
        Method = kw['Method']
    else:
        Method = 'Maximum'   

    if 'Instrument' in kw:
        Instrument = kw['Instrument']
    else:
        Instrument = 'GMOS'     
    
    if Method == 'Offset':
        if 'Offset' in kw:
            Offset = kw['Offset']
        
        Plate_Scale = 0.08*Bin;
        Offset
        
    
    
    
    if Method == 'Maximum':
        
        if Instrument == 'Deveny' or Instrument == 'Soar':
            
            xs = range(100,400)
            SS = np.median(data[100:400,1400:1430],axis=1)
            
            max_index, max_value = max(enumerate(SS), key=operator.itemgetter(1))
            
            p0 = [0,max(SS.flatten()),2.24,1.41,max_index+100] 
            
            coeff, fit, FWHM, Mask,XOut = Fit_MOFFAT(xs,SS,p0,SClip = False)
            
#            print(coeff)
            
            xs = np.array(xs)
        
#            f, axarr = plt.subplots(2)
#            axarr[0].plot(range(100,400),SS)
#            axarr[0].plot(xs,fit)
            
#            axarr[1].plot(xs,SS)
#            axarr[1].plot(xs,fit)
        
#            axarr[1].text((p0[4])+2*FWHM,p0[1]/2+50,'Center = ' + str(round(p0[4], 3)))
#            axarr[1].text((p0[4])+2*FWHM,p0[1]/2,'FWHM = ' + str(round(FWHM, 3)) + ' pixels' )
#            plt.show()
            return coeff[4]

        
        
        if Instrument == 'GMOSS' or 'GMOSN':
        
            Range = 50
            Trace = Cut_Image(data, Bin = Bin)
            
            Trace = np.array(Trace[0])
            
            print(Trace)
            
            max_index, max_value = max(enumerate(Trace), key=operator.itemgetter(1))
            
            xs = range(max_index-Range,max_index+Range)
            
            print(xs)
            
            p0 = [0,max_value,2.24,1.41,max_index] 
            
            coeff, fit, FWHM, Mask,XOut = Fit_MOFFAT(xs,Trace[xs],p0,SClip = False)
            
            print(coeff)
            
            xs = np.array(xs)
            
            
#            f, axarr = plt.subplots(2)
#            axarr[0].plot(range(In_y,Fin_y),Trace)
#            axarr[0].plot(xs+In_y,fit)
            
#            axarr[1].plot(xs+In_y,Trace[xs])
#            axarr[1].plot(xs+In_y,fit)
        
#            axarr[1].text((p0[4]+In_y)+2*FWHM,p0[1]/2+50,'Center = ' + str(round(p0[4]+In_y, 3)))
#            axarr[1].text((p0[4]+In_y)+2*FWHM,p0[1]/2,'FWHM = ' + str(round(FWHM, 3)) + ' pixels' )
#            plt.show()
            print(coeff[4]+800)
            return coeff[4]+800
    

def Lin_Interp(data_x,data_y,n_bin):
    data_x = np.array(data_x).astype(float)
    data_y = np.array(data_y).astype(float)
    
    New_vec_y = []
    New_vec_x = []
    for idx,elem in enumerate(data_x):
        if idx == 0:
            diff_y = 0
            diff_x = 0
            New_vec_x.append(data_x[idx])
            New_vec_y.append(data_y[idx])
            for elem2 in range(n_bin-1):
                New_vec_y.append(data_y[idx])
                New_vec_x.append(data_x[idx])
        else:
            diff_y = (data_y[idx]-data_y[idx-1])
            diff_x = (data_x[idx]-data_x[idx-1])                
            New_vec_x.append(data_x[idx-1])
            New_vec_y.append(data_y[idx-1]/n_bin)
            for elem2 in range(n_bin-1):
                New_vec_y.append(data_y[idx-1]/n_bin+((diff_y/n_bin))/n_bin*(elem2+1))
                New_vec_x.append(data_x[idx-1]+(diff_x/n_bin)*(elem2+1))

    return New_vec_x,New_vec_y
    


def Fit_Trace(data,Start,Range = 15, SClip = True, **kw):
    
    
    
    if 'Instrument' in kw:
        Instrument = kw['Instrument']
    else:
        Instrument = 'GMOS'   

    if 'Binning' in kw:
        Binning = kw['Binning']
    else:
        Binning = '2'
        
    if 'Gratting' in kw:
        Gratting = kw['Gratting']
    else:
        Gratting = 'R400+_G5325'   
        
    if 'Detector' in kw:
        Detector = kw['Detector']
    else:
        Detector = 'Det'
        
        
    print(Instrument)    


    if Instrument == 'Soar':
    
        x_size = np.shape(data)[1]
        MASK = np.array(range(x_size), dtype=bool)  
        MASK[:20] = False
        MASK[2030:] = False
        
        Trace = np.array(range(x_size), dtype=np.float)
        Trace[:] = 0.
    
        bkg = np.array(range(x_size), dtype=np.float)
        bkg[:] = 0.        
        
        xs = range(int(Start[1])-Range,int(Start[1])+Range)
        SS = data[xs,Start[0]]
        
        p0 = [0,np.max(SS),2.24,1.41,Start[1]]
        
        coeffIn, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0, SClip = SClip)

        Trace[Start[0]] = float(coeffIn[4])

        p0 = coeffIn
        for i in range(Start[0],15,-1):
            if i > 20 :
                xs = range(int(p0[4])-Range,int(p0[4])+Range)
                SS = np.median(data[xs,i-15:i+15],axis=1)
                print(i)
                coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn, SClip = SClip)
                if coeff[4] < 0 or coeff[4] > 2030:
                    coeff = p0
                MASK[i] = Mask
                Trace[i] = float(coeff[4])
                bkg[i] = coeff[0]
                if coeff[4] < 500 and coeff[4] > 12:
                    p0 = coeff

        p0 = coeffIn
        for i in range(Start[0],2030,1):
            xs = range(int(p0[4])-Range,int(p0[4])+Range)
            SS = np.median(data[xs,i-15:i+15],axis=1)
            print(i)
            coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn,SClip = SClip)
            if coeff[4] < 0 or coeff[4] > 2000:
                coeff = p0
            MASK[i] = Mask
            Trace[i] = float(coeff[4])
            bkg[i] = coeff[0]
            p0 = coeff

    #    MASK[CHIP_GAP[1][0]:CHIP_GAP[1][1]] = False
        Xind = np.array(range(x_size))    
        
        Pol = np.polyfit(Xind[MASK], Trace[MASK], 5)
        p = np.poly1d(Pol)
        Tr = p(range(x_size))
        
        plt.figure()
    
        plt.plot(Xind[MASK],Trace[MASK])
        plt.plot(Xind[MASK],Tr[MASK])
        plt.figure()
        plt.plot(Xind[MASK],bkg[MASK])

        return Trace, bkg, MASK 


    
    if Instrument == 'Deveny':
    
        x_size = np.shape(data)[1]
        MASK = np.array(range(x_size), dtype=bool)  
        MASK[:300] = False
        MASK[1840:] = False
        
        Trace = np.array(range(x_size), dtype=np.float)
        Trace[:] = 0.
    
        bkg = np.array(range(x_size), dtype=np.float)
        bkg[:] = 0.        
        
        xs = range(int(Start[1])-Range,int(Start[1])+Range)
        SS = data[xs,Start[0]]
        
        p0 = [0,np.max(SS),2.24,1.41,Start[1]]
        
        coeffIn, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0, SClip = SClip)

        Trace[Start[0]] = float(coeffIn[4])

        p0 = coeffIn
        for i in range(Start[0],15,-1):
            if i > 299 :
                xs = range(int(p0[4])-Range,int(p0[4])+Range)
                SS = np.median(data[xs,i-15:i+15],axis=1)
                print(i)
                coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn, SClip = SClip)
                if coeff[4] < 0 or coeff[4] > 2000:
                    coeff = p0
                MASK[i] = Mask
                Trace[i] = float(coeff[4])
                bkg[i] = coeff[0]
                if coeff[4] < 500 and coeff[4] > 12:
                    p0 = coeff

        p0 = coeffIn
        for i in range(Start[0],1840,1):
            xs = range(int(p0[4])-Range,int(p0[4])+Range)
            SS = np.median(data[xs,i-15:i+15],axis=1)
            coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn,SClip = SClip)
            if coeff[4] < 0 or coeff[4] > 2000:
                coeff = p0
            MASK[i] = Mask
            Trace[i] = float(coeff[4])
            bkg[i] = coeff[0]
            p0 = coeff

    #    MASK[CHIP_GAP[1][0]:CHIP_GAP[1][1]] = False
        Xind = np.array(range(x_size))    
        
        Pol = np.polyfit(Xind[MASK], Trace[MASK], 5)
        p = np.poly1d(Pol)
        Tr = p(range(x_size))
        
        plt.figure()
    
        plt.plot(Xind[MASK],Trace[MASK])
        plt.plot(Xind[MASK],Tr[MASK])
        plt.figure()
        plt.plot(Xind[MASK],bkg[MASK])

        return Trace, bkg, MASK 
        
    if Instrument == 'GMOSS' or Instrument == 'GMOSN':
 
        if Binning == '2':
            MASK = np.array(range(2250), dtype=bool)  
        
            MASK[0:751] = False
            MASK[2058:] = False
            MASK[CHIP_GAP_B2[0][0]:CHIP_GAP_B2[0][1]] = False
            
            Trace = np.array(range(2250), dtype=np.float)
            Trace[:] = 0.

            bkg = np.array(range(2250), dtype=np.float)
            bkg[:] = 0.
            
        if Binning == '4':
            if 'e2v DD CCD42-90' in Detector:
                MASK = np.array(range(1555), dtype=bool) 
                MASK[CHIP_GAP_B4[0][0]:CHIP_GAP_B4[0][1]] = False
        
                Trace = np.array(range(1555), dtype=np.float)
                Trace[:] = 0.
    
                bkg = np.array(range(1555), dtype=np.float)
                bkg[:] = 0. 
            if 'Hamamatsu' in Detector:   
                MASK = np.array(range(1562), dtype=bool) 
                MASK[CHIP_GAP_B4[0][0]:CHIP_GAP_B4[0][1]] = False
        
                Trace = np.array(range(1562), dtype=np.float)
                Trace[:] = 0.
    
                bkg = np.array(range(1562), dtype=np.float)
                bkg[:] = 0.
        
        New_xs = range(int(Start[1])-Range,int(Start[1])+Range)

        New_SS = data[New_xs,Start[0]]
        
        xs, SS = Lin_Interp(New_xs,New_SS,10)

        
        p0 = [0,np.max(SS),2.24,1.41,Start[1]]
        
        
        coeffIn, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0, SClip = SClip)
        
        Trace[Start[0]] = float(coeffIn[4])
        
        if Binning == '2':
            p0 = coeffIn
            for i in range(Start[0],15,-1):
                if i > 1050 or i < 1020 and i > 300 :
                    New_xs = range(int(p0[4])-Range,int(p0[4])+Range)
                    New_SS = np.median(data[New_xs,i-15:i+15],axis=1)
                    xs, SS = Lin_Interp(New_xs,New_SS,10)

                    print(i)
                    coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn, SClip = SClip)
                    if coeff[4] < 0 or coeff[4] > 2000:
                        coeff = p0
                    MASK[i] = Mask
                    Trace[i] = float(coeff[4])
                    bkg[i] = coeff[0]
                    p0 = coeff
             
            p0 = coeffIn
            for i in range(Start[0],2250,1):
                if i > 2100 or i < 2050:
                    print(i)
                    New_xs = range(int(p0[4])-Range,int(p0[4])+Range)
                    New_SS = np.median(data[New_xs,i-15:i+15],axis=1)
                    xs, SS = Lin_Interp(New_xs,New_SS,10)
                    coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn,SClip = SClip)
                    if coeff[4] < 0 or coeff[4] > 2000:
                        coeff = p0
                    MASK[i] = Mask
                    Trace[i] = float(coeff[4])
                    bkg[i] = coeff[0]
                    p0 = coeff
                    
        if Binning == '4':
            p0 = coeffIn
            for i in range(Start[0],15,-1):
                New_xs = range(int(p0[4])-Range,int(p0[4])+Range)  
                New_SS = np.median(data[New_xs,i-15:i+15],axis=1)
                xs, SS = Lin_Interp(New_xs,New_SS,10)
                print(i)
                coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn, SClip = SClip)
                if coeff[4] < 0 or coeff[4] > 2000:
                    coeff = p0
                MASK[i] = Mask
                Trace[i] = float(coeff[4])
                bkg[i] = coeff[0]
                p0 = coeff
             
            p0 = coeffIn
            for i in range(Start[0],1500,1):
                print(i)
                New_xs = range(int(p0[4])-Range,int(p0[4])+Range)
                New_SS = np.median(data[New_xs,i-15:i+15],axis=1)
                xs, SS = Lin_Interp(New_xs,New_SS,10)
                coeff, fit, FWHM, Mask, XOut = Fit_MOFFAT(xs,SS,p0,p_error = coeffIn,SClip = SClip)
                if coeff[4] < 0 or coeff[4] > 1050:
                    coeff = p0
                MASK[i] = Mask
                Trace[i] = float(coeff[4])
                bkg[i] = coeff[0]
                p0 = coeff    
    
    
        if Binning == '2':
            MASK[0:751] = False
            MASK[2058:] = False
            MASK[CHIP_GAP_B2[0][0]:CHIP_GAP_B2[0][1]] = False
        #    MASK[CHIP_GAP[1][0]:CHIP_GAP[1][1]] = False
            Xind = np.array(range(2250))    
            
       #     Pol = np.polyfit(Xind[MASK], Trace[MASK], 3)
       #     p = np.poly1d(Pol)
       #     Tr = p(range(2250))
        if Binning == '4':
            if 'e2v DD CCD42-90' in Detector:
                MASK[CHIP_GAP_B4[0][0]:CHIP_GAP_B4[0][1]] = False
            #    MASK[CHIP_GAP[1][0]:CHIP_GAP[1][1]] = False
                Xind = np.array(range(1555))    
                
                Pol = np.polyfit(Xind[MASK], Trace[MASK], 3)
                p = np.poly1d(Pol)
                Tr = p(range(1555))
            if 'Hamamatsu' in Detector:   
                MASK[CHIP_GAP_B4[0][0]:CHIP_GAP_B4[0][1]] = False
            #    MASK[CHIP_GAP[1][0]:CHIP_GAP[1][1]] = False
                Xind = np.array(range(1562))    
                
                Pol = np.polyfit(Xind[MASK], Trace[MASK], 3)
                p = np.poly1d(Pol)
                Tr = p(range(1562))
            
        plt.figure()
    
        plt.plot(Xind[MASK],Trace[MASK])
        #plt.plot(Xind[MASK],Tr[MASK])
        plt.figure()
        plt.plot(Xind[MASK],bkg[MASK])
        
        return Trace, bkg, MASK


def Extract_Spectrum(data,Trace,bkg,FWHM = 6, Mask = [],**kw):
    
    
    if 'Instrument' in kw:
        Instrument = kw['Instrument']
    else:
        Instrument = 'GMOS' 
        
    if 'Binning' in kw:
        Binning = kw['Binning']
    else:
        Binning = '2'
        
    if 'Gratting' in kw:
        Gratting = kw['Gratting']
    else:
        Gratting = 'R400+_G5325'     
    
    x_size = np.shape(data)[1]   
    
    
    
    if Instrument == 'Deveny' or Instrument == 'Soar':
        
        if len(Mask) == 0:
            Mask = np.ones(x_size, dtype=bool)
        
        
        Spec = []
        for i in range(x_size):
            Sec = data[Trace.astype(int)[i]-FWHM:Trace.astype(int)[i]+FWHM+1,i]- bkg[i]
            Spec.append(sum(Sec))
        
        plt.figure()
        Spec = np.array(Spec)
        plt.plot(Spec[Mask])    
        
        return Spec
            
    
    if Instrument == 'GMOSS' or Instrument == 'GMOSN':
    
        if len(Mask) == 0:
            if Binning == '2':
                Mask = np.ones(2250, dtype=bool)
            if Binning == '4':
                if Instrument == 'GMOSS':
                    Mask = np.ones(1562, dtype=bool)
                if Instrument == 'GMOSN':
                    Mask = np.ones(1562, dtype=bool)
        
        Spec = []
        if Binning == '2':
            for i in range(2250):
                Sec = data[Trace.astype(int)[i]-FWHM:Trace.astype(int)[i]+FWHM+1,i]- bkg[i]
                Spec.append(sum(Sec))
        if Binning == '4':
            if Instrument == 'GMOSS':
                for i in range(1499):
                    if Trace.astype(int)[i] < 10:
                        Trace[i] = 100
                    Sec = data[Trace.astype(int)[i]-(FWHM*2):Trace.astype(int)[i]+(FWHM*2)+1,i]- bkg[i]
                    New_xs, SS = Lin_Interp(range(Trace.astype(int)[i]-(FWHM*2),Trace.astype(int)[i]+(FWHM*2)+1),Sec,10)
                    cond1 = New_xs<Trace[i]+FWHM
                    cond2 = New_xs>Trace[i]-FWHM
                    Cond = cond1*cond2
                    Spec.append(sum(SS*Cond))    
            if Instrument == 'GMOSN':
                for i in range(1499):
                    Sec = data[Trace.astype(int)[i]-FWHM:Trace.astype(int)[i]+FWHM+1,i]- bkg[i]
                    Spec.append(sum(Sec))          
        
#        plt.figure()
#        Spec = np.array(Spec)
#        plt.plot(Spec[Mask])    
        
        return Spec
    
    
    
    




def group_consecutives(vals, step=1):
    """Return list of consecutive lists of numbers from vals (number list)."""
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result


def polyfit2d(x, y, f, deg):
    from numpy.polynomial import polynomial
    import numpy as np
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f)[0]
    return c.reshape(deg+1)


def Sigma_Clip(data, sig = 3):
    
    
    if type(data) == np.ma.core.MaskedArray:
        Cond1 = ma.median(data) - np.std(data)*sig
        Cond2 = ma.median(data) + np.std(data)*sig
        
        Cond_Tamp1 = 0
        Cond_Tamp2 = 0
        
        Counter = 0
        
        while Cond1 != Cond_Tamp1 or Cond2 != Cond_Tamp2:
            
            Counter +=1
            Cond_Tamp1 = Cond1
            Cond_Tamp2 = Cond2
            
            
            data = ma.masked_outside(data, Cond1, Cond2)
            
            Cond1 = ma.median(data) - np.std(data)*sig
            Cond2 = ma.median(data) + np.std(data)*sig

            if Counter > 100:
                return data            
        return data
        
    else:   
        SIG_CLIP = data > 0 
        SIG_CLIP[:] = True
    
        Cond1 = np.median(data[SIG_CLIP]) - np.std(data[SIG_CLIP])*sig
        Cond2 = np.median(data[SIG_CLIP]) + np.std(data[SIG_CLIP])*sig
    
        Cond_Tamp1 = 0
        Cond_Tamp2 = 0
        
        Counter = 0
        while Cond1 != Cond_Tamp1 or Cond2 != Cond_Tamp2:
            
            Counter +=1
            
            Cond_Tamp1 = Cond1
            Cond_Tamp2 = Cond2
            
            SIG_CLIP1 = data > Cond1
            SIG_CLIP2 = data < Cond2
            
            SIG_CLIP = SIG_CLIP1*SIG_CLIP2
            
            Cond1 = np.median(data[SIG_CLIP]) - np.std(data[SIG_CLIP])*sig
            Cond2 = np.median(data[SIG_CLIP]) + np.std(data[SIG_CLIP])*sig
            if Counter > 100:
                return SIG_CLIP
    
        return SIG_CLIP    

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
     import numpy as np
     from math import factorial     
     try:
         window_size = np.abs(np.int(window_size))
         order = np.abs(np.int(order))
     except ValueError:
         raise ValueError("window_size and order have to be of type int")
     if window_size % 2 != 1 or window_size < 1:
         raise TypeError("window_size size must be a positive odd number")
     if window_size < order + 2:
         raise TypeError("window_size is too small for the polynomials order")
     order_range = range(order+1)
     half_window = (window_size -1) // 2
     # precompute coefficients
     b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
     m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
     # pad the signal at the extremes with
     # values taken from the signal itself
     firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
     lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
     y = np.concatenate((firstvals, y, lastvals))
     return np.convolve( m[::-1], y, mode='valid')


def gauss(x, *p):
    A, mu, x0, sigma = p
    return A*np.exp(-((x-x0)-mu)**2/(2.*sigma**2))


def MOFFAT_NEW(x, S0,S1,a, b,x0 , Conv):
    aa = np.zeros(31)
    aa[15-Conv:15+Conv+1] = 1./(Conv*2+1)
    bb = S0+S1/(1+(x-x0)**2/a**2)**b
    cc = np.convolve(bb,aa,'same')
    
    return cc

def MOFFAT(x, *p):
    S0,S1,a, b, x0 = p
    bb = S0+S1/(1+(x-x0)**2/a**2)**b
    
    return bb

def Fit_MOFFAT(xs,ys, p0,p_error = [],SClip = True, Conv = 0):
#    coeff, var_matrix = curve_fit(MOFFAT,xs,ys,p0 )

    xs = np.array(xs)
    ys = np.array(ys)
    
    
    if len(p_error) == 0:
        p_error = p0
        
    try:
        coeff, var_matrix = curve_fit(MOFFAT,xs,ys,p0,maxfev = 100200)
        fit = MOFFAT(xs, *coeff)
        MASK = True
    except Exception:
        coeff = p_error
        fit = []
        MASK = False
        SClip = False
        pass

    ISINF = np.isinf(coeff)
#    fit = MOFFAT(xs, *coeff)
    for d in np.isinf(coeff):
        if d:
            coeff = p_error
            fit = []
            SClip = False
        
    FWHM = 2*coeff[2]*np.sqrt(2**(1/coeff[3])-1)
    
    if SClip == True:
        Res = ys - fit
        ResS = Sigma_Clip(Res,sig =5)
        LT = [i for i, x in enumerate(ResS) if not x]
        if len(LT) != 0:
            coeff, fit, FWHM, MASK, xs = Fit_MOFFAT(xs[ResS],ys[ResS], p0,p_error = [],SClip = True)
    
    
    return coeff, fit, FWHM, MASK, xs

def PolyXY(P,x,y,data):
    return data - P[0]+x*P[1]+y*P[2]+x*y*P[3]
