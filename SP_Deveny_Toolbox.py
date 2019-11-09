#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 15:43:21 2018

@author: maximedevogele

Deveny file processing

"""

import glob
from astropy.io import fits
import numpy as np
import SP_Toolbox as SP
import callhorizons
import ephem
from astropy.time import Time
from itertools import groupby
from operator import itemgetter

def Detect_Focus(ArcsList):
    Arcs_Number= []
    for ele in ArcsList:
        Arcs_Number.append(int(ele.split('.')[2]))
        
    Arcs_Group = []
    for k, g in groupby(enumerate(Arcs_Number), lambda (i, x): i-x):
        Arcs_Group.append(map(itemgetter(1), g))

    IsFocus = []
    for List in Arcs_Group:
        Focus = []
        for ele in List:
            files = glob.glob('*' + str(ele) + '.fits')
            hdulist = fits.open(files[0])
            Focus.append(hdulist[0].header['COLLFOC'])
            hdulist.close()
        if np.std(Focus) > 0.1:
            IsFocus.append(True)
        else:
            IsFocus.append(False)

    f = open("Fix_Header_Log",'a')
    d = open("Focus.fl",'w')
    s = open("Arcs.fl",'w')
    for List, IsF in zip(Arcs_Group,IsFocus): # Fix the header 
        for ele in List:
            files = glob.glob('*' + str(ele) + '.fits')
            hdulist = fits.open(files[0])
            if IsF:
                hdulist[0].header['OBSTYPE'] = 'COMPARISON'
                hdulist[0].header['IMAGETYP'] = 'COMPARISON'
                hdulist[0].header['OBJECT'] = 'Focus'
                hdulist[0].header['OBJNAME'] = 'Focus'
                f.write(files[0] + "\t COMPARISON \t is a focus sequence \t Header changed to Focus \n" )
                d.write(files[0] + '\n')
            else:                          
                hdulist[0].header['OBSTYPE'] = 'COMPARISON'
                hdulist[0].header['IMAGETYP'] = 'COMPARISON'
                hdulist[0].header['OBJECT'] = 'Arcs'
                hdulist[0].header['OBJNAME'] = 'Arcs' 
                f.write(files[0] + "\t COMPARISON \t is a arcs sequence \t Header changed to Arcs \n" )
                s.write(files[0] + '\n')
            hdulist.writeto(files[0],overwrite=True)
            hdulist.close()
    f.close()
    d.close()
    s.close()
    
            

def decdeg2dms(dd):
    negative = dd < 0
    dd = abs(dd)
    minutes,seconds = divmod(dd*3600,60)
    degrees,minutes = divmod(minutes,60)
    if negative:
        if degrees > 0:
            degrees = -degrees
        elif minutes > 0:
            minutes = -minutes
        else:
            seconds = -seconds
    return (degrees,minutes,seconds)

def RAdms2RAdeg(Deg,Min,Sec):

    Dec = (float(Deg)+float(Min)/60+float(Sec)/3600)*15

    return Dec


def Date2J2000(hdulist):
    
    T = ephem.Equatorial(hdulist[0].header['RA'],hdulist[0].header['DEC'], epoch= str(hdulist[0].header['EQUINOX']) )
    S = ephem.Equatorial(T, epoch='2000')

    return (S.ra*180/np.pi, S.dec*180/np.pi)

def Fix_Header():
    Files = glob.glob('./*.fits')
    f = open("Fix_Header_Log",'w')
    
    for files in Files:
        hdulist = fits.open(files)
        data = hdulist[0].data
        Type = hdulist[0].header['OBSTYPE']
        if hdulist[0].header['EXPTIME'] == 0:
            print(files + "\t" + Type + "\t is exposure of 0s \t Header changed to Bias \n")                
            f.write(files + "\t" + Type + "\t is exposure of 0s \t Header changed to Bias \n")
            hdulist[0].header['OBSTYPE'] = 'Bias'
            hdulist[0].header['IMAGETYP'] = 'Bias'
            hdulist[0].header['OBJNAME'] = 'Bias'
            hdulist[0].header['OBJECT'] = 'Bias'
            hdulist.writeto(files,overwrite=True)            
        else:    
            if hdulist[0].header['OBSTYPE'] == 'COMPARISON':
                if hdulist[0].header['LAMPCAL'] == '':
                    print("The lamps are off, this could no be a COMPARISON file")
                    if np.median(data) > 5000:
                        print("Probably a flat file")
                        hdulist[0].header['OBSTYPE'] = 'Flat'
                        hdulist[0].header['IMAGETYP'] = 'Flat'
                        hdulist[0].header['OBJECT'] = 'Flat'
                        hdulist[0].header['OBJNAME'] = 'Flat' 
                        f.write(files + "\t COMPARISON \t All lamps are switched off \t Probaly flat file Header changed to Flat \n")
                    else:
                        hdulist[0].header['OBSTYPE'] = 'Unknown'
                        hdulist[0].header['IMAGETYP'] = 'Unknown'
                        hdulist[0].header['OBJECT'] = 'Unknown'
                        hdulist[0].header['OBJNAME'] = 'Unknown' 
                        f.write(files + "\t COMPARISON \t All lamps are switched off \t Header changed to Unknown \n" )
            else:
                if not hdulist[0].header['LAMPCAL'] == '':
                    Type = hdulist[0].header['OBSTYPE']
                    if np.median(data) > 5000:
                        hdulist[0].header['OBSTYPE'] = 'Wrong FL'
                        hdulist[0].header['IMAGETYP'] = 'Wrong FL'
                        hdulist[0].header['OBJECT'] = 'Wrong FL'
                        hdulist[0].header['OBJNAME'] = 'Wrong FL'                 
                    else:
                        hdulist[0].header['OBSTYPE'] = 'COMPARISON'
                        hdulist[0].header['IMAGETYP'] = 'COMPARISON'
                        hdulist[0].header['OBJECT'] = 'arcs'
                        hdulist[0].header['OBJNAME'] = 'arcs'                
                    f.write(files + "\t" + Type + "\t Lamps are switched on \t Header changed to COMPARISON and arcs \n" )
                    hdulist.writeto(files,overwrite=True)
                    
            
#            if hdulist[0].header['OBSTYPE'] == 'BIAS':
                    
#                Median = np.median(data)
#                std = np.std(data)
#                if hdulist[0].header['LAMPCAL'] is not '':
#                    f.write(files + "\t" + Type + "\t Lamps are switched on \t Header changed to COMPARISON and arcs \n" )
#                    hdulist[0].header['OBSTYPE'] = 'COMPARISON'
#                    hdulist[0].header['IMAGETYP'] = 'COMPARISON'
#                    hdulist[0].header['OBJECT'] = 'arcs'
#                    hdulist[0].header['OBJNAME'] = 'arcs'                  
#                if Median < 2000 or Median > 2800 or std > 100:
#                    print("This file does not seems to be a BIAS frame")
#                    hdulist[0].header['OBSTYPE'] = 'Unknown'
#                    hdulist[0].header['IMAGETYP'] = 'Unknown'
#                    hdulist[0].header['OBJECT'] = 'Unknown'
#                    hdulist[0].header['OBJNAME'] = 'Unknown' 
#                    f.write(files + "\t BIAS \t the median values or std is out of range (" + str(Median) + " " + str(std) + ") \t Header changed to Unknown \n")
#                else:
#                    hdulist[0].header['OBJNAME'] = 'Bias'
#                    hdulist[0].header['OBJECT'] = 'Bias'
#                    hdulist.writeto(files,overwrite=True)
#                    f.write(files + "\t BIAS \t everything OK \n")
                    
            if hdulist[0].header['OBSTYPE'] == 'OBJECT':
                Median = np.median(data)
                std = np.std(data)
                Type = hdulist[0].header['OBSTYPE']
                if hdulist[0].header['LAMPCAL'] is not '' or hdulist[0].header['OBJNAME'] == 'arcs' or hdulist[0].header['OBJNAME'] == 'Arcs':
                    f.write(files + "\t" + Type + "\t Lamps are switched on \t Header changed to COMPARISON and arcs \n" )
                    hdulist[0].header['OBSTYPE'] = 'COMPARISON'
                    hdulist[0].header['IMAGETYP'] = 'COMPARISON'
                    hdulist[0].header['OBJECT'] = 'Arcs'
                    hdulist[0].header['OBJNAME'] = 'Arcs'
                    hdulist.writeto(files,overwrite=True)
                    print('Test')
                                        
                if hdulist[0].header['OBJNAME'] == 'focus' or hdulist[0].header['OBJNAME'] == 'Focus':
                    f.write(files + "\t" + Type + "\t Lamps are switched on \t Header changed to focus \n" )
                    hdulist[0].header['OBSTYPE'] = 'focus'
                    hdulist[0].header['IMAGETYP'] = 'focus'
                    hdulist[0].header['OBJECT'] = 'focus'
                    hdulist[0].header['OBJNAME'] = 'focus'
                    hdulist.writeto(files,overwrite=True)
                    
                
                '''
                else:
                    ObjName = hdulist[0].header['SCITARG']   
                    t = Time(hdulist[0].header['DATE'], format='isot', scale='utc')
                    call = callhorizons.query(ObjName)
                    call.set_discreteepochs([str(t.jd)])
                    try:
                        call.get_ephemerides('G37')
                        Asteroid = True
                    except ValueError:
                        print('Not an asteroid')
                        Asteroid = False
                    if Asteroid == True:
                        CoordsFits = np.array(Date2J2000(hdulist))
                        CoordsHor = np.array((float(call['RA']), float(call['DEC'])))
                        Dist = (ephem.separation(CoordsFits*np.pi/180,CoordsHor*np.pi/180)*180/np.pi)*3600
                        if Dist > 60:
                            print('bla')
        #                f.write(files + "\t OBJECT \t The object ("+ ObjName + ") is an asteroid \n")
                    else:
     #                   Res = Simbad.query_object(ObjName)
                        f.write(files + "\t OBJECT \t The nature of the object ("+ ObjName + ") is a star \n")
                        if Res == None:
                            f.write(files + "\t OBJECT \t The nature of the object ("+ ObjName + ") is unknown \n")
                '''
    f.close()

def Create_FileLists():

   Files = glob.glob('./*.fits')
   
   for files in Files:
       hdulist = fits.open(files)
#       if hdulist[0].header['IMAGETYP'] == 'BIAS':
#           Bias.append(files)
#       if hdulist[0].header['IMAGETYP'] == 'COMPARISON':
#           f = open(hdulist[0].header['OBJNAME'] + '.fl','a')
#           f.write(files + '\n')
#           f.close()
#       if hdulist[0].header['IMAGETYP'] == 'OBJECT':
#           Object.append(files)
#       if 'FLAT' in hdulist[0].header['IMAGETYP']
       
       if hdulist[0].header['OBSTYPE'] == 'OBJECT':
           Name = hdulist[0].header['SCITARG']
           Type = "Unknown"
           if Name[0:4].isdigit() or Name[0] == 'a':
               Type = "Ast"
           if Name[0:2] == "SA" or Name[0:3] == 'Hya':
               Type = "SA"
           if Type == "Unknown":
               print("SP could not determine if " + Name + " is an asteroid or a solar analog")
               print("SP could not determine if " + files + " contains an asteroid or a solar analog")
               Type = raw_input('Please enter "Ast" if this object is an asteroid or "SA" if it is a solar analog: ')
        
           f = open(Type + '_' + hdulist[0].header['SCITARG'] + '.fl','a')
           f.write(files + '\n')
           f.close()
           hdulist.close()
       elif "FLAT" in hdulist[0].header['OBSTYPE']:
           Target = hdulist[0].header['SCITARG']
           f = open(hdulist[0].header['OBSTYPE'] + '_' + Target +'.fl','a')
           f.write(files + '\n')
           f.close()
           hdulist.close() 
       else:
           f = open(hdulist[0].header['OBJNAME'] + '.fl','a')
           f.write(files + '\n')
           f.close()
           hdulist.close()          
       
def Get_BiasList(filename = 'Bias.fl'):

    try:
        with open(filename,'r') as f:
            List = f.read().splitlines()
    except IOError:
        print 'The file "' + filename + '" does not exist.'
        NameFile = raw_input("Please enter name of the file containing the Bias file list: ")
        List = Get_FileList(NameFile)
    
    return List

def Get_FlatList(filename = 'DomeFlats.fl'):

    try:
        with open(filename,'r') as f:
            List = f.read().splitlines()
    except IOError:
        print 'The file "' + filename + '" does not exist.'
        NameFile = raw_input("Please enter name of the file containing the Flat file list: ")
        List = Get_FileList(NameFile)
    
    return List

def Get_ArcsList(filename = 'Arcs.fl'):

    try:
        with open(filename,'r') as f:
            List = f.read().splitlines()
    except IOError:
        print 'The file "' + filename + '" does not exist.'
        NameFile = raw_input("Please enter name of the file containing the Arcs file list: ")
        List = Get_FileList(NameFile)
    
    return List


def Get_ObjectList(filename):

    try:
        with open(filename,'r') as f:
            List = f.read().splitlines()
    except IOError:
        print 'The file "' + filename + '" does not exist.'
        NameFile = raw_input("Please enter name of the file containing the Object file list: ")
        List = Get_FileList(NameFile)
    
    return List


def Get_FileList(FileName):

    try:
        with open(FileName,'r') as f:
            List = f.read().splitlines()
    except IOError:
        print 'The file "' + FileName + '" does not exist.'
        raise
    
    return List

def Get_ObservedAsteroids():

    Files = glob.glob('./*.fits')
    AstList = [] 
    for files in Files:
       hdulist = fits.open(files)

       
       if hdulist[0].header['OBSTYPE'] == 'OBJECT':
           Name = hdulist[0].header['SCITARG']
           if Name[0:4].isdigit() or Name[0] == 'a':
               AstList.append(Name)

    return set(AstList)

def Get_FileNumbers(FileName):

    
    List = Get_FileList(FileName)
    NumberList = []
    for files in List:
        NumberList.append(int(files.split('.')[-2]))
    
    return NumberList
        
    

def Get_ObservedSolAnal():

    Files = glob.glob('./*.fits')
    SAList = [] 
    for files in Files:
       hdulist = fits.open(files)

       
       if hdulist[0].header['OBSTYPE'] == 'OBJECT':
           Name = hdulist[0].header['OBJNAME']
           if Name[0:2] == "SA" or Name[0:3] == 'Hya':
               SAList.append(Name)

    return set(SAList)



def Bckg_Sub(FileName, Area = [250,350] ):
    hdulist = fits.open(FileName)
    image = hdulist[0].data
        
    Mask = np.zeros((516,2148))
    Mask[:,:] =  0
    
    Mask[0:Area[0],:] = 1
    Mask[Area[1]:,:] = 1
    Mask = Mask.astype(bool)
    
    image2 = np.ma.masked_array(image,mask = Mask)
    
    X = range(516)
    X = np.array(X)
        
    for i in range(2148):
        image3 = SP.Sigma_Clip(image2[:,i])
        XX = np.ma.masked_array(X,image3.mask)
        z = np.ma.polyfit(XX,image3 , 1)
        p = np.poly1d(z)
        image[:,i] = image[:,i] - p(X)
    
    hdulist.writeto( FileName.replace('.fits','') + '.bkgSub' + '.fits')

def Appaired_Ast_SA(ObsSum):
    Num1 = []
    for Ast in ObsSum['Ast']:
        Num1.append((ObsSum['Ast'][Ast][0],Ast))
    Sort_Ast = sorted(Num1)
    Num2 = []
    for SA in ObsSum['SA']:
        Num2.append((ObsSum['SA'][SA][0],SA))    
    Sort_SA = sorted(Num2)
    
    PairedSum=[]
    for Ast,SA in zip(Sort_Ast,Sort_SA):
        PairedSum.append((Ast[1],SA[1]))
        
    return PairedSum

def Get_Paired_Obj(Pair,Object):
    Ass = ''
    for tup in Pair:
        T = str(Object) in tup
        if T:
            if tup[0] == Object:
                Ass = tup[1]
            else:
                Ass = tup[0]
    
    return Ass
                    
    