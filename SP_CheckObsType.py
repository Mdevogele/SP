#!/usr/bin/env

""" SP_ObsTypeCheck - Script that check the observation type for each file and correct it if necessary
    v1.0: 2018-03-15, mdevogele@lowell.edu
    v1.1: 2018-04-24: add support of GMOS-S 
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
import numpy as np
import os
import shutil

from astropy.io import fits

from SP_CheckInstrument import CheckInstrument
from SP_Toolbox import *

from astroquery.jplhorizons import Horizons
from astroquery.simbad import Simbad

import astropy.coordinates as coord
from astropy import units as u

_SP_conf.filenames = []


def CheckObsType(filenames,telescope,obsparam):
    
    # List of usual solar analog ID for astroquery
    
    List_ID_SA = ['HD  11532', 'HD  28099', 'HD 292561', 'BD+00  2717','BD-00  2719','HD 139287',
              'BD+00  3383', 'TYC  447-508-1','BD-00  4074', 'BD-00  4251B','BD-00  4557']
    
    # List of usual solar analog ID for astroquery with linked to usual used designation
    
    List_ID_SA_Comp = {'HD  11532': 'SA93-101', 'HD  28099': 'Hya-64',
                       'HD 292561': 'SA98-978', 'BD+00  2717': 'SA102-1081',
                       'BD-00  2719': 'SA105-56', 'HD 139287': 'SA107-684',
                       'BD+00  3383': 'SA107-998', 'TYC  447-508-1': 'SA110-361',
                       'BD-00  4074': 'SA112-1333','BD-00  4251B': 'SA113-276' ,
                       'BD-00  4557': 'SA115-271'}    
    
#    telescope, obsparam = CheckInstrument(filenames)
    
    
    CollFoc = [];
    
    # Loop over all files to prepare 
    
    for idx, filename in enumerate(filenames):
        
        # If Gemini telescope use special function to open and stick together all the extentions
        try:
            if telescope == 'GMOSS' or telescope == 'GMOSN':
                hdulist = GMOS_open2([filename])
            else:
                hdulist = fits.open(filename ,ignore_missing_end=True)
        except IOError:
            logging.error('cannot open file %s' % filename)
            print('ERROR: cannot open file %s' % filename)
            filenames.pop(idx)
            continue

        
        # Check if the focus header exist (collimator focus), if not put a fake value. 
        # This is need to check for focus run expecialy for Deveny, The header was not present in early images.
        
        try: 
            hdulist[0].header[obsparam['focus']]
        except KeyError:
            hdulist[0].header[obsparam['focus']] = (12, 'SP: dummy focus value') # put dummy value is 'focus' header not found
            
        
        # Update header 
        print('Update header')
        hdulist[0].header['SP_PREPA'] = (str(True),'SP: Has this file been prepared?')
        hdulist[0].header['SP_PREPR'] = (str(False),'SP: Has this file been prepared?')
        hdulist[0].header['SP_COSMI'] = (str(False),'SP: Has this file been prepared?')
        hdulist[0].header['PROCTYPE'] = ('Prepared', 'SP: Type of processed file')
        
        
        
        
        # Still a special treatment for NOT, SHOULD BE FIXED
        if telescope == 'NOT':
            data = hdulist[1].data
            hdulist[0].data = data[800:-100,:].transpose()
            hdulist[1].data = []
        else:
            data = hdulist[0].data

        # Check if the images need to be cropped and transposed
        data = hdulist[0].data
        if obsparam['crop']:
            if obsparam['Y_crop'] and obsparam['X_crop']: # images need to be cropped in both X and Y axes
                data = data[obsparam['Y_crop'][0]:obsparam['Y_crop'][1], obsparam['X_crop'][0]:obsparam['X_crop'][1]]
                
                # Update header
                hdulist[0].header['HISTORY'] = 'The image has been cropped in both X and Y coordinates from :'
                hdulist[0].header['HISTORY'] = 'X axis: ' + str(obsparam['X_crop'][0]) + ' to ' + str(obsparam['X_crop'][1])
                hdulist[0].header['HISTORY'] = 'Y axis: ' + str(obsparam['Y_crop'][0]) + ' to ' + str(obsparam['Y_crop'][1])
                
                
            if obsparam['Y_crop'] and not obsparam['X_crop']: # images need to be cropped in Y axes only
                data = data[obsparam['Y_crop'][0]:obsparam['Y_crop'][1], :]      
                
                # Update header
                hdulist[0].header['HISTORY'] = 'The image has been cropped in coordinates from :'
                hdulist[0].header['HISTORY'] = 'Y axis: ' + str(obsparam['Y_crop'][0]) + ' to ' + str(obsparam['Y_crop'][1])                
                
                
            if obsparam['X_crop'] and not obsparam['Y_crop']: # images need to be cropped in X axes only
                data = data[:, obsparam['X_crop'][0]:obsparam['X_crop'][1]]               
        if obsparam['transpose']:
            data = data.transpose()
            hdulist[0].header['HISTORY'] = 'Image has been rotated'
            
        hdulist[0].data = data # update the data in the fits file


            
        # Get statistic information about the data 
        Med = np.nanmedian(data[10:-10,10:-10])
        Max = np.nanmax(data[10:-10,10:-10])
        std0 = np.nanstd(np.nanmedian(data[10:-10,10:-10],axis=0))
        std1 = np.nanstd(np.nanmedian(data[10:-10,10:-10],axis=1))
        std = np.nanstd(data[10:-10,10:-10])
        Mean =  np.nanmean(data[10:-10,10:-10])
        Median = np.nanmedian(data[10:-10,10:-10])
        
#        print(Mean)
#        print(std1)


        # Update the header with statistical information about the data
        
        hdulist[0].header['Mean'] = (Mean,'Mean of the images (SP)')
        hdulist[0].header['Median'] = (Median,'Median of the images (SP)')
        hdulist[0].header['Std'] = (std0,'Standard deviation of the image (SP)')        
        hdulist[0].header['StdX'] = (std1,'Standard deviation of the image along the x axis (SP)') 
        hdulist[0].header['StdY'] = (std,'Standard deviation of the image along the y axis (SP)') 
        
        
        if telescope =='NOT':
            index = filename.split('.')[0]
            print(index)
            if 'Open' in hdulist[0].header[obsparam['grating']]: # Acquisition files
                TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                Name= telescope + '_' + index + '_' + TIME + '_ACQ_' + EXPTIME + '.fits'
                hdulist.writeto(Name)
                hdulist.close()
#                shutil.copy(filename,Name)
                _SP_conf.filenames.append(Name)
                logging.info('%s changed to %s' % (filename, Name))
                print('%s changed to %s' % (filename, Name))            
            else:
                if std0 > 10000:
                    hdulist[0].header[obsparam['obstype']] = 'FLAT'
                    TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                    EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                    Name= telescope + '_' + index + '_' + TIME + '_FLAT_' + EXPTIME + '.fits'
                    hdulist.writeto(Name)
                    hdulist.close()
#                    shutil.copy(filename,Name)
                    _SP_conf.filenames.append(Name)
                    logging.info('%s changed to %s' % (filename, Name))
                    print('%s changed to %s' % (filename, Name))                

                else: # Biases
                    if std0 + std1 < 20 and int(hdulist[0].header[obsparam['exptime']]) == 0: 
                        hdulist[0].header[obsparam['obstype']] = 'BIAS'
                        TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                        EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                        Name= telescope + '_' + filename.split('.')[1] + '_' + TIME + '_BIAS_' + EXPTIME + '.fits'
                        hdulist.writeto(Name)
                        hdulist.close()
#                        shutil.copy(filename,Name)
                        _SP_conf.filenames.append(Name)
                        logging.info('%s changed to %s' % (filename, Name))
                        print('%s changed to %s' % (filename, Name))          

                    else:
                        
                        hdulist[0].header[obsparam['obstype']] = 'OBJECT'
                        TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                        EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                        OBJECT = hdulist[0].header[obsparam['object']].replace(' ','').replace('/','')
                        TYPE = 'Unknown'
                        if OBJECT.replace(' ',''):
                            try:
                                Horizons(id=OBJECT).ephemerides()
                                TYPE = 'Asteroid'
                            except ValueError:
                                try: 
                                    #result_table = Simbad.query_object(OBJECT)
                                    result_table = Simbad.query_region(coord.SkyCoord(str(hdulist[0].header[obsparam['ra']]) + ' ' + str(hdulist[0].header[obsparam['dec']]) ,unit=(u.deg,u.deg), frame='icrs'), radius='0d1m00s')
                                    T = result_table['MAIN_ID'][0]
                                    print(T)
                                    if result_table['MAIN_ID'][0] in List_ID_SA: 
                                        OBJECT = List_ID_SA_Comp[result_table['MAIN_ID'][0]]
                                        TYPE = 'SA'
                                        print(TYPE)
                                except:
                                    pass
                                if TYPE == 'Unknown':
        #                        except TypeError:
                                    if '(' in OBJECT and ')' in OBJECT:
                                        Number = OBJECT[OBJECT.find("(")+1:OBJECT.find(")")]
                                        if Number.isdigit():
                                            try:
                                                Horizons(id=Number).ephemerides()
                                                TYPE = 'Asteroid'
                                            except ValueError:
                                                TYPE = 'Unknown'
                                    if OBJECT[0:3].isdigit():
                                        if not ' ' in OBJECT:
                                            Name = OBJECT[0:4] + ' ' + OBJECT[4:]
                                            try:
                                                Horizons(id=Name).ephemerides()
                                                TYPE = 'Asteroid'
                                            except ValueError:
                                                TYPE = 'Unknown'
                                                
                        Name= telescope  + '_' + index + '_' + TIME + '_' + TYPE + '_' + OBJECT + '_' + EXPTIME + '.fits'
                        hdulist.writeto(Name,overwrite=True)
#                        shutil.copy(filename,Name)
                                
                        _SP_conf.filenames.append(Name)
                        logging.info('%s changed to %s' % (filename, Name))
                        print('%s changed to %s' % (filename, Name))
                        hdulist.close()  

                
#                print(hdulist[0].header[obsparam['grating']])      

        if telescope == 'SOAR':
            index = filename.split('_')[0]
            
            if hdulist[0].header[obsparam['grating']] == '<NO GRATING>': # Acquisition files
#                print('Acquisition files')
                TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                Name= telescope + '_' + index + '_' + TIME + '_ACQ_' + EXPTIME + '.fits'
                hdulist.close()
#                shutil.copy(filename,Name)
                _SP_conf.filenames.append(Name)
                logging.info('%s changed to %s' % (filename, Name))
                print('%s changed to %s' % (filename, Name))
            else:
                    
                if Med > Max/6 and std0 > 1000:
                    hdulist[0].header[obsparam['obstype']] = 'FLAT'
                    TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                    EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                    Name= telescope + '_' + index + '_' + TIME + '_FLAT_' + EXPTIME + '.fits'
                    hdulist.close()
#                    shutil.copy(filename,Name)
                    _SP_conf.filenames.append(Name)
                    logging.info('%s changed to %s' % (filename, Name))
                    print('%s changed to %s' % (filename, Name))
                
                else: # Biases
                    if std0 + std1 < 10 and int(hdulist[0].header[obsparam['exptime']]) == 0: 
                        hdulist[0].header[obsparam['obstype']] = 'BIAS'
                        TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                        EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                        Name= telescope + '_' + filename.split('.')[1] + '_' + TIME + '_BIAS_' + EXPTIME + '.fits'
                        hdulist.close()
#                        shutil.copy(filename,Name)
                        _SP_conf.filenames.append(Name)
                        logging.info('%s changed to %s' % (filename, Name))
                        print('%s changed to %s' % (filename, Name))                        
 
                    else:
                        hdulist[0].header[obsparam['obstype']] = 'OBJECT'
                        TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                        EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                        OBJECT = hdulist[0].header[obsparam['object']].replace(' ','').replace('/','')
                        TYPE = 'Unknown'
                        if OBJECT.replace(' ',''):
                            try:
                                Horizons(id=OBJECT).ephemerides()
                                TYPE = 'Asteroid'
                            except ValueError:
                                try: 
                                    #result_table = Simbad.query_object(OBJECT)
                                    result_table = Simbad.query_region(coord.SkyCoord(hdulist[0].header[obsparam['ra']] + ' ' + hdulist[0].header[obsparam['dec']] ,unit=(u.hourangle,u.deg), frame='icrs'), radius='0d1m00s')
                                    T = result_table['MAIN_ID'][0]
                                    if result_table['MAIN_ID'][0] in List_ID_SA: 
                                        OBJECT = List_ID_SA_Comp[result_table['MAIN_ID'][0]]
                                        TYPE = 'SA'
                                except:
                                    pass
                                if TYPE == 'Unknown':
        #                        except TypeError:
                                    if '(' in OBJECT and ')' in OBJECT:
                                        Number = OBJECT[OBJECT.find("(")+1:OBJECT.find(")")]
                                        if Number.isdigit():
                                            try:
                                                Horizons(id=Number).ephemerides()
                                                TYPE = 'Asteroid'
                                            except ValueError:
                                                TYPE = 'Unknown'
                                    if OBJECT[0:3].isdigit():
                                        if not ' ' in OBJECT:
                                            Name = OBJECT[0:4] + ' ' + OBJECT[4:]
                                            try:
                                                Horizons(id=Name).ephemerides()
                                                TYPE = 'Asteroid'
                                            except ValueError:
                                                TYPE = 'Unknown'
                        if std0/std1 > 50:
                            hdulist[0].header[obsparam['obstype']] = 'ARCS/FOCUS'                                           
                            TYPE = 'ARCS'
                            OBJECT = filename.split('_')[1]
                                            
                                    
                        Name= telescope  + '_' + index + '_' + TIME + '_' + TYPE + '_' + OBJECT + '_' + EXPTIME + '.fits'
    
#                        shutil.copy(filename,Name)
            hdulist[0].writeto(Name,overwrite=True)                        
            _SP_conf.filenames.append(Name)
            logging.info('%s changed to %s' % (filename, Name))
            print('%s changed to %s' % (filename, Name))
            hdulist.close()    
                       
        
        
        ######################################################################
        #                       Prepare DEVENY DATA
        ######################################################################
        
        if telescope == 'DEVENY':
            Arcs = False
            
            # Extract informations for header 
            
            TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
            EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')

            
            # get the index of the file
            index = filename.split('.')[1]
            
            # Identify FLAT images 
            if Med > Max/5 and std0 > 1000:
                # Update header
                hdulist[0].header['HISTORY'] = 'SP: %s changed to FLAT' % (str(hdulist[0].header[obsparam['obstype']]))
                hdulist[0].header[obsparam['obstype']] = 'FLAT'
                Name= telescope + '_' + index + '_' + TIME + '_FLAT_' + EXPTIME + '.fits'
                hdulist[0].header['HISTORY'] = 'SP: %s changed to %s' % (filename, Name)
                
            else:
                # Identify BIAS images
                if std0 + std1 < 10 and int(hdulist[0].header[obsparam['exptime']]) == 0: 
                    
                    # Update header
                    hdulist[0].header['HISTORY'] = 'SP: %s changed to BIAS' % (str(hdulist[0].header[obsparam['obstype']]))
                    hdulist[0].header[obsparam['obstype']] = 'BIAS'
                    Name= telescope + '_' + filename.split('.')[1] + '_' + TIME + '_BIAS_' + EXPTIME + '.fits'
                    hdulist[0].header['HISTORY'] = 'SP: %s changed to %s' % (filename, Name)
                    
                else:
                    # Identify ARCS images
                    if std0/std1 > 100:
                        
                        # Keep track of the focus value for arcs images
                        CollFoc.append((int(filename.split('.')[1]),hdulist[0].header['COLLFOC'],idx))                                                            
                        
#                        Name= telescope + '_' + filename.split('.')[1] + '_' + TIME + '_ARCS_' + EXPTIME + '.fits'

                        # Update header 
                        hdulist[0].header['HISTORY'] = 'SP: %s changed to ARCS/FOCUS' % (str(hdulist[0].header[obsparam['obstype']]))
                        hdulist[0].header[obsparam['obstype']] = 'ARCS/FOCUS' 
                        Arcs = True
                    # If not Flat, Bias, or Arcs then is a science acquisition            
                    else:
                        
                        # Update header                         
                        hdulist[0].header['HISTORY'] = 'SP: %s changed to OBJECT' % (str(hdulist[0].header[obsparam['obstype']]))
                        hdulist[0].header[obsparam['obstype']] = 'OBJECT'
                        
                        OBJECT = hdulist[0].header[obsparam['object']].replace(' ','').replace('/','')
                        
                        # Assume an unknown object as a first gest 
                        TYPE = 'Unknown'
                        print(obsparam['object'])
                        if OBJECT.replace(' ',''):
                            print(OBJECT)
                            try: # Try to retrieve the object using Horizons, if it find the object, assigned the object type to Asteroid
                                Horizons(id=OBJECT).ephemerides()
                                TYPE = 'Asteroid'
                            except ValueError: #if object not find in Horizons, check simbad for standard stars
                                try: 
                                    result_table = Simbad.query_region(coord.SkyCoord(hdulist[0].header[obsparam['ra']] + ' ' + hdulist[0].header[obsparam['dec']] ,unit=(u.hourangle,u.deg), frame='icrs'), radius='0d1m00s')
                                    T = result_table['MAIN_ID'][0]
                                    if result_table['MAIN_ID'][0] in List_ID_SA: # Compares the ID found with simbad with our list of standard stars, if there is a match, Type = SA
                                        OBJECT = List_ID_SA_Comp[result_table['MAIN_ID'][0]]
                                        TYPE = 'SA'
                                except: # If not found then object remains unknown
                                    pass
                                if TYPE == 'Unknown': #if object still unknown try to clean the object name
        #                        except TypeError:
                                    if '(' in OBJECT and ')' in OBJECT:
                                        Number = OBJECT[OBJECT.find("(")+1:OBJECT.find(")")]
                                        if Number.isdigit():
                                            try:
                                                Horizons(id=Number).ephemerides()
                                                TYPE = 'Asteroid'
                                            except ValueError:
                                                TYPE = 'Unknown'
                                    if OBJECT[0:3].isdigit():
                                        if not ' ' in OBJECT:
                                            Name = OBJECT[0:4] + ' ' + OBJECT[4:]
                                            try:
                                                Horizons(id=Name).ephemerides()
                                                TYPE = 'Asteroid'
                                            except ValueError:
                                                TYPE = 'Unknown'
                                 
                                
                                        
                                
                            Name= telescope  + '_' + index + '_' + TIME + '_' + TYPE + '_' + OBJECT + '_' + EXPTIME + '.fits'
                            hdulist[0].header['HISTORY'] = 'SP: %s changed to %s' % (filename, Name)

            if not Arcs:
                hdulist[0].writeto(Name,overwrite=True)
                _SP_conf.filenames.append(Name)
                print('%s changed to %s' % (filename, Name))
                logging.info('%s changed to %s' % (filename, Name))
                hdulist.close()        
        
        
        
        if telescope == 'GMOSS':
            index = filename.split('S')[2].replace('.fits','')
            TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
            EXPTIME = "{0:.1f}".format(round(hdulist[0].header['DARKTIME'],2)).replace('.','s')
            TYPE = hdulist[0].header[obsparam['obsclass']]            
            CENTWAVE  = str(int(hdulist[0].header['CENTWAVE']))
            OBJECT = hdulist[0].header[obsparam['object']].replace(' ','').replace('/','')
            Name= telescope + '_' + index + '_' + TIME + '_' + TYPE + '_' + OBJECT + '_' + CENTWAVE +'_' + EXPTIME + '.fits'
            hdulist.writeto(Name)
         
        if telescope == 'GMOSN':
            index = filename.split('S')[1].replace('.fits','')
            TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
            EXPTIME = "{0:.1f}".format(round(hdulist[0].header['DARKTIME'],2)).replace('.','s')
            TYPE = hdulist[0].header[obsparam['obsclass']]            
            CENTWAVE  = str(int(hdulist[0].header['CENTWAVE']))
            OBJECT = hdulist[0].header[obsparam['object']].replace(' ','').replace('/','')
            Name= telescope + '_' + index + '_' + TIME + '_' + TYPE + '_' + OBJECT + '_' + CENTWAVE +'_' + EXPTIME + '.fits'
            hdulist.writeto(Name)    
            
        

    if telescope == 'DEVENY':
        
        # In order to detect focus run and arcs run. During a focus run the focus 
        # is changing while is stays the same during a calibration run             
        Cons = Get_Consecutive(np.array(CollFoc)[:,0].astype(int))
        Focus = []
        Series = []
        for idx, elem in enumerate(Cons):
            for idx2, elem2 in enumerate(elem):
                Focus.append(CollFoc[Cons[idx][idx2][0]][1])
            if len(set(Focus)) == 1:
                Series.append('ARCS')
            else:
                Series.append('FOCUS')  
            Focus = []                 
                
        for idx, elem in enumerate(Cons):
            for idx2, elem2 in enumerate(elem):    
                try:
                    hdulist = fits.open(filenames[CollFoc[Cons[idx][idx2][0]][2]], mode='update' ,ignore_missing_end=True)
                except IOError:
                    logging.error('cannot open file %s' % filename)
                    print('ERROR: cannot open file %s' % filename)
                    filenames.pop(idx)
                    continue
                if Series[idx] == 'ARCS':
                    hdulist[0].header[obsparam['obstype']] = 'ARCS'
                    TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                    EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                    Name= telescope + '_' + filenames[CollFoc[Cons[idx][idx2][0]][2]].split('.')[1] + '_' + TIME + '_ARCS_' + '_' + EXPTIME + '.fits'
                    hdulist.close()
                    shutil.copy(filenames[CollFoc[Cons[idx][idx2][0]][2]],Name)
                    _SP_conf.filenames.append(Name)
                    logging.info('%s changed to %s' % (filenames[CollFoc[Cons[idx][idx2][0]][2]], Name))
                    print('%s changed to %s' % (filenames[CollFoc[Cons[idx][idx2][0]][2]], Name))
                else:
                    hdulist[0].header[obsparam['obstype']] = 'FOCUS'
                    TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').replace('.','')
                    EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
                    Name= telescope + '_' + filenames[CollFoc[Cons[idx][idx2][0]][2]].split('.')[1] + '_' +  TIME + '_FOCUS_' + '_' + EXPTIME + '.fits'
                    hdulist.close()
                    shutil.copy(filenames[CollFoc[Cons[idx][idx2][0]][2]],Name)
                    _SP_conf.filenames.append(Name)
                    logging.info('%s changed to %s' % (filenames[CollFoc[Cons[idx][idx2][0]][2]], Name))
                    print('%s changed to %s' % (filenames[CollFoc[Cons[idx][idx2][0]][2]], Name))
                
#    for elem in filenames:
#        os.remove(elem)                  
       
    return None
    
    