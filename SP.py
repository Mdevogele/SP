#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:31:15 2019

@author: maximedevogele
"""

from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askopenfilenames
from tkinter.messagebox import showerror

from astropy.io import fits

import logging
import datetime

from past.utils import old_div


from PIL import Image
from PIL import ImageTk
from PIL import ImageDraw

import argparse, shlex
import numpy as np

from scipy.ndimage import interpolation as interp

from scipy import misc


import time
from datetime import date
from datetime import timedelta

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os

import _SP_conf
Pipe_Path = _SP_conf.rootpath
module_logger = logging.getLogger(__name__ + '.A')
module_logger2 = logging.getLogger(__name__ + '.B')
module_logger_Bias = logging.getLogger(__name__ + '.C')
module_logger_Flat = logging.getLogger(__name__ + '.D')
module_logger_Preproc = logging.getLogger(__name__ + '.E')
module_logger_CosmCorr = logging.getLogger(__name__ + '.F')
module_logger_Bckgsub = logging.getLogger(__name__ + '.G')
module_logger_Extract = logging.getLogger(__name__ + '.H')
module_logger_WavCal = logging.getLogger(__name__ + '.I')



class simpleapp_tk(Tk):
    def __init__(self,parent):
        Tk.__init__(self,parent)
        
        self.zoom = 0.5
        self.parent = parent

        self.grid()
        
        self.frame_process=Frame(self, width=1000, height=400)
        self.frame_process.grid(column=0, row=0)

        self.frame_fits=Frame(self, width=1000, height=200)
        self.frame_fits.grid(column=0, row=1)
        
        self.TabControl = ttk.Notebook(self.frame_process)
        self.tab1 = Frame(self.TabControl)
        self.TabControl.add(self.tab1,text='Analysis')
        
        self.tab_Prepare = Frame(self.TabControl)
        self.TabControl.add(self.tab_Prepare,text='Prepare')

        self.tab_Bias = Frame(self.TabControl)
        self.TabControl.add(self.tab_Bias,text='Bias')

        self.tab_Flat = Frame(self.TabControl)
        self.TabControl.add(self.tab_Flat,text='Flat')

        self.tab_Preproc = Frame(self.TabControl)
        self.TabControl.add(self.tab_Preproc,text='Pre processing')

        self.tab_CosmCorr = Frame(self.TabControl)
        self.TabControl.add(self.tab_CosmCorr,text='Cosmic correction')

        self.tab_BckgSub = Frame(self.TabControl)
        self.TabControl.add(self.tab_BckgSub,text='Background subtraction')

        self.tab_Extract = Frame(self.TabControl)
        self.TabControl.add(self.tab_Extract,text='Spectra extraction')        

        self.tab_WavCal = Frame(self.TabControl)
        self.TabControl.add(self.tab_WavCal,text='Wavelength calibration')   
        
        self.TabControl.pack(expand=1, fill="both") 
        

        ### load dummy image for canvas 
        
        Manos_img= misc.imread(Pipe_Path + '/manos_splash.png')
        self.tkimage = ImageTk.PhotoImage(Image.fromarray(Manos_img), palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  

        im = Image.fromarray(Manos_img)
        self.tkimage.paste(im)


        self.LF_button = Button(self.tab1, text="Load files", width=10)
        self.LF_button.grid(row=1, column=0, sticky=W)
        self.LF_button.bind("<ButtonRelease-1>", self.load_file)        
        
        self.Plot_button = Button(self.tab1, text="Plot Spectrum", command=self.Plot_Spectrum, width=20)
        self.Plot_button.grid(row=2, column=0, sticky=W)   
        
        self.Fits_button = Button(self.tab1, text="Display fits", command=self.Display_fits, width=20)
        self.Fits_button.grid(row=3, column=0, sticky=W)   

        self.frame=Frame(self.tab1, width=1000, height=400)
        self.frame.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame.grid_propagate(False)

        self.frame.grid_rowconfigure(0, weight=1)
        self.frame.grid_columnconfigure(0, weight=1)

        self.frame.mytext2 = Text(self.frame, state="disabled")
        self.frame.mytext2.place(x=10, y=10, height=990, width=390)

        scrollb = Scrollbar(self.frame, command=self.frame.mytext2.yview)
        scrollb.grid(row=0, column=1, sticky='nsew')
        self.frame.mytext2['yscrollcommand'] = scrollb.set
#
#        # Running the pipeline:
#          
        self.Prepare_GUI()
        self.Bias_GUI()
        self.Flat_GUI()
        self.Preproc_GUI()
        self.CosmCorr_GUI()
        self.BckgSub_GUI()
        self.Extract_GUI()
        self.WavCal_GUI()

    def Prepare_GUI(self):
        
        self.Prepare_load_button = Button(self.tab_Prepare, text="Load files", width=20)
        self.Prepare_load_button.grid(row=0, column=0, sticky=W)  
        self.Prepare_load_button.bind("<ButtonRelease-1>", self.load_file_prepare)    
 
        self.Prepare_run_button = Button(self.tab_Prepare, text="Run prepare", width=20)
        self.Prepare_run_button.grid(row=1, column=0, sticky=W)  
        self.Prepare_run_button.bind("<ButtonRelease-1>", self.Prepare)   
       
        self.frame_Prepare=Frame(self.tab_Prepare, width=1000, height=400)
        self.frame_Prepare.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_Prepare.grid_propagate(False)

        self.frame_Prepare.grid_rowconfigure(0, weight=1)
        self.frame_Prepare.grid_columnconfigure(0, weight=1)

        self.frame_Prepare.mytext3 = Text(self.frame_Prepare, state="disabled")
        self.frame_Prepare.mytext3.place(x=10, y=10, height=990, width=390)       
        
    def Bias_GUI(self):
        
        self.Bias_load_button = Button(self.tab_Bias, text="Load files", width=20)
        self.Bias_load_button.grid(row=0, column=0, sticky=W)  
        self.Bias_load_button.bind("<ButtonRelease-1>", self.load_file_bias)    
 
        self.Bias_run_button = Button(self.tab_Bias, text="Create Master Bias", width=20)
        self.Bias_run_button.grid(row=1, column=0, sticky=W)  
        self.Bias_run_button.bind("<ButtonRelease-1>", self.Bias)   
       
        self.frame_Bias=Frame(self.tab_Bias, width=1000, height=400)
        self.frame_Bias.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_Bias.grid_propagate(False)

        self.frame_Bias.grid_rowconfigure(0, weight=1)
        self.frame_Bias.grid_columnconfigure(0, weight=1)

        self.frame_Bias.mytext_Bias = Text(self.frame_Bias, state="disabled")
        self.frame_Bias.mytext_Bias.place(x=10, y=10, height=990, width=390)               

        Label(self.tab_Bias, text='Master bias name:').grid(row=3,column = 0)
        self.MasterBiasName = Entry(self.tab_Bias) 
        self.MasterBiasName.grid(row=3, column=1) 
        self.MasterBiasName.insert(END,'MasterBias.fits')
        
    def Flat_GUI(self):
        
        self.Flat_load_button = Button(self.tab_Flat, text="Load files", width=20)
        self.Flat_load_button.grid(row=0, column=0, sticky=W)  
        self.Flat_load_button.bind("<ButtonRelease-1>", self.load_file_flat)    
 
        self.Flat_run_button = Button(self.tab_Flat, text="Create Master Flat", width=20)
        self.Flat_run_button.grid(row=1, column=0, sticky=W)  
        self.Flat_run_button.bind("<ButtonRelease-1>", self.Flat)   
       
        self.frame_Flat=Frame(self.tab_Flat, width=1000, height=400)
        self.frame_Flat.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_Flat.grid_propagate(False)

        self.frame_Flat.grid_rowconfigure(0, weight=1)
        self.frame_Flat.grid_columnconfigure(0, weight=1)

        self.frame_Flat.mytext_Flat = Text(self.frame_Flat, state="disabled")
        self.frame_Flat.mytext_Flat.place(x=10, y=10, height=990, width=390)               

        Label(self.tab_Flat, text='Master flat name:').grid(row=4,column = 0)
        self.MasterFlatName = Entry(self.tab_Flat) 
        self.MasterFlatName.grid(row=4, column=1) 
        self.MasterFlatName.insert(END,'MasterFlat.fits')

        Label(self.tab_Flat, text='Master bias name:').grid(row=3,column = 0)
        self.MasterBiasName = Entry(self.tab_Flat) 
        self.MasterBiasName.grid(row=3, column=1) 
        self.MasterBiasName.insert(END,'MasterBias.fits')


    def Preproc_GUI(self):
        
        self.Preproc_load_button = Button(self.tab_Preproc, text="Load files", width=20)
        self.Preproc_load_button.grid(row=0, column=0, sticky=W)  
        self.Preproc_load_button.bind("<ButtonRelease-1>", self.load_file_preproc)    
 
        self.Preproc_run_button = Button(self.tab_Preproc, text="Process files", width=20)
        self.Preproc_run_button.grid(row=1, column=0, sticky=W)  
        self.Preproc_run_button.bind("<ButtonRelease-1>", self.Preproc)   
       
        self.frame_Preproc=Frame(self.tab_Preproc, width=1000, height=400)
        self.frame_Preproc.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_Preproc.grid_propagate(False)

        self.frame_Preproc.grid_rowconfigure(0, weight=1)
        self.frame_Preproc.grid_columnconfigure(0, weight=1)

        self.frame_Preproc.mytext_Preproc = Text(self.frame_Preproc, state="disabled")
        self.frame_Preproc.mytext_Preproc.place(x=10, y=10, height=990, width=390)               

        Label(self.tab_Preproc, text='Master flat name:').grid(row=4,column = 0)
        self.MasterFlatName = Entry(self.tab_Preproc) 
        self.MasterFlatName.grid(row=4, column=1) 
        self.MasterFlatName.insert(END,'MasterFlat.fits')

        Label(self.tab_Preproc, text='Master bias name:').grid(row=3,column = 0)
        self.MasterBiasName = Entry(self.tab_Preproc) 
        self.MasterBiasName.grid(row=3, column=1) 
        self.MasterBiasName.insert(END,'MasterBias.fits')


    def CosmCorr_GUI(self):
        
        self.CosmCorr_load_button = Button(self.tab_CosmCorr, text="Load files", width=20)
        self.CosmCorr_load_button.grid(row=0, column=0, sticky=W)  
        self.CosmCorr_load_button.bind("<ButtonRelease-1>", self.load_file_cosmcorr)    
 
        self.CosmCorr_run_button = Button(self.tab_CosmCorr, text="Correct cosmics", width=20)
        self.CosmCorr_run_button.grid(row=1, column=0, sticky=W)  
        self.CosmCorr_run_button.bind("<ButtonRelease-1>", self.Cosmcorr)   
       
        self.frame_CosmCorr=Frame(self.tab_CosmCorr, width=1000, height=400)
        self.frame_CosmCorr.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_CosmCorr.grid_propagate(False)

        self.frame_CosmCorr.grid_rowconfigure(0, weight=1)
        self.frame_CosmCorr.grid_columnconfigure(0, weight=1)

        self.frame_CosmCorr.mytext_Cosmcorr = Text(self.frame_CosmCorr, state="disabled")
        self.frame_CosmCorr.mytext_Cosmcorr.place(x=10, y=10, height=990, width=390)   


    def BckgSub_GUI(self):
        
        self.BckgSub_load_button = Button(self.tab_BckgSub, text="Load files", width=20)
        self.BckgSub_load_button.grid(row=0, column=0, sticky=W)  
        self.BckgSub_load_button.bind("<ButtonRelease-1>", self.load_file_bckgsub)    
 
        self.BckgSub_run_button = Button(self.tab_BckgSub, text="Remove background", width=20)
        self.BckgSub_run_button.grid(row=1, column=0, sticky=W)  
        self.BckgSub_run_button.bind("<ButtonRelease-1>", self.BckgSub)   
       
        self.frame_BckgSub=Frame(self.tab_BckgSub, width=1000, height=400)
        self.frame_BckgSub.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_BckgSub.grid_propagate(False)

        self.frame_BckgSub.grid_rowconfigure(0, weight=1)
        self.frame_BckgSub.grid_columnconfigure(0, weight=1)

        self.frame_BckgSub.mytext_BckgSub = Text(self.frame_BckgSub, state="disabled")
        self.frame_BckgSub.mytext_BckgSub.place(x=10, y=10, height=990, width=390)  
              


    def Extract_GUI(self):
        
        self.Extract_load_button = Button(self.tab_Extract, text="Load files", width=20)
        self.Extract_load_button.grid(row=0, column=0, sticky=W)  
        self.Extract_load_button.bind("<ButtonRelease-1>", self.load_file_Extract)    
 
        self.Extract_run_button = Button(self.tab_Extract, text="Extract spectra", width=20)
        self.Extract_run_button.grid(row=1, column=0, sticky=W)  
        self.Extract_run_button.bind("<ButtonRelease-1>", self.Extract)   
       
        self.frame_Extract=Frame(self.tab_Extract, width=1000, height=400)
        self.frame_Extract.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_Extract.grid_propagate(False)

        self.frame_Extract.grid_rowconfigure(0, weight=1)
        self.frame_Extract.grid_columnconfigure(0, weight=1)

        self.frame_Extract.mytext_Extract = Text(self.frame_Extract, state="disabled")
        self.frame_Extract.mytext_Extract.place(x=10, y=10, height=990, width=390)  

    def WavCal_GUI(self):
        
        self.WavCal_load_button = Button(self.tab_WavCal, text="Load arc files", width=20)
        self.WavCal_load_button.grid(row=0, column=0, sticky=W)  
        self.WavCal_load_button.bind("<ButtonRelease-1>", self.load_file_WavCal)    

        self.WavCal_Spec_load_button = Button(self.tab_WavCal, text="Load spec file", width=20)
        self.WavCal_Spec_load_button.grid(row=2, column=0, sticky=W)  
        self.WavCal_Spec_load_button.bind("<ButtonRelease-1>", self.load_file_WavCal_Spec)  
 
        self.WavCal_run_button = Button(self.tab_WavCal, text="Calibrate wavelength", width=20)
        self.WavCal_run_button.grid(row=3, column=0, sticky=W)  
        self.WavCal_run_button.bind("<ButtonRelease-1>", self.WavCal)   
       
        self.frame_WavCal=Frame(self.tab_WavCal, width=1000, height=400)
        self.frame_WavCal.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_WavCal.grid_propagate(False)

        self.frame_WavCal.grid_rowconfigure(0, weight=1)
        self.frame_Extract.grid_columnconfigure(0, weight=1)

        self.frame_WavCal.mytext_Extract = Text(self.frame_WavCal, state="disabled")
        self.frame_WavCal.mytext_Extract.place(x=10, y=10, height=990, width=390)  



############################################################################## 
# SP_WavCal
##############################################################################


    def load_file_WavCal(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_WavCal=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_WavCal.append(elem)
        
        self.read_all_fits(self.files_WavCal)
        print(self.images[0])
        
        self.canvas.delete("all")
#        self.tkimage.forget()
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image
        
        im = self.images[0]
        self.tkimage.paste(im)
            
        return None

    def load_file_WavCal_Spec(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_WavCal_Spec=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_WavCal_Spec.append(elem)
        
#        self.read_all_fits(self.files_WavCal)
#        print(self.images[0])
#        
#        self.canvas.delete("all")
##        self.tkimage.forget()
#        self.canvas.forget()
#        
#        
#        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
#        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
#                             self.tkimage.width())
#        self.canvas.pack()
#        self.image = self.canvas.create_image(0, 0, anchor='nw',
#                                              image=self.tkimage)  
#        
#              # select first image
#        
#        im = self.images[0]
#        self.tkimage.paste(im)
            
        return None



    def WavCal(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        os.system('python ' + Pipe_Path + '/SP_WavCal.py ' +  + '-a ' + " ".join(self.files_WavCal))
        self.now = datetime.datetime.now()
        module_logger_WavCal.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None





############################################################################## 
# SP_Extract
##############################################################################


    def load_file_Extract(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_Extract=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_Extract.append(elem)
        
        self.read_all_fits(self.files_Extract)
        print(self.images[0])
        
        self.canvas.delete("all")
#        self.tkimage.forget()
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image
        
        im = self.images[0]
        self.tkimage.paste(im)
            
        return None

    def Extract(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        os.system('python ' + Pipe_Path + '/SP_Extract.py ' + " ".join(self.files_Extract))
        self.now = datetime.datetime.now()
        module_logger_Extract.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None



        

############################################################################## 
# SP_Bckgsub
##############################################################################


    def read_all_fits(self, filenames):
        """ read in all image data, scale images """
        self.data=[]
        self.Tx = []
        self.Ty = []
        self.images = []
        for idx, filename in enumerate(filenames):
            ## read image data
            hdulist = fits.open(filename, ignore_missing_end=True)
            imgdat = hdulist[0].data
            self.data.append(imgdat)

            median = np.median(imgdat[int(imgdat.shape[1]*0.25):
                                         int(imgdat.shape[1]*0.75),
                                         int(imgdat.shape[0]*0.25):
                                         int(imgdat.shape[0]*0.75)])
            std    = np.std(imgdat[int(imgdat.shape[1]*0.25):
                                      int(imgdat.shape[1]*0.75),
                                      int(imgdat.shape[0]*0.25):
                                      int(imgdat.shape[0]*0.75)])

            imgdat = old_div(np.clip(imgdat, median-0.5*std,
                                median+0.5*std),(old_div(std,256)))
            imgdat = imgdat - np.min(imgdat)

            imgdat = interp.zoom(imgdat, self.zoom)

            self.images.append(Image.fromarray(imgdat))


    def load_file_bckgsub(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_bckgsub=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_bckgsub.append(elem)
        
        self.read_all_fits(self.files_bckgsub)
        print(self.images[0])
        
        self.canvas.delete("all")
#        self.tkimage.forget()
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image
        
        im = self.images[0]
        self.tkimage.paste(im)
            
        return None

    def BckgSub(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        os.system('python ' + Pipe_Path + '/SP_BckgSub.py ' + " ".join(self.files_bckgsub))
        self.now = datetime.datetime.now()
        module_logger_bckgsub.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None


############################################################################## 
# SP_Cosmcorr
##############################################################################


    def load_file_cosmcorr(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_cosmcorr=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_cosmcorr.append(elem)

        self.read_all_fits(self.files_cosmcorr)
        print(self.images[0])
        
        self.canvas.delete("all")
#        self.tkimage.forget()
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image
        
        im = self.images[0]
        self.tkimage.paste(im)



        return None

    def Cosmcorr(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        os.system('python ' + Pipe_Path + '/SP_CosmCorr.py ' + " ".join(self.files_cosmcorr))
        self.now = datetime.datetime.now()
        module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None

############################################################################## 
# SP_Preproc
##############################################################################

    def load_file_preproc(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_preproc=[]
        for elem in self.fname:
            module_logger_Preproc.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_preproc.append(elem)

        self.read_all_fits(self.files_preproc)
        print(self.images[0])
        
        self.canvas.delete("all")
#        self.tkimage.forget()
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image
        
        im = self.images[0]
        self.tkimage.paste(im)


        return None


    def Preproc(self, event):
        now = datetime.datetime.now()
        module_logger_Preproc.info(now)
        os.system('python ' + Pipe_Path + '/SP_Preproc.py ' + " ".join(self.files_preproc) + ' -b ' + self.MasterBiasName.get() + ' -f ' + self.MasterFlatName.get())
        self.now = datetime.datetime.now()
        module_logger_Flat.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'End of flat processing' )
        
        return None



############################################################################## 
# SP_Flat
##############################################################################

    def load_file_flat(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_flat=[]
        for elem in self.fname:
            module_logger_Flat.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_flat.append(elem)

        self.read_all_fits(self.files_flat)
        print(self.images[0])
        
        self.canvas.delete("all")
#        self.tkimage.forget()
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image
        
        im = self.images[0]
        self.tkimage.paste(im)



        return None


    def Flat(self, event):
        now = datetime.datetime.now()
        module_logger_Flat.info(now)
        os.system('python ' + Pipe_Path + '/SP_Flat.py ' + " ".join(self.files_flat) + ' -b ' + self.MasterBiasName.get() + ' -o ' + self.MasterFlatName.get())
        self.now = datetime.datetime.now()
        module_logger_Flat.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'End of flat processing' )
        
        return None



############################################################################## 
# SP_Bias
##############################################################################

    def load_file_bias(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_bias=[]
        for elem in self.fname:
            module_logger_Bias.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_bias.append(elem)


        self.read_all_fits(self.files_bias)
        print(self.images[0])
        
        self.canvas.delete("all")
#        self.tkimage.forget()
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image
        
        im = self.images[0]
        self.tkimage.paste(im)


        return None


    def Bias(self, event):
        now = datetime.datetime.now()
        module_logger_Bias.info(now)
        print(self.MasterBiasName.get())
        print('python ' + Pipe_Path + '/SP_Bias.py ' + " ".join(self.files_bias) + ' -o ' + self.MasterBiasName.get())
        os.system('python ' + Pipe_Path + '/SP_Bias.py ' + " ".join(self.files_bias) + ' -o ' + self.MasterBiasName.get())
        self.now = datetime.datetime.now()
        module_logger_Bias.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'End of bias processing' )
        
        return None
############################################################################## 
# SP_Prepare
##############################################################################


    def load_file_prepare(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_prepare=[]
        for elem in self.fname:
            module_logger2.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_prepare.append(elem.split('/')[-1])


        return None

    def Prepare(self,event):
        now = datetime.datetime.now()
        module_logger.info(now)
        os.system('python ' + Pipe_Path + '/SP_Prepare.py ' + " ".join(self.files_prepare))
        self.now = datetime.datetime.now()
        module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'File preparation done' )
        
        return None

############################################################################## 
# Front tab
##############################################################################    
    
    def load_file(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files.append(elem.split('/')[-1])
        return None

    def Plot_Spectrum(self):
        plt.figure()
        for elem in self.fname:
            Spec = np.loadtxt(elem).transpose()
            plt.plot(Spec[0],Spec[1],'.')
        plt.show()
        
        return None
    
    def Display_fits(self):
        self.now = datetime.datetime.now()
        module_logger.info(self.now)
        os.system('python ' + Pipe_Path + '/SP_DisplayFits.py ' + " ".join(self.fname))
        
        return None






class MyHandlerText(logging.StreamHandler):
    def __init__(self, textctrl):
        logging.StreamHandler.__init__(self) # initialize parent
        self.textctrl = textctrl

    def emit(self, record):
        msg = self.format(record)
        self.textctrl.config(state="normal")
        self.textctrl.insert("end", msg + "\n")
        self.flush()
        self.textctrl.config(state="disabled")
        

if __name__ == "__main__":
    
    # create Tk object instance
    app = simpleapp_tk(None)
    app.title('my application')

    # setup logging handlers using the Tk instance created above
    # the pattern below can be used in other threads...
    # ...to allow other thread to send msgs to the gui
    # in this example, we set up two handlers just for demonstration (you could add a fileHandler, etc)
#    stderrHandler = logging.StreamHandler()  # no arguments => stderr
#    module_logger.addHandler(stderrHandler)
#    guiHandler = MyHandlerText(app.frame.mytext2)
#    module_logger.addHandler(guiHandler)
#    module_logger.setLevel(logging.INFO)
#    module_logger.info("from main")   
#
#    stderrHandler2 = logging.StreamHandler()  # no arguments => stderr
#    module_logger2.addHandler(stderrHandler2)
#    guiHandler2 = MyHandlerText(app.frame_Prepare.mytext3)    
#    module_logger2.addHandler(guiHandler2)
#    module_logger2.setLevel(logging.INFO)
#    module_logger2.info("from main")          
#
#
#    stderrHandler_Bias = logging.StreamHandler()  # no arguments => stderr
#    module_logger_Bias.addHandler(stderrHandler_Bias)
#    guiHandler_Bias = MyHandlerText(app.frame_Bias.mytext_Bias)
#    module_logger_Bias.addHandler(guiHandler_Bias)
#    module_logger_Bias.setLevel(logging.INFO)
#    module_logger_Bias.info("from main")      
#
#    stderrHandler_Flat = logging.StreamHandler()  # no arguments => stderr
#    module_logger_Flat.addHandler(stderrHandler_Flat)
#    guiHandler_Flat = MyHandlerText(app.frame_Flat.mytext_Flat)
#    module_logger_Flat.addHandler(guiHandler_Flat)
#    module_logger_Flat.setLevel(logging.INFO)
#    module_logger_Flat.info("from main")      
#
#    stderrHandler_Preproc = logging.StreamHandler()  # no arguments => stderr
#    module_logger_Preproc.addHandler(stderrHandler_Preproc)
#    guiHandler_Preproc = MyHandlerText(app.frame_Preproc.mytext_Preproc)
#    module_logger_Preproc.addHandler(guiHandler_Preproc)
#    module_logger_Preproc.setLevel(logging.INFO)
#    module_logger_Preproc.info("from main")      

#    stderrHandler_CosmCorr = logging.StreamHandler()  # no arguments => stderr
#    module_logger_CosmCorr.addHandler(stderrHandler_CosmCorr)
#    guiHandler_CosmCorr = MyHandlerText(app.frame_CosmCorr.mytext_CosmCorr)
#    module_logger_CosmCorr.addHandler(guiHandler_CosmCorr)
#    module_logger_CosmCorr.setLevel(logging.INFO)
#    module_logger_CosmCorr.info("from main")      

#    stderrHandler_Bckgsub = logging.StreamHandler()  # no arguments => stderr
#    module_logger_Bckgsub.addHandler(stderrHandler_Bckgsub)
#    guiHandler_Bckgsub = MyHandlerText(app.frame_Bckgsub.mytext_Bckgsub)
#    module_logger_Bckgsub.addHandler(guiHandler_Bckgsub)
#    module_logger_Bckgsub.setLevel(logging.INFO)
#    module_logger_Bckgsub.info("from main")   


    # start Tk
    app.mainloop()
 