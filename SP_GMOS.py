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
from tkinter import messagebox

from ttkwidgets import CheckboxTreeview

from astropy.visualization import PercentileInterval, ZScaleInterval

from astropy.io import fits

import logging
import datetime

import simplejson


from collections import defaultdict


from past.utils import old_div

import operator

import SP_Bias
import SP_Preproc
import SP_Prepare
import SP_CosmCorr
import SP_BckgSub
import SP_Extract
import SP_Combine
import SP_Flat
import SP_Subtract
import SP_DetectSpectra

from PIL import Image
from PIL import ImageTk
from PIL import ImageDraw

import argparse, shlex
import numpy as np

from glob import glob


from scipy.ndimage import interpolation as interp

from scipy import misc

from scipy.optimize import curve_fit

from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

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
module_logger_Combine = logging.getLogger(__name__ + '.I')
module_logger_WavCal = logging.getLogger(__name__ + '.J')
module_logger_TellCorr = logging.getLogger(__name__ + '.K')



class simpleapp_tk(Tk):
    def __init__(self,parent):
        Tk.__init__(self,parent)
        
        # Closing handeler
        
        self.protocol("WM_DELETE_WINDOW", self.on_closing)    
        
        # Manage log text file
        self.now = datetime.datetime.now()
        self.logger = logging.getLogger("my logger")
        logging.basicConfig(filename='Log_' + self.now.strftime("%Y_%m_%d_%H_%M_%S") +'.txt', filemode='a', format='%(asctime)s \t %(message)s',level=logging.INFO,datefmt='%Y-%m-%d %H:%M:%S')

        # zoom factor of the fits files.
        self.zoom = 0.25
        
        self.parent = parent

        # To use grid based layout
        self.grid()
        
        # Definition of the frame containing the processing tabs
        self.Process_Height = 200
        self.Process_Width = 1000
        self.frame_process=Frame(self, width=self.Process_Width, height=self.Process_Height )
        self.frame_process.grid(column=0, row=0)

        # Definition of the frame containing the fits images or graphs
        self.Fits_Height = 600
        self.frame_fits=Frame(self, width=self.Process_Width, height=self.Fits_Height)
        self.frame_fits.grid(column=0, row=1)
        
        
#        self.frame_tb = Frame(self,width = 400,heigh = 50)        
#        self.frame_tb.grid(column=0, row=2)
#        
        
        # Definition of the tab containing the different SP functions
  
        self.tab_index= {}

      
        self.TabControl = ttk.Notebook(self.frame_process)
        self.tab1 = Frame(self.TabControl)
        self.TabControl.add(self.tab1,text='Analysis')
        self.tab_index.update({0:{'Name':'Analysis','Clicked':False}})
        
        # SP_Prepare tab
        self.tab_Prepare = Frame(self.TabControl)
        self.TabControl.add(self.tab_Prepare,text='Prepare')
        self.tab_index.update({1:{'Name':'Prepare','Clicked':False}})

        # SP_Bias tab
        self.tab_Bias = Frame(self.TabControl)
        self.TabControl.add(self.tab_Bias,text='Bias')
        self.tab_index.update({2:{'Name':'Bias','Clicked':False}})

        # SP_Flat tab
        self.tab_Flat = Frame(self.TabControl)
        self.TabControl.add(self.tab_Flat,text='Flat')
        self.tab_index.update({3:{'Name':'Flat','Clicked':False}})

        # SP_Preproc tab
        self.tab_Preproc = Frame(self.TabControl)
        self.TabControl.add(self.tab_Preproc,text='Pre processing')
        self.tab_index.update({4:{'Name':'Preproc','Clicked':False}})

        # SP_Cosmcorr tab
        self.tab_CosmCorr = Frame(self.TabControl)
        self.TabControl.add(self.tab_CosmCorr,text='Cosmic correction')
        self.tab_index.update({5:{'Name':'Cosmcorr','Clicked':False}})


        # SP_DetectSpectra tab
        self.tab_DetectSpectra = Frame(self.TabControl)
        self.TabControl.add(self.tab_DetectSpectra,text='Detect spectra')
        self.tab_index.update({6:{'Name':'DetectSpectra','Clicked':False}})


        # SP_Subtract tab
        self.tab_Subtract = Frame(self.TabControl)
        self.TabControl.add(self.tab_Subtract,text='Image subtraction')
        self.tab_index.update({6:{'Name':'Bckgsub','Clicked':False}})

        # SP_BckgSub tab
        self.tab_BckgSub = Frame(self.TabControl)
        self.TabControl.add(self.tab_BckgSub,text='Background subtraction')
        self.tab_index.update({7:{'Name':'Bckgsub','Clicked':False}})

        # SP_Extract tab
        self.tab_Extract = Frame(self.TabControl)
        self.TabControl.add(self.tab_Extract,text='Spectra extraction')        
        self.tab_index.update({8:{'Name':'Extract','Clicked':False}})

        # SP_Combine tab
        self.tab_Combine = Frame(self.TabControl)
        self.TabControl.add(self.tab_Combine,text='Combine spectra')
        self.TabControl.bind(self.tab_Combine,"<ButtonRelease-1>", self.load_file)    
        self.tab_index.update({9:{'Name':'Combine','Clicked':False}})

        # SP_WavCal tab
        self.tab_WavCal = Frame(self.TabControl)
        self.TabControl.add(self.tab_WavCal,text='Wavelength calibration')   
        self.tab_index.update({10:{'Name':'Wavcal','Clicked':False}})

        # SP_TellCorr tab
        self.tab_TellCorr = Frame(self.TabControl)
        self.TabControl.add(self.tab_TellCorr,text='Telluric correction')   
        self.tab_index.update({11:{'Name':'TellCorr','Clicked':False}})
        
        
        
        self.TabControl.pack(expand=1, fill="both") 
        

        ### load dummy image for canvas 
        
        Manos_img= misc.imread(Pipe_Path + '/manos_splash.png')
        self.tkimage = ImageTk.PhotoImage(Image.fromarray(Manos_img), palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
        self.images = []
        self.images.append(Image.fromarray(Manos_img))
        self.tkimage.paste(self.images[0])


        self.LF_button = Button(self.tab1, text="Load files", width=10)
        self.LF_button.grid(row=1, column=0, sticky=W)
        self.LF_button.bind("<ButtonRelease-1>", self.load_file)        
        
        self.Plot_button = Button(self.tab1, text="Plot Spectrum", command=self.Plot_Spectrum, width=20)
        self.Plot_button.grid(row=2, column=0, sticky=W)   
        
        self.Fits_button = Button(self.tab1, text="Display fits", command=self.Display_fits, width=20)
        self.Fits_button.grid(row=3, column=0, sticky=W)   

        self.frame=Frame(self.tab1, width=1000, height=self.Process_Height)
        self.frame.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame.grid_propagate(False)

        self.frame.grid_rowconfigure(0, weight=1)
        self.frame.grid_columnconfigure(0, weight=1)

        self.frame.mytext2 = Text(self.frame, state="disabled")
        self.frame.mytext2.place(x=10, y=10, height=990, width=self.Process_Height-10)

        scrollb = Scrollbar(self.frame, command=self.frame.mytext2.yview)
        scrollb.grid(row=0, column=1, sticky='nsew')
        self.frame.mytext2['yscrollcommand'] = scrollb.set
        
        
        # Initialization of some needed variables
        
        self.files_Combine = []
        
        
#
#        # Running the pipeline:
#          
        self.Prepare_GUI()
        self.Bias_GUI()
        self.Flat_GUI()
        self.Preproc_GUI()
        self.CosmCorr_GUI()
        self.Subtract_GUI()
        self.DetectSpectra_GUI()
        self.BckgSub_GUI()
        self.Extract_GUI()
        self.Combine_GUI()
        self.WavCal_GUI()
        self.TellCorr_GUI()

        self.fig = plt.figure()
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.lastindex = []

        self.TabControl.bind("<<NotebookTabChanged>>", self.handle_tab_changed)

    def handle_tab_changed(self,event):

        selection = event.widget.select()
        tab = event.widget.index("current")
        
        if self.lastindex != self.TabControl.index(self.tab_Combine):
            print(self.lastindex)
            self.canvas.delete("all")
            self.canvas.forget()  
            
            if tab == self.TabControl.index(self.tab_Combine):
                self.canvas = FigureCanvasTkAgg(self.fig, self.frame_fits)
                self.canvas.show()
                self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1) 
            else: 
#                self.canvas.get_tk_widget().destroy()         
                self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
                self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
                self.canvas.pack()
                self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        else: 
            self.canvas.get_tk_widget().destroy()         
            if tab == self.TabControl.index(self.tab_Combine):
                self.canvas = FigureCanvasTkAgg(self.fig, self.frame_fits)
                self.canvas.show()
                self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1) 
            else: 
#                self.canvas.get_tk_widget().destroy()         
                self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
                self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
                self.canvas.pack()
                self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
                
                
        # Action to be performed when BIAS tab is clicked for the first time 
        if self.tab_index[int(tab)]['Name'] == 'Bias' and self.tab_index[int(tab)]['Clicked'] == False:
            self.files_bias_name = glob('*Bias*fits')
            for elem in self.files_bias_name:
                self.Pop_bias_tree(elem)
                
            Children = self.tree_bias.get_children()
            for elem in Children:
                self.tree_bias.change_state(elem, "checked")

            self.files_bias = []
            for elem in self.files_bias_name:
                self.files_bias.append(self.cdw + '/' + elem)                
        
        # Action to be performed when FLAT tab is clicked for the first time 
        if self.tab_index[int(tab)]['Name'] == 'Flat' and self.tab_index[int(tab)]['Clicked'] == False:
            self.files_flat_name = glob('*GCALflat*fits')

            for elem in self.files_flat_name:
                self.Pop_flat_tree(elem)
               
            Children = self.tree_flat.get_children()
            for elem in Children:
                self.tree_flat.change_state(elem, "checked")

            self.files_flat = []
            for elem in self.files_flat_name:
                self.files_flat.append(self.cdw + '/' + elem)


        # Action to be performed when Preproc tab is clicked for the first time 
        if self.tab_index[int(tab)]['Name'] == 'Preproc' and self.tab_index[int(tab)]['Clicked'] == False:
            print('First time clicked')
            self.files_pre_name = glob('*fits')
            Biases = []
            Flats = []
            ToProc = []
            Obj = []
            Object = {}
            Wavel = []
            
            Object = defaultdict(list)

            for elem in self.files_pre_name:
                header = fits.getheader(elem)
                
                try: 
                    if header['PROCTYPE'] == 'MASTER FLAT':
                        Flats.append(elem)
                    if header['PROCTYPE'] == 'MASTER BIAS':
                        Biases.append(elem)
                    if header['OBSTYPE'] == 'OBJECT' and header['PROCTYPE'] == 'Prepared':
                        Wavelength = header['WAVELENG']
                        ToProc.append(elem)
                        Object[header['Object']+'_'+str(int(Wavelength)/10)].append(elem)
                                     
                except:
                    pass


            for elem in Biases:
                self.Pop_Prep_bias_tree(elem)

            Children = self.tree_preproc_MBIAS.get_children()
            self.tree_preproc_MBIAS.change_state(Children[0], "checked")

            
            for elem in Flats:
                self.Pop_Prep_flat_tree(elem)

            Children = self.tree_preproc_MFLAT.get_children()
            self.tree_preproc_MFLAT.change_state(Children[0], "checked")

            for elem in Object:
                self.Create_Folder_toproc_tree(Object[elem],elem)


        if self.tab_index[int(tab)]['Name'] == 'Cosmcorr' and self.tab_index[int(tab)]['Clicked'] == False:
            print('First time clicked')
            self.files_cosm_name = glob('*fits')
#            Biases = []
#            Flats = []
#            ToProc = []
#            Obj = []
#            Object = {}
#            Wavel = []
#            
#            Object = defaultdict(list)

            for elem in self.files_cosm_name:
                header = fits.getheader(elem)
                try: 
                    if 'True' in header['SP_PREPR']:
                        print(elem)
                        self.Pop_Cosmcorr_tree(elem)
                except:
                    pass

            Children = self.tree_CosmCorr.get_children()
            for elem in Children:
                self.tree_CosmCorr.change_state(elem, "checked")

            
        self.lastindex = event.widget.index("current")       


    def Prepare_GUI(self):
        
#        self.Prepare_load_button = Button(self.tab_Prepare, text="Load files", width=20)
#        self.Prepare_load_button.grid(row=0, column=0, sticky=W)  
#        self.Prepare_load_button.bind("<ButtonRelease-1>", self.load_file_prepare)    
# 
#        self.Prepare_run_button = Button(self.tab_Prepare, text="Run prepare", width=20)
#        self.Prepare_run_button.grid(row=1, column=0, sticky=W)  
#        self.Prepare_run_button.bind("<ButtonRelease-1>", self.Prepare)   
#       
#        self.frame_Prepare=Frame(self.tab_Prepare, width=1000, height=self.Process_Height)
#        self.frame_Prepare.grid(column=2, row=0,rowspan=25,columnspan=10)
#
#        self.frame_Prepare.grid_propagate(False)
#
#        self.frame_Prepare.grid_rowconfigure(0, weight=1)
#        self.frame_Prepare.grid_columnconfigure(0, weight=1)
#
#        self.frame_Prepare.mytext3 = Text(self.frame_Prepare, state="disabled")
#        self.frame_Prepare.mytext3.place(x=10, y=10, height=990, width=390)       


        self.frame_Prepare_Button=Frame(self.tab_Prepare, width=100, height=400)
        self.frame_Prepare_Button.grid(column=0, row=0)

        self.Prepare_load_button = Button(self.frame_Prepare_Button, text="Load files", width=20)
        self.Prepare_load_button.grid(row=0, column=0, sticky=W)  
        self.Prepare_load_button.bind("<ButtonRelease-1>", self.load_file_prepare)    
 
        self.Prepare_run_button = Button(self.frame_Prepare_Button, text="Run prepare", width=20)
        self.Prepare_run_button.grid(row=1, column=0, sticky=W)  
        self.Prepare_run_button.bind("<ButtonRelease-1>", self.Prepare)   
       
        self.frame_Prepare=Frame(self.tab_Prepare, width=900, height=1000)
        self.frame_Prepare.grid(column=1, row=0)




        self.tree_prepare = CheckboxTreeview(self.frame_Prepare)
        self.tree_prepare["columns"]=("Type","ExpTime","Date","Grating","GratAng","Target","RotAng")
        self.tree_prepare.column("#0", width=250)
        self.tree_prepare.column("Type", width=100)
        self.tree_prepare.column("ExpTime", width=100)
        self.tree_prepare.column("Date", width=100)
        self.tree_prepare.column("Grating", width=100)
        self.tree_prepare.column("GratAng", width=100)
        self.tree_prepare.column("Target", width=100)
        self.tree_prepare.column("RotAng", width=100)
        
        self.tree_prepare.heading("#0",text="File name",anchor=W)
        self.tree_prepare.heading("Type",text="Type",anchor=CENTER)
        self.tree_prepare.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_prepare.heading("Date", text="Date",anchor=W)
        self.tree_prepare.heading("Grating", text="Grating",anchor=W)
        self.tree_prepare.heading("GratAng", text="Grat. ang.",anchor=W)        
        self.tree_prepare.heading("Target", text="Target",anchor=W)   
        self.tree_prepare.heading("RotAng", text="Rot. ang.",anchor=W)   
        
        self.tree_prepare.grid(column=1,row = 0)






        
    def Bias_GUI(self):
        

        self.Frame_Bias_Button = Frame(self.tab_Bias, width=200, height=400)
        self.Frame_Bias_Button.grid(row = 0,column =0)
        
        self.Frame_Bias_Info = Frame(self.tab_Bias, width=800, height=400)
        self.Frame_Bias_Info.grid(row = 0,column =1)
        
        Row_0_Bias = 0
        
        self.Bias_load_button = Button(self.Frame_Bias_Button, text="Load files", width=20)
        self.Bias_load_button.grid(row=Row_0_Bias, column=0, sticky=W)  
        self.Bias_load_button.bind("<ButtonRelease-1>", self.load_file_bias)   
        
        Row_0_Bias +=1

    
        self.menu_BiasCombine = StringVar(self.Frame_Bias_Button)
        self.menu_BiasCombine.set("Median") # default value   
        
        Label(self.Frame_Bias_Button, text='Method: ').grid(row=Row_0_Bias,column = 0)
        self.Bias_Method = OptionMenu(self.Frame_Bias_Button, self.menu_BiasCombine, "Median", "Mean")
        self.Bias_Method.grid(column=Row_0_Bias, row=1, sticky=W)
        
        
        Row_0_Bias +=1 
    
        self.Bias_run_button = Button(self.Frame_Bias_Button, text="Create Master Bias", width=20)
        self.Bias_run_button.grid(row=Row_0_Bias, column=0, sticky=W)  
        self.Bias_run_button.bind("<ButtonRelease-1>", self.Bias)   
       
        Row_0_Bias +=1

        Label(self.Frame_Bias_Button, text='Master bias name:').grid(row=Row_0_Bias,column = 0)
        self.MasterBiasName_Bias = Entry(self.Frame_Bias_Button) 
        self.MasterBiasName_Bias.grid(row=Row_0_Bias, column=1) 
        self.MasterBiasName_Bias.insert(END,'MasterBias.fits')

       
        
        self.tree_bias = CheckboxTreeview(self.Frame_Bias_Info)
        self.tree_bias["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_bias.column("#0", width=400)
        self.tree_bias.column("Type", width=50)
        self.tree_bias.column("ExpTime", width=50)
        self.tree_bias.column("Mean", width=70)
        self.tree_bias.column("Median", width=70)
        self.tree_bias.column("std", width=70)
        self.tree_bias.column("Date", width=200)
        
        self.tree_bias.heading("#0",text="File name",anchor=W)
        self.tree_bias.heading("Type",text="Type",anchor=CENTER)
        self.tree_bias.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_bias.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_bias.heading("Median", text="Median",anchor=CENTER)
        self.tree_bias.heading("std", text="std",anchor=CENTER)
        self.tree_bias.heading("Date", text="Date",anchor=W)
        
        self.tree_bias.grid(column=1,row = 0)          
        
        
        
        
        
        
    def Flat_GUI(self):
        
#        self.Flat_load_button = Button(self.tab_Flat, text="Load files", width=20)
#        self.Flat_load_button.grid(row=0, column=0, sticky=W)  
#        self.Flat_load_button.bind("<ButtonRelease-1>", self.load_file_flat)    
# 
#        self.Flat_run_button = Button(self.tab_Flat, text="Create Master Flat", width=20)
#        self.Flat_run_button.grid(row=1, column=0, sticky=W)  
#        self.Flat_run_button.bind("<ButtonRelease-1>", self.Flat)   
#       
#        self.frame_Flat=Frame(self.tab_Flat, width=1000, height=self.Process_Height)
#        self.frame_Flat.grid(column=2, row=0,rowspan=25,columnspan=10)
#
#        self.frame_Flat.grid_propagate(False)
#
#        self.frame_Flat.grid_rowconfigure(0, weight=1)
#        self.frame_Flat.grid_columnconfigure(0, weight=1)
#
#        self.frame_Flat.mytext_Flat = Text(self.frame_Flat, state="disabled")
#        self.frame_Flat.mytext_Flat.place(x=10, y=10, height=990, width=390)               
#
#        Label(self.tab_Flat, text='Master flat name:').grid(row=4,column = 0)
#        self.MasterFlatName_Flat = Entry(self.tab_Flat) 
#        self.MasterFlatName_Flat.grid(row=4, column=1) 
#        self.MasterFlatName_Flat.insert(END,'MasterFlat.fits')
#
#        Label(self.tab_Flat, text='Master bias name:').grid(row=3,column = 0)
#        self.MasterBiasName_Flat = Entry(self.tab_Flat) 
#        self.MasterBiasName_Flat.grid(row=3, column=1) 
#        self.MasterBiasName_Flat.insert(END,'MasterBias.fits')
        
        self.Frame_Flat_Button = Frame(self.tab_Flat, width=200, height=400)
        self.Frame_Flat_Button.grid(row = 0,column =0)
        
        self.Frame_Flat_Info = Frame(self.tab_Flat, width=800, height=400)
        self.Frame_Flat_Info.grid(row = 0,column =1)        
        
        
        self.Flat_load_button = Button(self.Frame_Flat_Button, text="Load files", width=20)
        self.Flat_load_button.grid(row=0, column=0, sticky=W)  
        self.Flat_load_button.bind("<ButtonRelease-1>", self.load_file_flat)    
 
        self.Flat_run_button = Button(self.Frame_Flat_Button, text="Create Master Flat", width=20)
        self.Flat_run_button.grid(row=1, column=0, sticky=W)  
        self.Flat_run_button.bind("<ButtonRelease-1>", self.Flat)   
                     

        Label(self.Frame_Flat_Button, text='Master flat name:').grid(row=4,column = 0)
        self.MasterFlatName_Flat = Entry(self.Frame_Flat_Button) 
        self.MasterFlatName_Flat.grid(row=4, column=1) 
        self.MasterFlatName_Flat.insert(END,'MasterFlat.fits')

        Label(self.Frame_Flat_Button, text='Master bias name:').grid(row=3,column = 0)
        self.MasterBiasName = Entry(self.Frame_Flat_Button) 
        self.MasterBiasName.grid(row=3, column=1) 
        self.MasterBiasName.insert(END,'MasterBias.fits')        
        
        
        self.tree_flat = CheckboxTreeview(self.Frame_Flat_Info)
        self.tree_flat["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_flat.column("#0", width=400)
        self.tree_flat.column("Type", width=50)
        self.tree_flat.column("ExpTime", width=50)
        self.tree_flat.column("Mean", width=70)
        self.tree_flat.column("Median", width=70)
        self.tree_flat.column("std", width=70)
        self.tree_flat.column("Date", width=200)
        
        self.tree_flat.heading("#0",text="File name",anchor=W)
        self.tree_flat.heading("Type",text="Type",anchor=CENTER)
        self.tree_flat.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_flat.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_flat.heading("Median", text="Median",anchor=CENTER)
        self.tree_flat.heading("std", text="std",anchor=CENTER)
        self.tree_flat.heading("Date", text="Date",anchor=W)

        self.tree_flat.grid(column=1,row = 0)   


    def Preproc_GUI(self):
        
#        self.Preproc_load_button = Button(self.tab_Preproc, text="Load files", width=20)
#        self.Preproc_load_button.grid(row=0, column=0, sticky=W)  
#        self.Preproc_load_button.bind("<ButtonRelease-1>", self.load_file_preproc)    
# 
#        self.Preproc_run_button = Button(self.tab_Preproc, text="Process files", width=20)
#        self.Preproc_run_button.grid(row=1, column=0, sticky=W)  
#        self.Preproc_run_button.bind("<ButtonRelease-1>", self.Preproc)   
#       
#        self.frame_Preproc=Frame(self.tab_Preproc, width=1000, height=self.Process_Height)
#        self.frame_Preproc.grid(column=2, row=0,rowspan=25,columnspan=10)
#
#        self.frame_Preproc.grid_propagate(False)
#
#        self.frame_Preproc.grid_rowconfigure(0, weight=1)
#        self.frame_Preproc.grid_columnconfigure(0, weight=1)
#
#        self.frame_Preproc.mytext_Preproc = Text(self.frame_Preproc, state="disabled")
#        self.frame_Preproc.mytext_Preproc.place(x=10, y=10, height=990, width=390)               
#
#        Label(self.tab_Preproc, text='Master flat name:').grid(row=4,column = 0)
#        self.MasterFlatName = Entry(self.tab_Preproc) 
#        self.MasterFlatName.grid(row=4, column=1) 
#        self.MasterFlatName.insert(END,'MasterFlat.fits')
#
#        Label(self.tab_Preproc, text='Master bias name:').grid(row=3,column = 0)
#        self.MasterBiasName = Entry(self.tab_Preproc) 
#        self.MasterBiasName.grid(row=3, column=1) 
#        self.MasterBiasName.insert(END,'MasterBias.fits')


        self.Frame_Preproc_Button = Frame(self.tab_Preproc, width=200, height=400)
        self.Frame_Preproc_Button.grid(row = 0,column =0)
        
        self.Frame_Preproc_Info = Frame(self.tab_Preproc, width=800, height=400)
        self.Frame_Preproc_Info.grid(row = 0,column =1)  

        
        self.Preproc_load_button = Button(self.Frame_Preproc_Button, text="Load files", width=20)
        self.Preproc_load_button.grid(row=0, column=0, sticky=W)  
        self.Preproc_load_button.bind("<ButtonRelease-1>", self.load_file_preproc)    
 
        self.Preproc_run_button = Button(self.Frame_Preproc_Button, text="Process files", width=20)
        self.Preproc_run_button.grid(row=1, column=0, sticky=W)  
        self.Preproc_run_button.bind("<ButtonRelease-1>", self.Preproc)   
       
        Label(self.Frame_Preproc_Button, text='Master flat name:').grid(row=4,column = 0)
        self.MasterFlatName = Entry(self.Frame_Preproc_Button) 
        self.MasterFlatName.grid(row=4, column=1) 
        self.MasterFlatName.insert(END,'MasterFlat.fits')

        Label(self.Frame_Preproc_Button, text='Master bias name:').grid(row=3,column = 0)
        self.MasterBiasName = Entry(self.Frame_Preproc_Button) 
        self.MasterBiasName.grid(row=3, column=1) 
        self.MasterBiasName.insert(END,'MasterBias.fits')



        self.tree_preproc_MBIAS = CheckboxTreeview(self.Frame_Preproc_Info, height =3)
        self.tree_preproc_MBIAS["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_preproc_MBIAS.column("#0", width=400)
        self.tree_preproc_MBIAS.column("Type", width=50)
        self.tree_preproc_MBIAS.column("ExpTime", width=50)
        self.tree_preproc_MBIAS.column("Mean", width=70)
        self.tree_preproc_MBIAS.column("Median", width=70)
        self.tree_preproc_MBIAS.column("std", width=70)
        self.tree_preproc_MBIAS.column("Date", width=200)
        
        self.tree_preproc_MBIAS.heading("#0",text="File name",anchor=W)
        self.tree_preproc_MBIAS.heading("Type",text="Type",anchor=CENTER)
        self.tree_preproc_MBIAS.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_preproc_MBIAS.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_preproc_MBIAS.heading("Median", text="Median",anchor=CENTER)
        self.tree_preproc_MBIAS.heading("std", text="std",anchor=CENTER)
        self.tree_preproc_MBIAS.heading("Date", text="Date",anchor=W)
        
        self.tree_preproc_MBIAS.grid(column=1,row = 0)    

        self.tree_preproc_MFLAT = CheckboxTreeview(self.Frame_Preproc_Info, height =3)
        self.tree_preproc_MFLAT["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_preproc_MFLAT.column("#0", width=400)
        self.tree_preproc_MFLAT.column("Type", width=50)
        self.tree_preproc_MFLAT.column("ExpTime", width=50)
        self.tree_preproc_MFLAT.column("Mean", width=70)
        self.tree_preproc_MFLAT.column("Median", width=70)
        self.tree_preproc_MFLAT.column("std", width=70)
        self.tree_preproc_MFLAT.column("Date", width=200)
        
        self.tree_preproc_MFLAT.heading("#0",text="File name",anchor=W)
        self.tree_preproc_MFLAT.heading("Type",text="Type",anchor=CENTER)
        self.tree_preproc_MFLAT.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_preproc_MFLAT.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_preproc_MFLAT.heading("Median", text="Median",anchor=CENTER)
        self.tree_preproc_MFLAT.heading("std", text="std",anchor=CENTER)
        self.tree_preproc_MFLAT.heading("Date", text="Date",anchor=W)
        
        self.tree_preproc_MFLAT.grid(column=1,row = 1)    

        self.tree_preproc_TOPROC = CheckboxTreeview(self.Frame_Preproc_Info, height = 10)
        self.tree_preproc_TOPROC["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_preproc_TOPROC.column("#0", width=400)
        self.tree_preproc_TOPROC.column("Type", width=50)
        self.tree_preproc_TOPROC.column("ExpTime", width=50)
        self.tree_preproc_TOPROC.column("Mean", width=70)
        self.tree_preproc_TOPROC.column("Median", width=70)
        self.tree_preproc_TOPROC.column("std", width=70)
        self.tree_preproc_TOPROC.column("Date", width=200)
        
        self.tree_preproc_TOPROC.heading("#0",text="File name",anchor=W)
        self.tree_preproc_TOPROC.heading("Type",text="Type",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("Median", text="Median",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("std", text="std",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("Date", text="Date",anchor=W)
        
        self.tree_preproc_TOPROC.grid(column=1,row = 2)    



    def CosmCorr_GUI(self):
        
#        self.CosmCorr_load_button = Button(self.tab_CosmCorr, text="Load files", width=20)
#        self.CosmCorr_load_button.grid(row=0, column=0, sticky=W)  
#        self.CosmCorr_load_button.bind("<ButtonRelease-1>", self.load_file_cosmcorr)    
# 
#        self.CosmCorr_run_button = Button(self.tab_CosmCorr, text="Correct cosmics", width=20)
#        self.CosmCorr_run_button.grid(row=1, column=0, sticky=W)  
#        self.CosmCorr_run_button.bind("<ButtonRelease-1>", self.Cosmcorr)   
#       
#        self.frame_CosmCorr=Frame(self.tab_CosmCorr, width=1000, height=self.Process_Height)
#        self.frame_CosmCorr.grid(column=2, row=0,rowspan=25,columnspan=10)
#
#        self.frame_CosmCorr.grid_propagate(False)
#
#        self.frame_CosmCorr.grid_rowconfigure(0, weight=1)
#        self.frame_CosmCorr.grid_columnconfigure(0, weight=1)
#
#        self.frame_CosmCorr.mytext_Cosmcorr = Text(self.frame_CosmCorr, state="disabled")
#        self.frame_CosmCorr.mytext_Cosmcorr.place(x=10, y=10, height=990, width=390)   

        self.Frame_CosmCorr_Button = Frame(self.tab_CosmCorr, width=200, height=400)
        self.Frame_CosmCorr_Button.grid(row = 0,column =0)
        
        self.Frame_CosmCorr_Info = Frame(self.tab_CosmCorr, width=800, height=400)
        self.Frame_CosmCorr_Info.grid(row = 0,column =1)          
        
        
        self.CosmCorr_load_button = Button(self.Frame_CosmCorr_Button, text="Load files", width=20)
        self.CosmCorr_load_button.grid(row=0, column=0, sticky=W)  
        self.CosmCorr_load_button.bind("<ButtonRelease-1>", self.load_file_cosmcorr)    
 
        self.CosmCorr_run_button = Button(self.Frame_CosmCorr_Button, text="Correct cosmics", width=20)
        self.CosmCorr_run_button.grid(row=1, column=0, sticky=W)  
        self.CosmCorr_run_button.bind("<ButtonRelease-1>", self.Cosmcorr)   

        self.tree_CosmCorr = CheckboxTreeview(self.Frame_CosmCorr_Info)
        self.tree_CosmCorr["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_CosmCorr.column("#0", width=400)
        self.tree_CosmCorr.column("Type", width=50)
        self.tree_CosmCorr.column("ExpTime", width=50)
        self.tree_CosmCorr.column("Mean", width=70)
        self.tree_CosmCorr.column("Median", width=70)
        self.tree_CosmCorr.column("std", width=70)
        self.tree_CosmCorr.column("Date", width=200)
        
        self.tree_CosmCorr.heading("#0",text="File name",anchor=W)
        self.tree_CosmCorr.heading("Type",text="Type",anchor=CENTER)
        self.tree_CosmCorr.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_CosmCorr.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_CosmCorr.heading("Median", text="Median",anchor=CENTER)
        self.tree_CosmCorr.heading("std", text="std",anchor=CENTER)
        self.tree_CosmCorr.heading("Date", text="Date",anchor=W)
        
        self.tree_CosmCorr.grid(column=1,row = 0)            
        
        

    def DetectSpectra_GUI(self):
        
        self.DetectSpectra_load_button = Button(self.tab_DetectSpectra, text="Load files", width=20)
        self.DetectSpectra_load_button.grid(row=0, column=0, sticky=W)  
        self.DetectSpectra_load_button.bind("<ButtonRelease-1>", self.load_file_detectspectra)   

        self.DetectSpectra_run_button = Button(self.tab_DetectSpectra, text="Detect spectra", width=20)
        self.DetectSpectra_run_button.grid(row=1, column=0, sticky=W)  
        self.DetectSpectra_run_button.bind("<ButtonRelease-1>", self.Detect_Spectra)        
        
        Label(self.tab_DetectSpectra, text='Detection sensitivity:').grid(row=2,column = 0)
        self.Det_Sens = Entry(self.tab_DetectSpectra) 
        self.Det_Sens.grid(row=2, column=1) 
        self.Det_Sens.insert(END,'3')   
 
        Label(self.tab_DetectSpectra, text='Output file name:').grid(row=3,column = 0)
        self.Output_Name_Det = Entry(self.tab_DetectSpectra) 
        self.Output_Name_Det.grid(row=3, column=1) 
        self.Output_Name_Det.insert(END,'DetSpec1.fits')        

    def Subtract_GUI(self):
        
        # Add a button to load the fits files to subtract from each others
        self.Subtract_load_button = Button(self.tab_Subtract, text="Load files", width=20)
        self.Subtract_load_button.grid(row=0, column=0, sticky=W)  
        self.Subtract_load_button.bind("<ButtonRelease-1>", self.load_file_subtract)   

        self.Subtract_run_button = Button(self.tab_Subtract, text="Subtract images", width=20)
        self.Subtract_run_button.grid(row=1, column=0, sticky=W)  
        self.Subtract_run_button.bind("<ButtonRelease-1>", self.Subtract_images)          
 
        Label(self.tab_Subtract, text='Output file name:').grid(row=3,column = 0)
        self.Output_Name_Sub = Entry(self.tab_Subtract) 
        self.Output_Name_Sub.grid(row=3, column=1) 
        self.Output_Name_Sub.insert(END,'Sub1.fits')       

    def BckgSub_GUI(self):
        
        self.BckgSub_load_button = Button(self.tab_BckgSub, text="Load files", width=20)
        self.BckgSub_load_button.grid(row=0, column=0, sticky=W)  
        self.BckgSub_load_button.bind("<ButtonRelease-1>", self.load_file_bckgsub)    

        self.BckgSub_LocSpec_button = Button(self.tab_BckgSub, text="Loc spec file", width=20)
        self.BckgSub_LocSpec_button.grid(row=1, column=0, sticky=W)  
        self.BckgSub_LocSpec_button.bind("<ButtonRelease-1>", self.load_LocSpec)  

 
        self.BckgSub_run_button = Button(self.tab_BckgSub, text="Remove background", width=20)
        self.BckgSub_run_button.grid(row=2, column=0, sticky=W)  
        self.BckgSub_run_button.bind("<ButtonRelease-1>", self.BckgSub)   
       
        self.frame_BckgSub=Frame(self.tab_BckgSub, width=1000, height=self.Process_Height)
        self.frame_BckgSub.grid(column=3, row=0,rowspan=25,columnspan=10)

        self.frame_BckgSub.grid_propagate(False)

        self.frame_BckgSub.grid_rowconfigure(0, weight=1)
        self.frame_BckgSub.grid_columnconfigure(0, weight=1)

        self.frame_BckgSub.mytext_BckgSub = Text(self.frame_BckgSub, state="disabled")
        self.frame_BckgSub.mytext_BckgSub.place(x=10, y=10, height=990, width=390)  
              
        
        self.Live_BckgSub = IntVar()
        self.Check_Live_BckgSub = Checkbutton(self.tab_BckgSub, text="Live processing", variable=self.Live_BckgSub).grid(column = 0,row=3)


    def Extract_GUI(self):
        
        # using index for rows to allow to easily add one row if needed
        # the suffix correspond to the respetive column
        # the row index is always incremented once it was assigned 
        row_index_0 = 0
        row_index_1 = 0
        row_index_2 = 0
        
        self.Extract_load_button = Button(self.tab_Extract, text="Load files", width=20)
        self.Extract_load_button.grid(row=row_index_0, column=0, sticky=W)  
        self.Extract_load_button.bind("<ButtonRelease-1>", self.load_file_Extract)    
        row_index_0 += 1 # 1

        self.Extract_LocSpec_button = Button(self.tab_Extract, text="Loc spec file", width=20)
        self.Extract_LocSpec_button.grid(row=row_index_0, column=0, sticky=W)  
        self.Extract_LocSpec_button.bind("<ButtonRelease-1>", self.load_LocSpec)   
        row_index_0 += 1 # 2    

        self.Extract_Bckg_button = Button(self.tab_Extract, text="Background file", width=20)
        self.Extract_Bckg_button.grid(row=row_index_0, column=0, sticky=W)  
        self.Extract_Bckg_button.bind("<ButtonRelease-1>", self.load_Bckg)   
        row_index_0 += 1 # 3    
    
        self.Extract_run_button = Button(self.tab_Extract, text="Extract spectra", width=20)
        self.Extract_run_button.grid(row=row_index_0, column=0, sticky=W)  
        self.Extract_run_button.bind("<ButtonRelease-1>", self.Extract)   
        row_index_0 += 1 # 4       
        
        self.frame_Extract=Frame(self.tab_Extract, width=1000, height=self.Process_Height)
        self.frame_Extract.grid(column=2, row=row_index_2,rowspan=25,columnspan=10)
        row_index_2 += 1 #1

        self.frame_Extract.grid_propagate(False)

        self.frame_Extract.grid_rowconfigure(0, weight=1)
        self.frame_Extract.grid_columnconfigure(0, weight=1)

        self.frame_Extract.mytext_Extract = Text(self.frame_Extract, state="disabled")
        self.frame_Extract.mytext_Extract.place(x=10, y=10, height=990, width=390)  

        Label(self.tab_Extract, text='# of pixel to extract:').grid(row=row_index_0,column = 0)
        row_index_0 += 4 #1     
        self.FWHM = Entry(self.tab_Extract) 
        self.FWHM.grid(row=row_index_0 - 1, column=1) #  row_index_0 - 1 To have same index as the Label 
        self.FWHM.insert(END,'6')

        self.Live_Extract = IntVar()
        self.Check_Live_Extract = Checkbutton(self.tab_Extract, text="Live processing", variable=self.Live_Extract).grid(column = 0,row=4)




    def Combine_GUI(self):
        
        # Add a button to load the spectra to combine
        self.Combine_load_button = Button(self.tab_Combine, text="Load files", width=20)
        self.Combine_load_button.grid(row=0, column=0, sticky=W)  
        self.Combine_load_button.bind("<ButtonRelease-1>", self.load_file_Combine)    
 
        # Add a button to combine the selected spectra
        self.Combine_run_button = Button(self.tab_Combine, text="Combine spectra", width=20)
        self.Combine_run_button.grid(row=1, column=0, sticky=W)  
        self.Combine_run_button.bind("<ButtonRelease-1>", self.Combine)   


#        self.Combine_plot_button = Button(self.tab_Combine, text="Combine spectra", width=20)
#        self.Combine_plot_button.grid(row=2, column=0, sticky=W)  
#        self.Combine_plot_button.bind("<ButtonRelease-1>", self.Combine)  
       
        # Add the frame which contains the plot of the spectra 
        self.frame_Combine=Frame(self.tab_Combine, width=1000, height=self.Process_Height)
        self.frame_Combine.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_Combine.grid_propagate(False)

        self.frame_Combine.grid_rowconfigure(0, weight=1)
        self.frame_Combine.grid_columnconfigure(0, weight=1)

        self.frame_Combine.mytext_Extract = Text(self.frame_Combine, state="disabled")
        self.frame_Combine.mytext_Extract.place(x=10, y=10, height=990, width=390)  


        self.menu_Combine = StringVar(self.tab_Combine)
        self.menu_Combine.set("Median") # default value
#
        # Add a menu to select the combination method (Median or Mean)
        self.w = OptionMenu(self.tab_Combine, self.menu_Combine, "Median", "Mean")
        self.w.bind("<ButtonRelease-1>", self.Update_Spectra)
        self.w.grid(column=0, row=3, sticky=W)

        # Add a label to choice the name of the output file 
        Label(self.tab_Combine, text='Combine spectrum name:').grid(row=2,column = 0)
        self.CombineName = Entry(self.tab_Combine) 
        self.CombineName.grid(row=2, column=1) 
        self.CombineName.insert(END,'OutSpectrum.spec')





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
       
        self.frame_WavCal=Frame(self.tab_WavCal, width=1000, height=self.Process_Height)
        self.frame_WavCal.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_WavCal.grid_propagate(False)

        self.frame_WavCal.grid_rowconfigure(0, weight=1)
        self.frame_WavCal.grid_columnconfigure(0, weight=1)

        self.frame_WavCal.mytext_WaCal = Text(self.frame_WavCal, state="disabled")
        self.frame_WavCal.mytext_WaCal.place(x=10, y=10, height=990, width=390)  
        
        Label(self.tab_WavCal, text='Wavelength calibrated spectrum name:').grid(row=4,column = 0)
        self.WVName = Entry(self.tab_WavCal) 
        self.WVName.grid(row=4, column=1) 
        self.WVName.insert(END,'OutSpectrum.spec')


    def TellCorr_GUI(self):
        
        self.TellCorr_load_button = Button(self.tab_TellCorr, text="Load asteroid spectrum", width=20)
        self.TellCorr_load_button.grid(row=0, column=0, sticky=W)  
        self.TellCorr_load_button.bind("<ButtonRelease-1>", self.load_file_TellCorr)    

        self.TellCorr_Spec_load_button = Button(self.tab_TellCorr, text="Load standard spectrum", width=20)
        self.TellCorr_Spec_load_button.grid(row=2, column=0, sticky=W)  
        self.TellCorr_Spec_load_button.bind("<ButtonRelease-1>", self.load_file_TellCorr_Spec)  
 
        self.TellCorr_run_button = Button(self.tab_TellCorr, text="Correct", width=20)
        self.TellCorr_run_button.grid(row=3, column=0, sticky=W)  
        self.TellCorr_run_button.bind("<ButtonRelease-1>", self.TellCorr)   
       
        self.frame_TellCorr=Frame(self.tab_WavCal, width=1000, height=self.Process_Height)
        self.frame_TellCorr.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_TellCorr.grid_propagate(False)

        self.frame_TellCorr.grid_rowconfigure(0, weight=1)
        self.frame_TellCorr.grid_columnconfigure(0, weight=1)

        self.frame_TellCorr.mytext_TellCorr = Text(self.frame_TellCorr, state="disabled")
        self.frame_TellCorr.mytext_TellCorr.place(x=10, y=10, height=990, width=390)  
        
        Label(self.tab_TellCorr, text='Output spectrum name:').grid(row=4,column = 0)
        self.SpName = Entry(self.tab_TellCorr) 
        self.SpName.grid(row=4, column=1) 
        self.SpName.insert(END,'OutSpectrum.spec')




############################################################################## 
# SP_TellCorr
##############################################################################


    def load_file_TellCorr(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_TellCorr=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_TellCorr.append(elem)
        
#        self.read_all_fits(self.files_TellCorr)
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

    def load_file_TellCorr_Spec(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_TellCorr_Spec=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_TellCorr_Spec.append(elem)
        
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



    def TellCorr(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        print('python ' + Pipe_Path + '/SP_TellCorr.py ' + " ".join(self.files_TellCorr)  + ' -s ' + " ".join(self.files_TellCorr_Spec)  + ' -i Deveny')
        os.system('python ' + Pipe_Path + '/SP_TellCorr.py ' + " ".join(self.files_TellCorr)  + ' -s ' + " ".join(self.files_TellCorr_Spec) + ' -i Deveny')
        self.now = datetime.datetime.now()
        module_logger_WavCal.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Wavelength calibration done' )
        
        return None



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

        self.index = 0        
        im = self.images[self.index]
        self.tkimage.paste(im)
        
        # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)           

            
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
        print('python ' + Pipe_Path + '/SP_WavCal.py ' + " ".join(self.files_WavCal_Spec)  + '-a ' + " ".join(self.files_WavCal)  + ' -o ' + os.path.split(self.files_WavCal[0])[0] + '/' + self.WVName.get())
        os.system('python ' + Pipe_Path + '/SP_WavCal.py ' + " ".join(self.files_WavCal_Spec)  + ' -a ' + " ".join(self.files_WavCal) + ' -m auto' + ' -o ' + os.path.split(self.files_WavCal[0])[0] + '/' + self.WVName.get())
        self.now = datetime.datetime.now()
        module_logger_WavCal.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Wavelength calibration done' )
        
        return None




############################################################################## 
# SP_Combine
##############################################################################


    def Update_Spectra(self, event):
        # Function that allow to update the plot of the combined spectrum 
        if self.files_Combine:
            # Retrieve the method selected by the user
            Method = self.menu_Combine.get()
            
            if Method == 'Median': # Use median to combine the spectra
                self.Combined = np.nanmedian(self.spec,axis=0)
            if Method == 'Mean': # Use mean to combine the spectra
                self.Combined = np.nanmean(self.spec,axis=0)
            
            # Plot the combine spectrum
            self.ax.plot(self.Combined[0])
            self.canvas.draw()
            
            
            


    def Plot_Spectra(self):        
        # Function that allows to plot the combined spectrum for the first time
        
        if self.files_Combine: #  Check if there are spectra to combine
            # Initialize a few variables 
            self.spec = []
            self.FileToCombine = []
            
            # Loop to retrieve only the spectra selected by the user to combine them 
            for elem,elem2 in zip(self.files_Combine,self.fname):
                if self.Spec_List[elem2].get():
                    self.spec.append(np.loadtxt(elem).transpose())
                    self.FileToCombine.append(elem2)
            
            self.spec = np.array(self.spec)
            
            # Retrieve the method selected by the user to combine the spectra
            Method = self.menu_Combine.get()
            
            # Combine the spectra using the selected method
            if Method == 'Median':
                self.Combined = np.nanmedian(self.spec,axis=0)
            if Method == 'Mean':
                self.Combined = np.nanmean(self.spec,axis=0)
            
            # Clean the figure axes
            self.ax.remove()
            self.ax = self.fig.add_axes([0, 0, 1, 1])
            
            # Loop to plot the individual spectra using dots
            for elem in self.spec:
                self.ax.plot(elem[0,:],'.')

            # Plot the combined spectrum using a continuous line
            self.ax.plot(self.Combined[0])
#            print(self.Combined)
            
            # Use customs axes limit to avoid outliers
            self.ax.set_ylim([0,np.nanmedian(self.Combined[0])+2.5*np.nanstd(self.Combined[0])])
            
            # 
            self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

            self.canvas.draw()



    def load_file_Combine(self, event):
        # Function to load the spectra to combine
        
        # Open a file dialog
        self.fname = askopenfilenames()
        
        # Get the current time for logging
        self.now = datetime.datetime.now()
        
        # Initialize a few variables
        self.files_Combine=[]
        self.Spec_List = {}
        
        for elem in self.fname: # Loop over the selected files
            # Logging process
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1] )
            # Append the selected files into one variable  # why do that, is self.fname already that?
            self.files_Combine.append(elem)
            # Dictionarry for the check boxes 
            self.Spec_List[elem] = 0
            
#        print(self.Spec_List[self.fname[0]])
#        self.canvas.delete("all")
#        self.canvas.forget()



        # This loop allows to create the check boxes based on the number of spectra loaded
        for idx,elem in enumerate(self.Spec_List):
            self.Spec_List[elem] = IntVar(value=1)
            self.Check_Select_Spec = Checkbutton(self.tab_Combine, text="Spectrum " + str(idx+1), variable=self.Spec_List[elem],command=self.Plot_Spectra).grid(column = 2,row=idx)

        print(self.Spec_List)

        self.Plot_Spectra()
        return None



    def Combine(self, event):
        
        print(os.path.split(self.files_Combine[0])[0])
        
        now = datetime.datetime.now()
        module_logger.info(now)
        SP_Combine.Combine_Spectra(self.FileToCombine,os.path.split(self.files_Combine[0])[0] + '/' + self.CombineName.get(),False)
#        os.system('python ' + Pipe_Path + '/SP_Combine.py ' + " ".join(self.files_Combine) + ' -o ' + os.path.split(self.files_Combine[0])[0] + '/' + self.CombineName.get())
        self.now = datetime.datetime.now()
        module_logger_Extract.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Spectra combined' )
        
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
        
#               select first image
        
        self.index = 0
        im = self.images[self.index]
        self.tkimage.paste(im)
        
        
        # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)
        
        
            
        return None

    def load_Bckg(self, event):
        # Ask the user to select one background file
        self.fname = askopenfilename()
        
        # Get current time for logging purpuses 
        self.now = datetime.datetime.now()
        
        self.file_Bckg=self.fname
        module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading background file ' + self.file_Bckg.split('/')[-1])
        
        return None


    def Extract(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        
        print(self.Live_Extract.get())
        print(bool(self.Live_Extract.get()))
        SP_Extract.Extract_Spectrum(self.files_Extract,False,self.LocSpec_file,False,int(self.FWHM.get()),None,Live = self,Live2 = bool(self.Live_Extract.get()), Bckg = self.file_Bckg)
#        os.system('python ' + Pipe_Path + '/SP_Extract.py ' + " ".join(self.files_Extract))
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

            median = np.nanmedian(imgdat[int(imgdat.shape[1]*0.25):
                                         int(imgdat.shape[1]*0.75),
                                         int(imgdat.shape[0]*0.25):
                                         int(imgdat.shape[0]*0.75)])
            std    = np.nanstd(imgdat[int(imgdat.shape[1]*0.25):
                                      int(imgdat.shape[1]*0.75),
                                      int(imgdat.shape[0]*0.25):
                                      int(imgdat.shape[0]*0.75)])

            
            where_are_NaNs = np.isnan(imgdat)
            imgdat[where_are_NaNs] = 0
    
            imgdat = old_div(np.clip(imgdat, median-0.5*std,
                                median+0.5*std),(old_div(std,256)))
            
            imgdat = imgdat - np.nanmin(imgdat)

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

        self.index = 0        
        im = self.images[self.index]
        self.tkimage.paste(im)
        
        # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)           
        return None

    def BckgSub(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        SP_BckgSub.BckgSub(self.files_bckgsub,True,'auto','Bckg',self.LocSpec_file,'True',Area = [0,0],test = self,Live = bool(self.Live_BckgSub.get()))
#        os.system('python ' + Pipe_Path + '/SP_BckgSub.py ' + " ".join(self.files_bckgsub))
        self.now = datetime.datetime.now()
        module_logger_Bckgsub.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None


    def load_LocSpec(self,event):
        self.fname = askopenfilenames()
        
        self.loc_spec_files = []
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.loc_spec_files.append(elem)
        print(self.loc_spec_files)
        self.LocSpec_file = self.loc_spec_files[0].split('.fits')[0]
        print(self.LocSpec_file)
        
        


############################################################################## 
# SP_DetectSpectra
##############################################################################

    def load_file_detectspectra(self, event):
        
        # open a file dialog
        self.fname = askopenfilenames()
        
        # Get current time for logging
        self.now = datetime.datetime.now()
        
        # Initialize a few variables for the loop
        self.files_detectspectra=[]
        for elem in self.fname: # Loop over the selected files
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_detectspectra.append(elem)
            # Dictionarry for the check boxes 


        self.read_all_fits(self.files_detectspectra)
        print(self.images[0])
        
        self.canvas.delete("all")
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
              # select first image

        self.index = 0        
        im = self.images[self.index]
        self.tkimage.paste(im)
        
        # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)           
        return None
    
    def Detect_Spectra(self,event):
        
        now = datetime.datetime.now()
        module_logger.info(now)
        SP_DetectSpectra.Detect_Spectra(self.files_detectspectra,self.Output_Name_Det.get(),float(self.Det_Sens.get()))
        self.now = datetime.datetime.now()
        module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        self.read_all_fits([self.Output_Name_Det.get() + '.fits'])
        print(self.images[0])
        
        self.canvas.delete("all")
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)

        OffFile = self.Output_Name_Det.get() + '_Offset.txt'        
        with open(OffFile,'r') as f:
            Offset = simplejson.load(f)
        SpecFile = self.Output_Name_Det.get() + '_Spec_loc.txt'    
        with open(SpecFile,'r') as f:
            Spec_Loc = simplejson.load(f)
            
        
        for elem in Spec_Loc:
            self.canvas.create_rectangle(0,elem*self.zoom, self.tkimage.width(),elem*self.zoom -1,outline="#fb0", fill="#fb0")
            
        
                
        
        
        
        return None        
        

############################################################################## 
# SP_Subtract
##############################################################################

    def load_file_subtract(self, event):
        
        # open a file dialog
        self.fname = askopenfilenames()
        
        # Get current time for logging
        self.now = datetime.datetime.now()
        
        # Initialize a few variables for the loop
        self.files_subtract=[]
        self.fits_list = {}
        for elem in self.fname: # Loop over the selected files
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_subtract.append(elem)
            # Dictionarry for the check boxes 
            self.fits_list[elem] = 0 

        self.Check_Select_Spec = []
        self.image_one = StringVar(value =self.files_subtract[0])
        self.image_two = StringVar(value =self.files_subtract[1])
        
        Label(self.tab_Subtract, text='Image 1 -').grid(row=0,column = 2)
        Label(self.tab_Subtract, text='Image 2',justify='right',anchor=W).grid(row=0,column = 3,sticky="nesw")


        for idx,elem in enumerate(self.files_subtract):
            self.Check_Select_Spec.append(Radiobutton(self.tab_Subtract,value =elem, text="", variable=self.image_one,command=self.Display_Subtract_images).grid(column = 2,row=idx+1))
            self.Check_Select_Spec.append(Radiobutton(self.tab_Subtract,value =elem, text=elem.split("/")[-1], variable=self.image_two,command=self.Display_Subtract_images).grid(column = 3,row=idx+1))
            
            
        self.Display_Subtract_images()    


        # Clear the the previous canvas to display the new fits files
        self.canvas.delete("all")
        self.canvas.forget()
        
        # Display the new fits files
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
        
        # select first image
        self.index = 0       
        im = self.images[self.index]
        self.tkimage.paste(im)

        # events (allow the uses to click on the images) 
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)


        return None
    
    def Subtract_images(self,event):
        Image1 = self.image_one.get()
        Image2 = self.image_two.get()
        
        hdulist1 = fits.open(Image1)
        hdulist2 = fits.open(Image2)
                
        Image1_data = hdulist1[0].data 
        Image2_data = hdulist2[0].data 

        hdulist1[0].data = imgdat = Image1_data-Image2_data       
        
        hdulist1.writeto(self.Output_Name_Sub.get(),overwrite = True)
        

    def Display_Subtract_images(self):
        
        Image1 = self.image_one.get()
        Image2 = self.image_two.get()
        
        hdulist1 = fits.open(Image1)
        hdulist2 = fits.open(Image2)
                
        Image1_data = hdulist1[0].data 
        Image2_data = hdulist2[0].data 
        
        where_are_NaNs = np.isnan(Image1_data)
        Image1_data[where_are_NaNs] = 0

        where_are_NaNs = np.isnan(Image2_data)
        Image2_data[where_are_NaNs] = 0
        
        
        imgdat = Image1_data-Image2_data
    
        median = np.nanmedian(imgdat[int(imgdat.shape[1]*0.25):
                                     int(imgdat.shape[1]*0.75),
                                     int(imgdat.shape[0]*0.25):
                                     int(imgdat.shape[0]*0.75)])
        std    = np.nanstd(imgdat[int(imgdat.shape[1]*0.25):
                                  int(imgdat.shape[1]*0.75),
                                  int(imgdat.shape[0]*0.25):
                                  int(imgdat.shape[0]*0.75)])


        imgdat = old_div(np.clip(imgdat, median-0.5*std,
                            median+0.5*std),(old_div(std,256)))
        imgdat = imgdat - np.min(imgdat)

        imgdat = interp.zoom(imgdat, self.zoom)
        
        
        self.images = []
        
        self.images.append(Image.fromarray(imgdat))      
        
 
        self.canvas.delete("all")
        self.canvas.forget()
        
        # Display the new fits files
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)  
    

        # events (allow the uses to click on the images) 
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)
       
    

############################################################################## 
# SP_Cosmcorr
##############################################################################


    def Pop_Cosmcorr_tree(self,elem):
        hdulist = fits.open(elem)
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header['OBSTYPE']
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.CosmCorr_Tree = self.tree_CosmCorr.insert("",
                                    'end',
                                    "",
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))


    def load_file_cosmcorr(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_cosmcorr=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_cosmcorr.append(elem)
            self.Pop_Cosmcorr_tree(elem)

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
 
        self.index = 0       
        im = self.images[self.index]
        self.tkimage.paste(im)

                # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)

        return None

    def Cosmcorr(self, event):


        now = datetime.datetime.now()
        module_logger.info(now)
        
        files_to_use = []
        Children = self.tree_CosmCorr.get_children()
        for elem in Children:
            if self.tree_CosmCorr.tag_has("checked", elem):
                files_to_use.append(self.tree_CosmCorr.item(elem)['text'])        
        
        SP_CosmCorr.Cosmic(files_to_use,False)
        
#        os.system('python ' + Pipe_Path + '/SP_CosmCorr.py ' + " ".join(self.files_cosmcorr))
        self.now = datetime.datetime.now()
        module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None

############################################################################## 
# SP_Preproc
##############################################################################

    def Pop_Prep_bias_tree(self,elem):
        hdulist = fits.open(elem)
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header['OBSTYPE']
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Bias_Tree = self.tree_preproc_MBIAS.insert("",
                                    'end',
                                    "",
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))

    def Pop_Prep_flat_tree(self,elem):
        hdulist = fits.open(elem)
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header['OBSTYPE']
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Bias_Tree = self.tree_preproc_MFLAT.insert("",
                                    'end',
                                    "",
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))



    def Pop_Prep_toproc_tree(self,elem):
        hdulist = fits.open(elem)
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header['OBSTYPE']
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Toproc_Tree = self.tree_preproc_TOPROC.insert("",
                                    'end',
                                    "",
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))
        
    def Create_Folder_toproc_tree(self,files,foldername):
        Toproc_Tree_Folder = self.tree_preproc_TOPROC.insert("",
                                                            'end'
                                                             "",
                                                             text = foldername)
        for elem in files:
            hdulist = fits.open(elem)
            ExpTime = hdulist[0].header['EXPTIME']
            Type = hdulist[0].header['OBSTYPE']
            Date = hdulist[0].header['DATE-OBS']
            Mean = np.mean(hdulist[0].data)
            Median = np.median(hdulist[0].data)
            Std = np.std(hdulist[0].data)            
            Toproc_Tree = self.tree_preproc_TOPROC.insert(Toproc_Tree_Folder,
                                                          'end',
                                                          "",
                                                          text= elem.split('/')[-1],
                                                          values=(Type,
                                                                  str(ExpTime),
                                                                  str("{0:.1f}".format(Mean)),
                                                                  str("{0:.0f}".format(Median)),
                                                                  str("{0:.1f}".format(Std)),
                                                                  Date))                                                          
            


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
        
        self.index = 0;
        im = self.images[self.index]
        self.tkimage.paste(im)

                # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)

        return None


    def Preproc(self, event):
        now = datetime.datetime.now()
        module_logger_Preproc.info(now)
        
        files_to_use = []
        Children = self.tree_preproc_TOPROC.get_children()
        for elem2 in Children:
            if not self.tree_preproc_TOPROC.tag_has("unchecked", elem2):
                Children2 = self.tree_preproc_TOPROC.get_children(elem2)
                for elem in Children2:
                    if self.tree_preproc_TOPROC.tag_has("checked", elem):
                        print(self.tree_preproc_TOPROC.item(elem)['text'])
                        files_to_use.append(self.tree_preproc_TOPROC.item(elem)['text']) 
             
        Hold = False
        Children = self.tree_preproc_MBIAS.get_children()
        Bias_to_use = []
        Bias_Selected = False
        for elem in Children:
            if self.tree_preproc_MBIAS.tag_has("checked",elem):
                if not Bias_Selected:
                    Bias_to_use = self.tree_preproc_MBIAS.item(elem)['text']
                    Bias_Selected = True
                else:
                    print('More than one Master bias selected!')
                    Hold = True
 
        Children = self.tree_preproc_MFLAT.get_children()
        Flat_to_use = []
        Flat_Selected = False
        for elem in Children:
            if self.tree_preproc_MFLAT.tag_has("checked",elem):
                if not Flat_Selected:
                    Flat_to_use = self.tree_preproc_MFLAT.item(elem)['text']
                    Flat_Selected = True
                else:
                    print('More than one Master flat selected!')
               
                
        if not Hold:
            SP_Preproc.Preproc(files_to_use,Bias_to_use,Flat_to_use,False,'Procc',False)
            
    #        os.system('python ' + Pipe_Path + '/SP_Preproc.py ' + " ".join(self.files_preproc) + ' -b ' + self.MasterBiasName.get() + ' -f ' + self.MasterFlatName.get())
            self.now = datetime.datetime.now()
            module_logger_Flat.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'End of flat processing' )
        
        
        return None



############################################################################## 
# SP_Flat
##############################################################################

    def Pop_flat_tree(self,elem):
        hdulist = fits.open(elem)
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header['OBSTYPE']
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.nanmean(hdulist[0].data)
        Median = np.nanmedian(hdulist[0].data)
        Std = np.nanstd(hdulist[0].data)
        self.Flat_Tree = self.tree_flat.insert("",
                                    'end',
                                    "",
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))


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

        self.index = 0        
        im = self.images[self.index]
        self.tkimage.paste(im)

                # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)

        return None


    def Flat(self, event):
        now = datetime.datetime.now()
        module_logger_Flat.info(now)
        
        files_to_use = []
        Children = self.tree_flat.get_children()
        for elem in Children:
            if self.tree_flat.tag_has("checked", elem):
                files_to_use.append(self.tree_flat.item(elem)['text'])
        
        SP_Flat.Create_Flat(files_to_use,self.MasterFlatName_Flat.get(),True,self.MasterBiasName.get(),'none',False)
        #os.system('python ' + Pipe_Path + '/SP_Flat.py ' + " ".join(self.files_flat) + ' -b ' + self.MasterBiasName.get() + ' -o ' + self.MasterFlatName.get())
        self.now = datetime.datetime.now()
        module_logger_Flat.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'End of flat processing' )
        
        return None



############################################################################## 
# SP_Bias
##############################################################################
    def Pop_bias_tree(self,elem):
        hdulist = fits.open(elem)
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header['OBSTYPE']
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.nanmean(hdulist[0].data)
        Median = np.nanmedian(hdulist[0].data)
        Std = np.nanstd(hdulist[0].data)
        self.Bias_Tree = self.tree_bias.insert("",
                                    'end',
                                    "",
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))




    def load_file_bias(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_bias=[]
        self.files_bias_name = []
        
        # Clear bias tree
        self.tree_bias.delete(*self.tree_bias.get_children())
        
        for elem in self.fname:
            module_logger_Bias.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_bias.append(elem)
            self.files_bias_name.append(elem.split('/')[-1])
            # Populate bias tree
            self.Pop_bias_tree(elem)

        Children = self.tree_bias.get_children()
        for elem in Children:
            self.tree_bias.change_state(elem, "checked")

        self.menu_BiasFiles = StringVar(self.frame_fits)
        self.menu_BiasFiles.set(self.files_bias_name[0]) # default value  

        self.BiasFiles = OptionMenu(self.frame_fits, self.menu_BiasFiles, *self.files_bias_name)
        self.BiasFiles.pack()   

        self.read_all_fits(self.files_bias)
        
        self.canvas.delete("all")
        self.canvas.forget()
        
        
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(self.frame_fits, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.canvas.pack()
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)       
              # select first image
        
        self.index = 0
        im = self.images[self.index]
        self.tkimage.paste(im)
        
                        # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)


        return None


    def Bias(self, event):
        now = datetime.datetime.now()
        module_logger_Bias.info(now)
        
        files_to_use = []
        Children = self.tree_bias.get_children()
        for elem in Children:
            if self.tree_bias.tag_has("checked", elem):
                files_to_use.append(self.tree_bias.item(elem)['text'])        
        
        SP_Bias.Create_Bias(files_to_use,self.MasterBiasName_Bias.get(),False,False,self.menu_BiasCombine.get())
        
        
#        os.system('python ' + Pipe_Path + '/SP_Bias.py ' + " ".join(self.files_bias) + ' -o ' + self.MasterBiasName.get())
        self.now = datetime.datetime.now()
        module_logger_Bias.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'End of bias processing' )
        
        return None
############################################################################## 
# SP_Prepare
##############################################################################


    def load_file_prepare(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.tree_prepare.delete(*self.tree_prepare.get_children())
        self.files_prepare=[]
        for elem in self.fname:
            header = fits.getheader(elem) # get the header without opening data
            Type = header['OBSTYPE']
            ExpTime = header['EXPTIME']
            Date = header['DATE-OBS']
            Grat = header['GRATING']
            GratAng = header['WAVELENG']
            Target = header['OBJECT']
            RotAng = header['CRPA']
            self.Prepare_Tree = self.tree_prepare.insert("",
                                    'end',
                                    "",
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            Date,
                                            str(Grat),
                                            str(GratAng),
                                            str(Target),
                                            str(RotAng)))
            self.files_prepare.append(elem.split('/')[-1])

        Children = self.tree_prepare.get_children()
        for elem in Children:
            self.tree_prepare.change_state(elem, "checked")

        return None


    def Prepare(self,event):
        now = datetime.datetime.now()
        module_logger.info(now)

        files_to_prepare = []
        Children = self.tree_prepare.get_children()
        for elem in Children:
            if self.tree_prepare.tag_has("checked", elem):
                files_to_prepare.append(self.tree_prepare.item(elem)['text'])

        SP_Prepare.Prepare(files_to_prepare)        
        return None

############################################################################## 
# Front tab
##############################################################################    
    
    def load_file(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files=[]
        for elem in self.fname:
#            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
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
#        module_logger.info(self.now)
        os.system('python ' + Pipe_Path + '/SP_DisplayFits.py ' + " ".join(self.fname))
        
        return None



### Event listener
    def left_click(self,event):
        print(event.x,event.y)
        x, y = old_div(event.x,self.zoom), old_div(event.y,self.zoom) 
        print(x,y)
        self.x = x
        self.y = y

    def right_click(self, event):
        """ select source """
        x, y = old_div(event.x,self.zoom), old_div(event.y,self.zoom)

    def key(self, event):
        """ keyboard events """
        if event.char == 'a':
            # previous frame
            self.nextframe()
        elif event.char == 'q':
            # quit
            self.top.quit()
            
        elif event.char == 'f':
            self.Fit_Spectrum(self.data[self.index],self.x,self.y,30)


    def nextframe(self):
        """ display frame using iterator i"""

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
        
                # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)
        
        
        self.index += 1
        if self.index>=len(self.images):
            self.index -= len(self.images)
        print(self.index)
        im = self.images[self.index]
        self.tkimage.paste(im)



    def Fit_Spectrum(self,image,x,y,box):


        def MOFFAT(x, *p):
            S0,S1,a, b, x0 = p
            bb = S0+S1/(1+(x-x0)**2/a**2)**b
            
            return bb
        
        X = []
        Y = []
        
        print(x,y)
        D = image[int(y)-box:int(y)+box,int(x)-box:int(x)+box]
        
        Line = np.median(D,axis=1)
        
        
        xs = range(len(Line))
        ys = Line
        
        max_index, max_value = max(enumerate(Line), key=operator.itemgetter(1))
        
        p0 = [0,max_value,2.24,1.41,max_index] 
        Res = []

        coeff, var_matrix = curve_fit(lambda x, S0,S1,a, b,x0: MOFFAT(x, S0,S1,a, b,x0),xs,ys,p0,maxfev = 100000)
        fit = MOFFAT(xs, *coeff)
        print(fit)
        
        FWHM = 2*coeff[2]*np.sqrt(2**(1/coeff[3])-1)
        
        region = sum(Line>(coeff[0]+50))/2
        
        self.Yloc = max_index+int(y)-box
        print(self.Yloc)
        
        fig, ax = plt.subplots()
        plt.plot(Line,Label = 'Spectrum')
        plt.plot(fit,Label = 'Fit')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.1, 0.9,'FWHM = ' + str(round(FWHM,1)) + ' pix',transform=ax.transAxes,fontsize=14,verticalalignment='top', bbox=props)
        ax.text(0.1, 0.8,'Trace = ' + str(region) + ' pix',transform=ax.transAxes,fontsize=14,verticalalignment='top', bbox=props)
        plt.legend()
        plt.show()



    # Closing command 
    def on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
#            self.f.close()
            self.quit()
            self.destroy()




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
 