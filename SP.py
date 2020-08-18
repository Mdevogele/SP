#!/usr/bin/env python
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

try:
    from astropy.visualization import PercentileInterval, ZScaleInterval
    from astropy.io import fits
except: 
    print('Astropy module was not found. try: pip install astropy')
    
    
import logging
import datetime

from past.utils import old_div

import operator

import SP_Bias
import SP_Prepare
import SP_BckgSub
import SP_Extract
import SP_Combine
import SP_Flat
import SP_Preproc
import SP_WavCal
from SP_CheckInstrument import CheckInstrument


import SP_Toolbox as tb


from PIL import Image
from PIL import ImageTk
from PIL import ImageDraw

import argparse, shlex
import numpy as np

from glob import glob

from scipy.ndimage import interpolation as interp

import imageio as misc


# try:
#     from scipy import misc
# except:
#     import imageio as misc

from scipy.optimize import curve_fit


if sys.version_info <= (3,0):
    # if python 2
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    import matplotlib.backends.tkagg as tkagg
else: 
    # if python 3
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg

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



class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text, Lab):
        "Display text in tooltip window"
        self.text = text
#        if self.tipwindow or not self.text:
#            return
#        x, y, cx, cy = self.widget.bbox("insert")
#        x = x + self.widget.winfo_rootx() + 57
#        y = y + cy + self.widget.winfo_rooty() +27
#        self.tipwindow = tw = Toplevel(self.widget)
#        tw.wm_overrideredirect(1)
#        tw.wm_geometry("+%d+%d" % (x, y))
#        Lab['text'] = text
        print(Lab)
        Lab.set(text)
#        label = Label(tw, text=self.text, justify=LEFT,
#                      background="#ffffe0", relief=SOLID, borderwidth=1,
#                      font=("tahoma", "8", "normal"))
#        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()


class simpleapp_tk(Tk):
    def __init__(self,parent):
        Tk.__init__(self,parent)
        
        self.cdw = os.getcwd()
        
        # Closing handeler
        
        self.protocol("WM_DELETE_WINDOW", self.on_closing)    
        
        # Manage log text file
        self.now = datetime.datetime.now()
        self.logger = logging.getLogger("my logger")
#        self.f = open('Log_' + self.now.strftime("%Y_%m_%d_%H_%M_%S") +'.txt','w')
        logging.basicConfig(filename='Log_' + self.now.strftime("%Y_%m_%d_%H_%M_%S") +'.txt', filemode='a', format='%(asctime)s \t %(message)s',level=logging.INFO,datefmt='%Y-%m-%d %H:%M:%S')

        # zoom factor of the fits files.
        self.zoom = 0.5
        
        self.parent = parent

        # To use grid based layout
        self.grid()
        
        # Definition of the frame containing the processing tabs
        self.frame_process=Frame(self, width=1000, height=400)
        self.frame_process.grid(column=0, row=0)

        # Definition of the frame containing the fits images or graphs
        self.frame_Tips=Frame(self, width=1000, height=200)
        self.frame_Tips.grid(column=0, row=1)
        
        self.TipText = StringVar(self)
        self.Label_Tips = Label(self.frame_Tips, textvariable=self.TipText).pack()
        self.TipText.set('test')
        
        # Definition of the frame containing the fits images or graphs
        self.frame_fits=Frame(self, width=1000, height=200)
        self.frame_fits.grid(column=0, row=2)
        
        # Definition of the frame containing the plot toolbar
        self.frame_tb = Frame(self,width = 400,heigh = 50)        
        self.frame_tb.grid(column=0, row=3)
        
        
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

        # SP_BckgSub tab
        self.tab_BckgSub = Frame(self.TabControl)
        self.TabControl.add(self.tab_BckgSub,text='Background subtraction')
        self.tab_index.update({6:{'Name':'Bckgsub','Clicked':False}})

        # SP_Extract tab
        self.tab_Extract = Frame(self.TabControl)
        self.TabControl.add(self.tab_Extract,text='Spectra extraction')        
        self.tab_index.update({7:{'Name':'Extract','Clicked':False}})

        # SP_Combine tab
        self.tab_Combine = Frame(self.TabControl)
        self.TabControl.add(self.tab_Combine,text='Combine spectra')
        self.TabControl.bind(self.tab_Combine,"<ButtonRelease-1>", self.load_file)    
        self.tab_index.update({8:{'Name':'Combine','Clicked':False}})

        # SP_WavCal tab
        self.tab_WavCal = Frame(self.TabControl)
        self.TabControl.add(self.tab_WavCal,text='Wavelength calibration')   
        self.tab_index.update({9:{'Name':'Wavcal','Clicked':False}})

        # SP_TellCorr tab
        self.tab_TellCorr = Frame(self.TabControl)
        self.TabControl.add(self.tab_TellCorr,text='Telluric correction')   
        self.tab_index.update({10:{'Name':'TellCorr','Clicked':False}})
        
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
        
        
        # Initialization of some needed variables
        
        self.files_Combine = []
        self.Yloc = False
        
        
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
        self.Combine_GUI()
        self.WavCal_GUI()
        self.TellCorr_GUI()

        self.fig = plt.figure()
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.lastindex = []

        self.TabControl.bind("<<NotebookTabChanged>>", self.handle_tab_changed)


    def CreateToolTip(self,widget, text, Lab = 'lab'):
        toolTip = ToolTip(widget)
        def enter(event):
            toolTip.showtip(text,Lab)
        def leave(event):
            toolTip.hidetip()
        widget.bind('<Enter>', enter)
        widget.bind('<Leave>', leave)


    def handle_tab_changed(self,event):

        selection = event.widget.select()
        tab = event.widget.index("current")
        
        if self.lastindex != self.TabControl.index(self.tab_Combine):
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
            self.files_bias_name = glob('*BIAS*fits')
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
            self.files_flat_name = glob('*FLAT*fits')

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
            for elem in self.files_pre_name:
                header = fits.getheader(elem)
                
                try: 
                    if header['PROCTYPE'] == 'MASTER FLAT':
                        Flats.append(elem)
                    if header['PROCTYPE'] == 'MASTER BIAS':
                        Biases.append(elem)
                    if header[obsparam['obstype']] == 'OBJECT' and header['PROCTYPE'] == 'Prepared':
                        ToProc.append(elem)
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

            for elem in ToProc:
                self.Pop_Prep_toproc_tree(elem)           

            Children = self.tree_preproc_TOPROC.get_children()
            for elem in Children:
                self.tree_preproc_TOPROC.change_state(elem, "checked")



        # Action to be performed when WAVCAL tab is clicked for the first time 
        if self.tab_index[int(tab)]['Name'] == 'Wavcal' and self.tab_index[int(tab)]['Clicked'] == False:
            self.files_WavCal_name = glob('*ARCS*fits')

            self.Pop_WavCal_tree(self.files_WavCal_name)
         
#            Children = self.tree_WavCal.get_children()
#            for elem in Children:
#                self.tree_WavCal.change_state(elem, "checked")

            self.files_WavCal = []
            for elem in self.files_WavCal_name:
                self.files_WavCal.append(self.cdw + '/' + elem)





        self.lastindex = event.widget.index("current")
        self.tab_index[int(tab)]['Clicked'] = True


    def Prepare_GUI(self):
        
        
        self.frame_Prepare_Button=Frame(self.tab_Prepare, width=100, height=400)
        self.frame_Prepare_Button.grid(column=0, row=0)

        self.Prepare_load_button = Button(self.frame_Prepare_Button, text="Load files", width=20)
        self.Prepare_load_button.grid(row=0, column=0, sticky=W)  
        self.Prepare_load_button.bind("<ButtonRelease-1>", self.load_file_prepare)    
 
        self.Prepare_run_button = Button(self.frame_Prepare_Button, text="Run prepare", width=20)
        self.Prepare_run_button.grid(row=1, column=0, sticky=W)  
        self.Prepare_run_button.bind("<ButtonRelease-1>", self.Prepare)   
       
        self.frame_Prepare=Frame(self.tab_Prepare, width=900, height=400)
        self.frame_Prepare.grid(column=1, row=0)

        
        self.tree_prepare = CheckboxTreeview(self.frame_Prepare)
        self.tree_prepare["columns"]=("Type","ExpTime","Date","Grating","GratAng","Target","RotAng","Collfoc")
        self.tree_prepare.column("#0", width=250)
        self.tree_prepare.column("Type", width=100)
        self.tree_prepare.column("ExpTime", width=100)
        self.tree_prepare.column("Date", width=100)
        self.tree_prepare.column("Grating", width=100)
        self.tree_prepare.column("GratAng", width=100)
        self.tree_prepare.column("Target", width=100)
        self.tree_prepare.column("RotAng", width=100)
        self.tree_prepare.column("Collfoc", width=100)
        
        self.tree_prepare.heading("#0",text="File name",anchor=W)
        self.tree_prepare.heading("Type",text="Type",anchor=CENTER)
        self.tree_prepare.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_prepare.heading("Date", text="Date",anchor=W)
        self.tree_prepare.heading("Grating", text="Grating",anchor=W)
        self.tree_prepare.heading("GratAng", text="Grat. ang.",anchor=W)        
        self.tree_prepare.heading("Target", text="Target",anchor=W)   
        self.tree_prepare.heading("RotAng", text="Rot. ang.",anchor=W)   
        self.tree_prepare.heading("Collfoc", text="Coll. foc.",anchor=W)   
        
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

        self.Frame_Preproc_Button = Frame(self.tab_Preproc, width=200, height=400)
        self.Frame_Preproc_Button.grid(row = 0,column =0)
        
        self.Frame_Preproc_Info = Frame(self.tab_Preproc, width=800, height=400)
        self.Frame_Preproc_Info.grid(row = 0,column =1)  

        
        self.Preproc_load_button = Button(self.Frame_Preproc_Button, text="Load files", width=20)
        self.Preproc_load_button.grid(row=0, column=0, sticky=W)  
        self.Preproc_load_button.bind("<ButtonRelease-1>", self.load_file_preproc)    
        self.CreateToolTip(self.Preproc_load_button,'Tips: Load the file that you want to remove the bias and correct from the flat from',Lab = self.TipText)

        # Button for manually load the masterbias
        self.Preproc_MasterBias_Load_button = Button(self.Frame_Preproc_Button, text="Load master bias", width=20)
        self.Preproc_MasterBias_Load_button.grid(row=1, column=0, sticky=W)  
        self.Preproc_MasterBias_Load_button.bind("<ButtonRelease-1>", self.load_file_preproc_bias)   
        self.CreateToolTip(self.Preproc_MasterBias_Load_button,'Tips: Manually select the master bias',Lab = self.TipText)
        

        # Button for manually load the masterflat
        self.Preproc_MasterFlat_Load_button = Button(self.Frame_Preproc_Button, text="Load master flat", width=20)
        self.Preproc_MasterFlat_Load_button.grid(row=2, column=0, sticky=W)  
        self.Preproc_MasterFlat_Load_button.bind("<ButtonRelease-1>", self.load_file_preproc_flat)   
        self.CreateToolTip(self.Preproc_MasterFlat_Load_button,'Tips: Manually select the master flat',Lab = self.TipText)


 
        self.Preproc_run_button = Button(self.Frame_Preproc_Button, text="Process files", width=20)
        self.Preproc_run_button.grid(row=3, column=0, sticky=W)  
        self.Preproc_run_button.bind("<ButtonRelease-1>", self.Preproc)   
        self.CreateToolTip(self.Preproc_run_button,'Tips: Launch the preprocessing',Lab = self.TipText)


       
#        Label(self.Frame_Preproc_Button, text='Master flat name:').grid(row=4,column = 0)
#        self.MasterFlatName = Entry(self.Frame_Preproc_Button) 
#        self.MasterFlatName.grid(row=4, column=1) 
#        self.MasterFlatName.insert(END,'MasterFlat.fits')

#        Label(self.Frame_Preproc_Button, text='Master bias name:').grid(row=3,column = 0)
#        self.MasterBiasName = Entry(self.Frame_Preproc_Button) 
#        self.MasterBiasName.grid(row=3, column=1) 
#        self.MasterBiasName.insert(END,'MasterBias.fits')

        self.tree_preproc_MBIAS = CheckboxTreeview(self.Frame_Preproc_Info, height =3)
        self.tree_preproc_MBIAS["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_preproc_MBIAS.column("#0", width=400)
        self.tree_preproc_MBIAS.column("Type", width=50)
        self.tree_preproc_MBIAS.column("ExpTime", width=50)
        self.tree_preproc_MBIAS.column("Mean", width=70)
        self.tree_preproc_MBIAS.column("Median", width=70)
        self.tree_preproc_MBIAS.column("std", width=70)
        self.tree_preproc_MBIAS.column("Date", width=200)
        self.CreateToolTip(self.tree_preproc_MBIAS,'Tips: List of master biases. Select one file to use to for processing',Lab = self.TipText)
        
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
        self.CreateToolTip(self.tree_preproc_MFLAT,'Tips: List of master flats. Select one file to use to for processing',Lab = self.TipText)
        
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
        self.CreateToolTip(self.tree_preproc_TOPROC,'Tips: List of files to be processed. Only selected files will be processed',Lab = self.TipText)
        
        self.tree_preproc_TOPROC.heading("#0",text="File name",anchor=W)
        self.tree_preproc_TOPROC.heading("Type",text="Type",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("Median", text="Median",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("std", text="std",anchor=CENTER)
        self.tree_preproc_TOPROC.heading("Date", text="Date",anchor=W)
        
        self.tree_preproc_TOPROC.grid(column=1,row = 2)    







    def CosmCorr_GUI(self):

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

        


    def Extract_GUI(self):


        self.Frame_Extract_Button = Frame(self.tab_Extract, width=200, height=400)
        self.Frame_Extract_Button.grid(row = 0,column =0)
        
        self.Frame_Extract_Info = Frame(self.tab_Extract, width=800, height=400)
        self.Frame_Extract_Info.grid(row = 0,column =1)     

        
        self.Extract_load_button = Button(self.Frame_Extract_Button, text="Load files", width=20)
        self.Extract_load_button.grid(row=0, column=0, sticky=W)  
        self.Extract_load_button.bind("<ButtonRelease-1>", self.load_file_Extract)    
 
        self.Extract_run_button = Button(self.Frame_Extract_Button, text="Extract spectra", width=20)
        self.Extract_run_button.grid(row=1, column=0, sticky=W)  
        self.Extract_run_button.bind("<ButtonRelease-1>", self.Extract)   
  

        Label(self.Frame_Extract_Button, text='# of pixel to extract:').grid(row=3,column = 0)
        self.FWHM = Entry(self.Frame_Extract_Button) 
        self.FWHM.grid(row=3, column=1) 
        self.FWHM.insert(END,'6')

        self.Live_Extract = IntVar()
        self.Check_Live_Extract = Checkbutton(self.Frame_Extract_Button, text="Live processing", variable=self.Live_Extract).grid(column = 0,row=4)


        self.tree_Extract = CheckboxTreeview(self.Frame_Extract_Info)
        self.tree_Extract["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_Extract.column("#0", width=400)
        self.tree_Extract.column("Type", width=50)
        self.tree_Extract.column("ExpTime", width=50)
        self.tree_Extract.column("Mean", width=70)
        self.tree_Extract.column("Median", width=70)
        self.tree_Extract.column("std", width=70)
        self.tree_Extract.column("Date", width=200)
        
        self.tree_Extract.heading("#0",text="File name",anchor=W)
        self.tree_Extract.heading("Type",text="Type",anchor=CENTER)
        self.tree_Extract.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_Extract.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_Extract.heading("Median", text="Median",anchor=CENTER)
        self.tree_Extract.heading("std", text="std",anchor=CENTER)
        self.tree_Extract.heading("Date", text="Date",anchor=W)
        
        self.tree_Extract.grid(column=1,row = 0)   



    def Combine_GUI(self):
        
        self.Combine_load_button = Button(self.tab_Combine, text="Load files", width=20)
        self.Combine_load_button.grid(row=0, column=0, sticky=W)  
        self.Combine_load_button.bind("<ButtonRelease-1>", self.load_file_Combine)    
 
        self.Combine_run_button = Button(self.tab_Combine, text="Combine spectra", width=20)
        self.Combine_run_button.grid(row=1, column=0, sticky=W)  
        self.Combine_run_button.bind("<ButtonRelease-1>", self.Combine)   

       
        self.frame_Combine=Frame(self.tab_Combine, width=1000, height=400)
        self.frame_Combine.grid(column=2, row=0,rowspan=25,columnspan=10)

        self.frame_Combine.grid_propagate(False)

        self.frame_Combine.grid_rowconfigure(0, weight=1)
        self.frame_Combine.grid_columnconfigure(0, weight=1)

        self.frame_Combine.mytext_Extract = Text(self.frame_Combine, state="disabled")
        self.frame_Combine.mytext_Extract.place(x=10, y=10, height=990, width=390)  


        self.menu_Combine = StringVar(self.tab_Combine)
        self.menu_Combine.set("Median") # default value
#
        self.w = OptionMenu(self.tab_Combine, self.menu_Combine, "Median", "Mean")
        self.w.bind("<ButtonRelease-1>", self.Update_Spectra)
        self.w.grid(column=0, row=3, sticky=W)

        Label(self.tab_Combine, text='Combine spectrum name:').grid(row=2,column = 0)
        self.CombineName = Entry(self.tab_Combine) 
        self.CombineName.grid(row=2, column=1) 
        self.CombineName.insert(END,'OutSpectrum.spec')





    def WavCal_GUI(self):
        
        self.Frame_WavCal_Button = Frame(self.tab_WavCal, width=200, height=400)
        self.Frame_WavCal_Button.grid(row = 0,column =0)
        
        self.Frame_WavCal_Info = Frame(self.tab_WavCal, width=200, height=400)
        self.Frame_WavCal_Info.grid(row = 0,column =1)     

        
        self.WavCal_load_button = Button(self.Frame_WavCal_Button, text="Load arc files", width=20)
        self.WavCal_load_button.grid(row=0, column=0, sticky=W)  
        self.WavCal_load_button.bind("<ButtonRelease-1>", self.load_file_WavCal)    

        self.WavCal_Spec_load_button = Button(self.Frame_WavCal_Button, text="Load spec file", width=20)
        self.WavCal_Spec_load_button.grid(row=2, column=0, sticky=W)  
        self.WavCal_Spec_load_button.bind("<ButtonRelease-1>", self.load_file_WavCal_Spec)  
 
        self.WavCal_run_button = Button(self.Frame_WavCal_Button, text="Calibrate wavelength", width=20)
        self.WavCal_run_button.grid(row=3, column=0, sticky=W)  
        self.WavCal_run_button.bind("<ButtonRelease-1>", self.WavCal)   
 

       
        Label(self.Frame_WavCal_Button, text='Wavelength calibrated spectrum name:').grid(row=4,column = 0)
        self.WVName = Entry(self.Frame_WavCal_Button) 
        self.WVName.grid(row=4, column=1) 
        self.WVName.insert(END,'OutSpectrum.spec')

        
        self.menu_WavCal_Type = StringVar(self.Frame_WavCal_Button)
        self.menu_WavCal_Type.set("2D") # default value   
        
        Label(self.Frame_WavCal_Button, text='Type: ').grid(row=5,column = 0)
        self.WavCal_Type = OptionMenu(self.Frame_WavCal_Button, self.menu_WavCal_Type, "1D", "2D")
        self.WavCal_Type.grid(column=1, row=5, sticky=W)        


        self.menu_WavCal_Method = StringVar(self.Frame_WavCal_Button)
        self.menu_WavCal_Method.set("Auto") # default value   
        
        Label(self.Frame_WavCal_Button, text='Method: ').grid(row=6,column = 0)
        self.WavCal_Method = OptionMenu(self.Frame_WavCal_Button, self.menu_WavCal_Method, "Auto", "Template")
        self.WavCal_Method.grid(column=1, row=6, sticky=W)   



        
        self.tree_WavCal = CheckboxTreeview(self.Frame_WavCal_Info)
        self.tree_WavCal["columns"]=("Type","ExpTime","Grating","Lamp1","Lamp2","Lamp3","Lamp4")
        self.tree_WavCal.column("#0", width=400)
        self.tree_WavCal.column("Type", width=50)
        self.tree_WavCal.column("ExpTime", width=50)
        self.tree_WavCal.column("Grating", width=50)
        self.tree_WavCal.column("Lamp1", width=30)
        self.tree_WavCal.column("Lamp2", width=30)
        self.tree_WavCal.column("Lamp3", width=30)
        self.tree_WavCal.column("Lamp4", width=30)
        
        self.tree_WavCal.heading("#0",text="File name",anchor=W)
        self.tree_WavCal.heading("Type",text="Type",anchor=CENTER)
        self.tree_WavCal.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_WavCal.heading("Grating", text="Grating", anchor=CENTER)
        self.tree_WavCal.heading("Lamp1", text="Ar",anchor=CENTER)
        self.tree_WavCal.heading("Lamp2", text="Ne",anchor=CENTER)
        self.tree_WavCal.heading("Lamp3", text="Hg",anchor=CENTER)
        self.tree_WavCal.heading("Lamp4", text="Cd",anchor=W)
        
        self.tree_WavCal.grid(column=1,row = 0)   


        
        


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
       
        self.frame_TellCorr=Frame(self.tab_TellCorr, width=1000, height=400)
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


    def Pop_WavCal_tree(self,files_list):

        telescope, obsparam = CheckInstrument([files_list[0]])
        
        if telescope != 'NOT':
            self.tree_WavCal.delete(*self.tree_WavCal.get_children())
            Split = []
            for elem in files_list:
    #                self.Pop_WavCal_tree(elem)
                Split.append(int(elem.split('/')[-1].split('_')[1]))
              
            Conseq = tb.Get_Consecutive(Split)    
        else:
            Conseq = tb.Get_Consecutive(range(len(files_list)))
 
        for elem in Conseq:
            
            telescope, obsparam = CheckInstrument([files_list[elem[0][0]]])
            hdulist = fits.open(files_list[elem[0][0]])
            ExpTime = hdulist[0].header['EXPTIME']
            Type = hdulist[0].header[obsparam['obstype']]                
            Date = hdulist[0].header['DATE-OBS']
            Grating = hdulist[0].header[obsparam['grating']]
            
            if obsparam['HDR_lamp']:
                Lamp1_Text = obsparam['lamp1_name']
                Lamp2_Text = obsparam['lamp2_name']
                Lamp3_Text = obsparam['lamp3_name']
                Lamp4_Text = obsparam['lamp4_name']
    
                print(hdulist[0].header[obsparam['lamp1']])
                try:
                    if 'True' in hdulist[0].header[obsparam['lamp1']]:
                        Lamp1_OnOff = u'\u2713'
                    else:
                        Lamp1_OnOff = X
                except: 
                    Lamp1_OnOff = 'Ukn'
                    
                try:
                    if 'True' in hdulist[0].header[obsparam['lamp2']]:
                        Lamp2_OnOff = u'\u2713'
                    else:
                        Lamp2_OnOff = X
                except: 
                    Lamp2_OnOff = 'Ukn'
                    
                try:
                    if 'True' in hdulist[0].header[obsparam['lamp3']]:
                        Lamp3_OnOff = u'\u2713'
                    else:
                        Lamp3_OnOff = X
                except: 
                    Lamp3_OnOff = 'Ukn'
                    
                try:
                    if 'True' in hdulist[0].header[obsparam['lamp4']]:
                        Lamp4_OnOff = u'\u2713'
                    else:
                        Lamp4_OnOff = X
                except:
                    Lamp4_OnOff = 'Ukn'
            else:
                Lamp1_OnOff = 'Ukn'
                Lamp2_OnOff = 'Ukn'
                Lamp3_OnOff = 'Ukn'
                Lamp4_OnOff = 'Ukn'
                Lamp1_Text = 'Ukn'
                Lamp2_Text = 'Ukn'
                Lamp3_Text = 'Ukn'
                Lamp4_Text = 'Ukn'
                

            self.tree_WavCal.heading("Lamp1", text=Lamp1_Text,anchor=CENTER)
            self.tree_WavCal.heading("Lamp2", text=Lamp2_Text,anchor=CENTER)
            self.tree_WavCal.heading("Lamp3", text=Lamp3_Text,anchor=CENTER)
            self.tree_WavCal.heading("Lamp4", text=Lamp4_Text,anchor=W)
            
            #Create the folder containing the individual files
            Folder1 = self.tree_WavCal.insert("",
                                'end',
                                text= str(elem),
                                values=(Type,
                                        str(ExpTime),
                                        str(Grating),
                                        Lamp1_OnOff,
                                        Lamp2_OnOff,
                                        Lamp3_OnOff,
                                        Lamp4_OnOff))
                                
            # populate the individual files 
            print(self.files_WavCal_name)
            for elem2 in elem:
                files = self.files_WavCal_name[elem2[0]]
            
                hdulist = fits.open(files)
                ExpTime = hdulist[0].header['EXPTIME']
                Type = hdulist[0].header[obsparam['obstype']]
                Date = hdulist[0].header['DATE-OBS']
                Mean = np.mean(hdulist[0].data)
                Median = np.median(hdulist[0].data)
                Std = np.std(hdulist[0].data)
                self.WavCal_Tree = self.tree_WavCal.insert(Folder1,
                                'end',
                                text= 'bla',
                                values=(Type,
                                        str(ExpTime),
                                        str("{0:.1f}".format(Mean)),
                                        str("{0:.0f}".format(Median)),
                                        str("{0:.1f}".format(Std)),
                                        Date))

            Children = self.tree_WavCal.get_children()
            for elem in Children:
                self.tree_WavCal.change_state(elem, "checked")






    def load_file_WavCal(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_WavCal=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_WavCal.append(elem)
            self.files_WavCal_name.append(elem)
        
        self.Pop_WavCal_tree(self.files_WavCal)


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
        SP_WavCal.WavCal(self.files_WavCal_Spec,self.files_WavCal,os.path.split(self.files_WavCal[0])[0] + '/' + self.WVName.get(),False,self.menu_WavCal_Method.get(),250,False)
        self.now = datetime.datetime.now()
        module_logger_WavCal.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Wavelength calibration done' )
        
        return None




############################################################################## 
# SP_Combine
##############################################################################


    def Update_Spectra(self, event):
        print('bla')
        if self.files_Combine:

            Method = self.menu_Combine.get()
            
            if Method == 'Median':
                self.Combined = np.nanmedian(self.spec,axis=0)
            if Method == 'Mean':
                self.Combined = np.nanmean(self.spec,axis=0)
                
            self.ax.plot(self.Combined[0])
            self.canvas.draw()
            
            
            


    def Plot_Spectra(self):        
        if self.files_Combine:
            self.spec = []
            self.FileToCombine = []
#            plt.figure()
            for elem,elem2 in zip(self.files_Combine,self.fname):
                if self.Spec_List[elem2].get():
                    self.spec.append(np.loadtxt(elem).transpose()/np.nanmedian(np.loadtxt(elem).transpose()[0,1100:1600]))
                    self.FileToCombine.append(elem2)
            
            Method = self.menu_Combine.get()
            
            self.spec = np.array(self.spec)
            if Method == 'Median':
                self.Combined = np.nanmedian(self.spec,axis=0)
            if Method == 'Mean':
                self.Combined = np.nanmean(self.spec,axis=0)
            self.ax.remove()
            self.ax = self.fig.add_axes([0, 0, 1, 1])
            for elem in self.spec:
#                if self.Spec_List[elem2].get():
#                    print(elem2)
                    self.ax.plot(elem[0,:],'.')
            
            self.ax.plot(self.Combined[0])
            print(self.Combined)
            
            self.ax.set_ylim([np.nanmedian(self.Combined[0])-2*np.nanstd(self.Combined[0]),np.nanmedian(self.Combined[0])+1*np.nanstd(self.Combined[0])])
            
#            self.canvas = FigureCanvasTkAgg(fig, self.frame_fits)
#            self.canvas.show()
#            self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1) 
    
            toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame_tb)
            toolbar.update()
            self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

            self.canvas.draw()



    def load_file_Combine(self, event):
        
#        self.fig = plt.figure()
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_Combine=[]
        self.Spec_List = {}
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1] )
            self.files_Combine.append(elem)
            self.Spec_List[elem] = 0  #Dictionarry for the check boxes 
        print(self.Spec_List[self.fname[0]])
#        self.canvas.delete("all")
#        self.canvas.forget()


        for idx,elem in enumerate(self.Spec_List):
            self.Spec_List[elem] = IntVar(value=1)
            self.Check_Select_Spec = Checkbutton(self.tab_Combine, text="Spectrum " + str(idx+1), variable=self.Spec_List[elem],command=self.Plot_Spectra).grid(column = 2,row=idx)

        print(self.Spec_List)

        self.Plot_Spectra()
        return None



    def Combine(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        SP_Combine.Combine_Spectra(self.FileToCombine,os.path.split(self.files_Combine[0])[0] + '/' + self.CombineName.get(),False)
        np.savetxt(self.CombineName.get(),self.Combined.transpose())
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

    def Extract(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        
        print(self.Live_Extract.get())
        print(bool(self.Live_Extract.get()))
        SP_Extract.Extract_Spectrum(self.files_Extract,False,'bla',False,int(self.FWHM.get()),self.Yloc,Live = self,Live2 = bool(self.Live_Extract.get()))
#        os.system('python ' + Pipe_Path + '/SP_Extract.py ' + " ".join(self.files_Extract))
        self.now = datetime.datetime.now()
        module_logger_Extract.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None



        

############################################################################## 
# SP_Bckgsub
##############################################################################


    # GUI
    
    def BckgSub_GUI(self):

        self.Frame_BckgSub_Button = Frame(self.tab_BckgSub, width=200, height=400)
        self.Frame_BckgSub_Button.grid(row = 0,column =0)
        
        self.Frame_BckgSub_Info = Frame(self.tab_BckgSub, width=800, height=400)
        self.Frame_BckgSub_Info.grid(row = 0,column =1)     

        
        self.BckgSub_load_button = Button(self.Frame_BckgSub_Button, text="Load files", width=20)
        self.BckgSub_load_button.grid(row=0, column=0, sticky=W)  
        self.BckgSub_load_button.bind("<ButtonRelease-1>", self.load_file_bckgsub)    
 
        self.BckgSub_run_button = Button(self.Frame_BckgSub_Button, text="Remove background", width=20)
        self.BckgSub_run_button.grid(row=1, column=0, sticky=W)  
        self.BckgSub_run_button.bind("<ButtonRelease-1>", self.BckgSub)   
       
        Label(self.Frame_BckgSub_Button, text='# pixels:').grid(row=2,column = 0)
        self.NPIX = Entry(self.Frame_BckgSub_Button) 
        self.NPIX.grid(row=2, column=1) 
        self.NPIX.insert(END,'100')  
        
        self.SelectPix_BckgSub = IntVar()
        self.Check_Live_BckgSub = Checkbutton(self.Frame_BckgSub_Button, text="Use selected pixel", variable=self.SelectPix_BckgSub).grid(column = 0,row=3)   
        
        self.SelectPix_Entry = Entry(self.Frame_BckgSub_Button) 
        self.SelectPix_Entry.grid(row=3, column=1) 
        self.SelectPix_Entry.insert(END,250)  


      
        self.Live_BckgSub = IntVar()
        self.Check_Live_BckgSub = Checkbutton(self.Frame_BckgSub_Button, text="Live processing", variable=self.Live_BckgSub).grid(column = 0,row=4)
        
        self.tree_BckgSub = CheckboxTreeview(self.Frame_BckgSub_Info)
        self.tree_BckgSub["columns"]=("Type","ExpTime","Mean","Median","std","Date")
        self.tree_BckgSub.column("#0", width=400)
        self.tree_BckgSub.column("Type", width=50)
        self.tree_BckgSub.column("ExpTime", width=50)
        self.tree_BckgSub.column("Mean", width=70)
        self.tree_BckgSub.column("Median", width=70)
        self.tree_BckgSub.column("std", width=70)
        self.tree_BckgSub.column("Date", width=200)
        
        self.tree_BckgSub.heading("#0",text="File name",anchor=W)
        self.tree_BckgSub.heading("Type",text="Type",anchor=CENTER)
        self.tree_BckgSub.heading("ExpTime", text="Exp time",anchor=CENTER)
        self.tree_BckgSub.heading("Mean", text="Mean",anchor=CENTER)
        self.tree_BckgSub.heading("Median", text="Median",anchor=CENTER)
        self.tree_BckgSub.heading("std", text="std",anchor=CENTER)
        self.tree_BckgSub.heading("Date", text="Date",anchor=W)
        
        self.tree_BckgSub.grid(column=1,row = 0)   



    def Pop_Bckg_tree(self,elem):
        hdulist = fits.open(elem)
        telescope, obsparam = CheckInstrument([elem])   
        ExpTime = 'NA'
        Type = hdulist[0].header[obsparam['obstype']]
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Bckg_Tree = self.tree_BckgSub.insert("",
                                    'end',
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))

        Children = self.tree_BckgSub.get_children()
        self.tree_BckgSub.change_state(Children[-1], "checked")


    # Function that read the fits files
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
            imgdat[np.isnan(imgdat)] = 0
            imgdat[np.isinf(imgdat)] = 0
            self.data.append(imgdat)
            
            imgdat = self.Zscale(imgdat)
            
#            print(imgdat)

            imgdat = interp.zoom(imgdat, self.zoom)

            self.images.append(Image.fromarray(imgdat))


    # Function that load the background files
    def load_file_bckgsub(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_bckgsub=[]
        for elem in self.fname:
            module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_bckgsub.append(elem)
            self.Pop_Bckg_tree(elem)
        
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


    # Correct the background
    def BckgSub(self, event):
        now = datetime.datetime.now()
        module_logger.info(now)
        
        print(self.SelectPix_BckgSub.get())
        if bool(self.SelectPix_BckgSub.get()):
            SP_BckgSub.BckgSub(self.files_bckgsub,True,'range','Bckg','bla','True',Area = [int(self.SelectPix_Entry.get())-int(self.NPIX.get()),int(self.SelectPix_Entry.get())+int(self.NPIX.get())],test = self,Live = bool(self.Live_BckgSub.get()))
        else:
            SP_BckgSub.BckgSub(self.files_bckgsub,True,'auto','Bckg','bla','True',test = self,Live = bool(self.Live_BckgSub.get()))            
#        os.system('python ' + Pipe_Path + '/SP_BckgSub.py ' + " ".join(self.files_bckgsub))
        self.now = datetime.datetime.now()
        module_logger_Bckgsub.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None


############################################################################## 
# SP_Cosmcorr
##############################################################################


    def Pop_Cosmcorr_tree(self,elem):
        hdulist = fits.open(elem)
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header[obsparam['obstype']]
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.CosmCorr_Tree = self.tree_CosmCorr.insert("",
                                    'end',
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
        os.system('python ' + Pipe_Path + '/SP_CosmCorr.py ' + " ".join(self.files_cosmcorr))
        self.now = datetime.datetime.now()
        module_logger.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +': ' + 'Cosmic correction done' )
        
        return None

############################################################################## 
# SP_Preproc
##############################################################################

    def Pop_Prep_bias_tree(self,elem):
        hdulist = fits.open(elem)
        telescope, obsparam = CheckInstrument([elem])        
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header[obsparam['obstype']]
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Bias_Tree = self.tree_preproc_MBIAS.insert("",
                                    'end',
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))

    def Pop_Prep_flat_tree(self,elem):
        hdulist = fits.open(elem)
        telescope, obsparam = CheckInstrument([elem])   
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header[obsparam['obstype']]
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Bias_Tree = self.tree_preproc_MFLAT.insert("",
                                    'end',
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))



    def Pop_Prep_toproc_tree(self,elem):
        hdulist = fits.open(elem)
        telescope, obsparam = CheckInstrument([elem])   
        ExpTime = 'NA'
        Type = hdulist[0].header[obsparam['obstype']]
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Toproc_Tree = self.tree_preproc_TOPROC.insert("",
                                    'end',
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            str("{0:.1f}".format(Mean)),
                                            str("{0:.0f}".format(Median)),
                                            str("{0:.1f}".format(Std)),
                                            Date))

        Children = self.tree_preproc_TOPROC.get_children()
        self.tree_preproc_TOPROC.change_state(Children[-1], "checked")




    def load_file_preproc(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_preproc=[]
        for elem in self.fname:
            module_logger_Preproc.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_preproc.append(elem)
            self.Pop_Prep_toproc_tree(elem)

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


    def load_file_preproc_bias(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_bias_name = []
        for elem in self.fname:
            module_logger_Preproc.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_bias_name.append(elem)
            self.Pop_Prep_bias_tree(elem)

#        self.read_all_fits(self.files_preproc)
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
#        self.index = 0;
#        im = self.images[self.index]
#        self.tkimage.paste(im)
#
#                # events
#        self.canvas.focus_set()
#        self.canvas.bind("<Key>", self.key)
#        self.canvas.bind("<Button 1>", self.left_click)
#        self.canvas.bind("<Button 2>", self.right_click)

        return None


    def load_file_preproc_flat(self, event):
        self.fname = askopenfilenames()
        self.now = datetime.datetime.now()
        self.files_flat_name = []
        for elem in self.fname:
            module_logger_Preproc.info(str(self.now.strftime("%Y-%m-%d %H:%M:%S")) +':' + ' loading file ' + str(elem).split('/')[-1])
            self.files_flat_name.append(elem)
            self.Pop_Prep_flat_tree(elem)
            
        return None


    def Preproc(self, event):
        now = datetime.datetime.now()
        module_logger_Preproc.info(now)
        
        files_to_use = []
        Children = self.tree_preproc_TOPROC.get_children()
        for elem in Children:
            if self.tree_preproc_TOPROC.tag_has("checked", elem):
                print(self.tree_preproc_TOPROC.item(elem)['text'])
                files_to_use.append(self.tree_preproc_TOPROC.item(elem)['text']) 
#                self.Pop_bias_tree(elem)
                
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
        telescope, obsparam = CheckInstrument([elem])        
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header[obsparam['obstype']]
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Flat_Tree = self.tree_flat.insert("",
                                    'end',
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
            self.Pop_flat_tree(elem)

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
        telescope, obsparam = CheckInstrument([elem])        
        ExpTime = hdulist[0].header['EXPTIME']
        Type = hdulist[0].header[obsparam['obstype']]
        Date = hdulist[0].header['DATE-OBS']
        Mean = np.mean(hdulist[0].data)
        Median = np.median(hdulist[0].data)
        Std = np.std(hdulist[0].data)
        self.Bias_Tree = self.tree_bias.insert("",
                                    'end',
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
        
        telescope, obsparam = CheckInstrument([self.fname[0]])
        
        for elem in self.fname:
            header = fits.getheader(elem) # get the header without opening data
            Type = header[obsparam['obstype']]
            ExpTime = header['EXPTIME']
            Date = header['DATE-OBS']
            Grat = header[obsparam['grating']]
            try:
                GratAng = header[obsparam['grat_ang']]
            except:
                GratAng = 'NA'
                
            Target = header[obsparam['object']]
            RotAng = header[obsparam['posangle']]
            
            try:
                CollFoc = header[obsparam['focus']]
            except:
                CollFoc = 'NA'
                
            self.Prepare_Tree = self.tree_prepare.insert("",
                                    'end',
                                    text= elem.split('/')[-1],
                                    values=(Type,
                                            str(ExpTime),
                                            Date,
                                            str(Grat),
                                            str(GratAng),
                                            str(Target),
                                            str(RotAng),
                                            str(CollFoc)))
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
                print(self.tree_prepare.item(elem)['text'])
                files_to_prepare.append(self.tree_prepare.item(elem)['text'])

        print(self.files_prepare)
        print(files_to_prepare)
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
        
        self.SelectPix_Entry.delete(0,END)
        self.SelectPix_Entry.insert(0,self.Yloc)
        
        fig, ax = plt.subplots()
        plt.plot(Line,Label = 'Spectrum')
        plt.plot(fit,Label = 'Fit')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.1, 0.9,'FWHM = ' + str(round(FWHM,1)) + ' pix',transform=ax.transAxes,fontsize=14,verticalalignment='top', bbox=props)
        ax.text(0.1, 0.8,'Trace = ' + str(region) + ' pix',transform=ax.transAxes,fontsize=14,verticalalignment='top', bbox=props)
        plt.legend()
        plt.show()
        
    def Zscale(self,imgdat):
        
        Intervals = ZScaleInterval()
        Limits = Intervals.get_limits(imgdat)

        imgdat = (np.clip(imgdat, Limits[0],Limits[1])-Limits[0])/(old_div(Limits[1]-Limits[0],256))

        return imgdat

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


    # start Tk
    app.mainloop()
 