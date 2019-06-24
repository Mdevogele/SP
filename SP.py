#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:31:15 2019

@author: maximedevogele
"""

from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askopenfilenames
from tkinter.messagebox import showerror

import argparse, shlex
import numpy as np

import time
from datetime import date
from datetime import timedelta

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os

import _SP_conf
Pipe_Path = _SP_conf.rootpath


class MyFrame(Frame):
    def __init__(self):
        Frame.__init__(self)
        self.master.title("Spectroscopy pipeline")
        self.master.rowconfigure(5, weight=1)
        self.master.columnconfigure(5, weight=1)
        self.grid(sticky=W+E+N+S)

        self.LF_button = Button(self, text="Load files", command=self.load_file, width=10)
        self.LF_button.grid(row=1, column=0, sticky=W)
                
        
        self.Plot_button = Button(self, text="Plot Spectrum", command=self.Plot_Spectrum, width=10)
        self.Plot_button.grid(row=2, column=0, sticky=W)   
        
        self.Fits_button = Button(self, text="Display fits", command=self.Display_fits, width=10)
        self.Fits_button.grid(row=3, column=0, sticky=W)   

    def load_file(self):
        self.fname = askopenfilenames()
        return None

    def Plot_Spectrum(self):
        plt.figure()
        for elem in self.fname:
            Spec = np.loadtxt(elem).transpose()
            plt.plot(Spec[0],Spec[1],'.')
        plt.show()
        
        return None
    
    def Display_fits(self):
    
        os.system('python ' + Pipe_Path + '/SP_DisplayFits.py ' + " ".join(self.fname))
        
        return None
        

if __name__ == "__main__":
    MyFrame().mainloop()







#
#def OpenFile(master):
#    
#    master.filename = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
#
#
#if __name__ == '__main__':
#
#
#    master = Tk(screenName='Spectroscopy pipeline',  baseName=None,  className='Tk',  useTk=1)
#    
#    button = Button(master, text='Get objects', width=25, command=lambda: OpenFile()).grid(row=1,columnspan=1,rowspan=1)
#    
#    master.filename = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
#    print(master.filename)
#    
#    master.mainloop()
#
#    pass