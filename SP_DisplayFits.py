#!/usr/bin/env python

""" SP_DisplayFits - Display fits files 
    v1.0: 2010-06-24, mdevogele@lowell.edu
"""

from __future__ import print_function
from __future__ import division

from past.utils import old_div
import os, sys
import numpy
import warnings
from tkinter import *
from PIL import Image
from PIL import ImageTk
from PIL import ImageDraw
import argparse
from astropy.io import fits
from scipy.ndimage import interpolation as interp
from scipy.interpolate import InterpolatedUnivariateSpline

from astropy.stats import sigma_clipped_stats


from photutils.datasets import make_4gaussians_image
from photutils import centroid_com, centroid_1dg, centroid_2dg

import operator

import time

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


import SP_Toolbox as tb


# self.image = number of the displayed image




class Clicker(object):

    def __init__(self, master, zoom, filelist):
        self.top = master
        self.files = filelist
        self.zoom = zoom
        self.target_index = [None for i in range(len(self.files))]
        self.interp_index = [None for i in range(len(self.files))]
        self.index = 0
        self.images = []
        self.ldac   = []
        self.mjd    = []
        
        self.redcircle = []

        self.JD = []

        ### load image data
        print('please wait, loading images...', end=' ')
        sys.stdout.flush()

        self.read_all_fits(self.files)

        print('done!')

        # create title bar
        self.title = Label(text='%s (%d/%d)' %
                           (self.images[0],
                            self.index+1, len(self.files)))
        self.title.pack()


        # select first image
        self.tkimage = ImageTk.PhotoImage(self.images[0], palette=256)
        self.canvas = Canvas(master, height=self.tkimage.height(), width=
                             self.tkimage.width())
        self.image = self.canvas.create_image(0, 0, anchor='nw',
                                              image=self.tkimage)    

        # frame counter variable
        self.evar = IntVar()
        self.evar.set(1)

        self.canvas.pack(side='top', expand=1)

        # display image
        self.nextframe()

        # events
        self.canvas.focus_set()
        self.canvas.bind("<Key>", self.key)
        self.canvas.bind("<Button 1>", self.left_click)
        self.canvas.bind("<Button 2>", self.right_click)


    def left_click(self,event):
        x, y = old_div(event.x,self.zoom), old_div(event.y,self.zoom) 
        self.x = x
        self.y = y

    def right_click(self, event):
        """ select source """
        x, y = old_div(event.x,self.zoom), old_div(event.y,self.zoom)
#        
#        
#        
        self.Select_Target(self.data[self.index],x,y,30)


    def Select_Target(self,image,x,y,box):
        
        X = []
        Y = []
        
        D = image[int(y)-box:int(y)+box,int(x)-box:int(x)+box]
        
        Line = numpy.median(D,axis=1)
        
        
        xs = range(len(Line))
        ys = Line
        
        max_index, max_value = max(enumerate(Line), key=operator.itemgetter(1))
        
        p0 = [0,max_value,2.24,1.41,max_index] 
        Res = []

        coeff, var_matrix = curve_fit(lambda x, S0,S1,a, b,x0: tb.MOFFAT(x, S0,S1,a, b,x0),xs,ys,p0,maxfev = 100000)
        fit = tb.MOFFAT(xs, *coeff)
        print(fit)
        
        FWHM = 2*coeff[2]*numpy.sqrt(2**(1/coeff[3])-1)
        
        region = sum(Line>(coeff[0]+50))/2
        
        
        fig, ax = plt.subplots()
        plt.plot(Line,Label = 'Spectrum')
        plt.plot(fit,Label = 'Fit')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.1, 0.9,'FWHM = ' + str(round(FWHM,1)) + ' pix',transform=ax.transAxes,fontsize=14,verticalalignment='top', bbox=props)
        ax.text(0.1, 0.8,'Trace = ' + str(region) + ' pix',transform=ax.transAxes,fontsize=14,verticalalignment='top', bbox=props)
        plt.legend()
        plt.show()


    def key(self, event):
        """ keyboard events """
        if event.char == 'a':
            # previous frame
            self.nextframe(-1)
        elif event.char == 'd':
            # next frame
            self.nextframe(1)
        elif event.char == 'q':
            # quit
            self.top.quit()
        elif event.char == '+':
            # zoom in
            if self.zoom < 4.:
                self.zoom *= 2
            self.nextframe(0)
        elif event.char == '-':
            # zoom out
            if self.zoom > 0.25:
                self.zoom /= 2
            self.nextframe(0)
        elif event.char == 'z':
            self.Select_Target(self.data[self.index],self.x,self.y,30)
        
        elif event.char == 'x':
            self.Show_Targets()
        elif event.char == 'p':
            self.Photometrie()    

#    def right_click(self, event):
#        """ next frame """
#        self.nextframe(1)

    def read_all_fits(self, filenames, zoom=0.5):
        """ read in all image data, scale images """
        self.data=[]
        self.Tx = []
        self.Ty = []
        for idx, filename in enumerate(filenames):
            if idx > 0:
                print('\b\b\b\b%3d' % (idx+1), end=' ')
            else:
                print('%3d' % (idx+1), end=' ')
            sys.stdout.flush()

            ## read image data
            hdulist = fits.open(filename, ignore_missing_end=True)
            imgdat = hdulist[0].data
            self.data.append(imgdat)

            median = numpy.median(imgdat[int(imgdat.shape[1]*0.25):
                                         int(imgdat.shape[1]*0.75),
                                         int(imgdat.shape[0]*0.25):
                                         int(imgdat.shape[0]*0.75)])
            std    = numpy.std(imgdat[int(imgdat.shape[1]*0.25):
                                      int(imgdat.shape[1]*0.75),
                                      int(imgdat.shape[0]*0.25):
                                      int(imgdat.shape[0]*0.75)])

            imgdat = old_div(numpy.clip(imgdat, median-0.5*std,
                                median+0.5*std),(old_div(std,256)))
            imgdat = imgdat - numpy.min(imgdat)

            imgdat = interp.zoom(imgdat, self.zoom)

            self.images.append(Image.fromarray(imgdat))

    def nextframe(self,i=1, imgnum=-1):
        """ display frame using iterator i"""

        if imgnum == -1:
            self.index += i
        else:
            self.index = imgnum - 1
        if self.index >= len(self.files):
            self.index = 0
        elif self.index < 0:
            self.index = len(self.files) - 1
        filename = self.files[self.index]
        if not os.path.exists(filename):
            print("Unable to find %s" % filename)
            self.top.quit()
        self.evar.set(self.index+1)

        self.title.configure(text='%s (%d/%d)' %
                           (os.path.basename(filename),
                            self.index+1, len(self.files)))

        im = self.images[self.index]

        self.tkimage.paste(im)
        
        #delete previous circles
        for elem in self.redcircle:
            self.canvas.delete(elem)
            



if __name__ == '__main__':
    
    
    # define command line arguments
    parser = argparse.ArgumentParser(description='manual target identification')
    parser.add_argument('-zoom', help='image zoom factor', default=0.5)
    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    zoom = float(args.zoom)
    filenames = args.images

    root = Tk()
    app = Clicker(root, zoom, filenames)
    root.mainloop()

    pass


