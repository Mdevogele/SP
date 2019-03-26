#!/usr/bin/env python2
# -*- coding: utf-8 -*-

""" DIAGNOSTICS - diagnostic routines for photometry pipeline
    v1.0: 2018-04-10, mdevogele@lowell.edu
"""

from __future__ import print_function
from __future__ import division

# Spectroscopuy Pipeline 
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

### diagnostics guidelines
#
# - diagnostics.html goes into the data directory
# - all supplementary data and website go into diagroot (see _pp_init.py)
# - if there are sub-directories with data, create diagnostics.html in 
#   each directory with data; summary.html links to other directories
###

import os

import logging

import _SP_conf

from astropy.io import fits
import numpy as np


try:
    from scipy.misc import toimage # requires Pillow
    from scipy.misc import imresize # requires Pillow
    from scipy.misc import bytescale
except ImportError:
    print('Modules scipy or pillow not found. Please install with: pip '
          'install scipy pillow')
    sys.exit()
    
from astropy.visualization import (astropy_mpl_style, ZScaleInterval,
                                   ImageNormalize, LogStretch, LinearStretch)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt
    matplotlib.rcdefaults() # restore default parameters
    plt.style.use(astropy_mpl_style) # use astropy style
except ImportError:
    print('Module matplotlib not found. Please install with: pip install '
          'matplotlib')
    sys.exit()



from SP_CheckInstrument import CheckInstrument
import SP_Deveny_Toolbox as SP_Dev_Tool
import SP_Toolbox as tb


def create_website(filename, content=''):
    """
    create empty website for diagnostics output
    """

    html  = "<!DOCTYPE html PUBLIC '-//W3C//DTD HTML 4.01//EN'>\n"
    html += "<HTML>\n"
    html += "<HEAD>\n"
    html += "<TITLE>Photometry Pipeline - Diagnostics</TITLE>\n"
    html += "</HEAD>\n"
    html += "<BODY>\n"
    html += content
    html += "</BODY>\n"
    html += "</HTML>\n"

    outf = open(filename, 'w')
    outf.writelines(html)
    outf.close()

    return None


def append_website(filename, content, insert_at='</BODY>', 
                   replace_below='X?!do not replace anything!?X'):
    """ 
    append content to an existing website: 
    insert content before line that contains `insert_at`
    replace lines between `replace_below` and `insert_at` (by 
    default, nothing is replaced)
    """
    # read existing code
    existing_html = open(filename, 'r').readlines()

    # insert content into existing html
    outf = open(filename, 'w')
    delete = False
    for line in existing_html:
        if replace_below in line:
            delete = True
            continue
        if insert_at in line:
            outf.writelines(content)
            delete = False
        if delete:
            continue
        outf.writelines(line)
    outf.close()

    return None


### add File summary

def create_index(filenames, directory,
                 display=False, imagestretch='linear'):
    """
    create index.html
    diagnostic root website for one pipeline process
    """

    if display:
        print('create frame index table and frame images')
    logging.info('create frame index table and frame images')

    # obtain grating from first image file
    telescope, obsparam = CheckInstrument([filenames[0]])
    refheader = fits.open(filenames[0], ignore_missing_end=True)[0].header
    grating = [refheader[obsparam['grating']]]

    del(refheader)

    html = "<H2>data directory: %s</H2>\n" % directory

    html += ("<H1>%s/%s-band - Diagnostic Output</H1>\n" + \
               "%d frames total, see full pipeline " + \
               "<A HREF=\"%s\">log</A> for more information\n") % \
               (obsparam['telescope_instrument'], grating,
                len(filenames), 
                '.diagnostics/' + 
                _SP_conf.log_filename.split('.diagnostics/')[1])

    ### create frame information table
    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += "<TH>Idx</TH><TH>Filename</TH><TH>Date</TH>" + \
            "<TH>Objectname</TH><TH>Grating</TH>" + \
            "<TH>Airmass</TH><TH>Exptime (s)</TH>" + \
            "<TH>Type</TH>\n</TR>\n"

    # fill table and create frames
    filename = filenames
    for idx, filename in enumerate(filenames):

        ### fill table
        hdulist = fits.open(filename, ignore_missing_end=True)
        header = hdulist[0].header

        # read out image binning mode
        binning = tb.get_binning(header, obsparam)

        #framefilename = _pp_conf.diagroot + '/' + filename + '.png'
        framefilename = '.diagnostics/' + filename + '.png'

        try:
            objectname = header[obsparam['object']]
        except KeyError:
            objectname ='Unknown Target'

        html += ("<TR><TD>%d</TD><TD><A HREF=\"%s\">%s</A></TD>" + \
                 "<TD>%s</TD><TD>%s</TD>" + \
                 "<TD>%s</TD><TD>%4.2f</TD><TD>%.1f</TD>" + \
                 "<TD>%s</TD>\n</TR>\n") % \
            (idx+1, framefilename, filename, header[obsparam['date_keyword']], 
             objectname,
             header[obsparam['grating']],
             float(header[obsparam['airmass']]),
             float(header[obsparam['exptime']]),
             header[obsparam['obstype']])

   
        ### create frame image
        imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
        imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
        #resize image larger than 1000px on one side

        norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

        plt.close()
        hdulist.close()
        del(imgdat)

    html += '</TABLE>\n'
   
    create_website(_SP_conf.index_filename, html)


    ### add to summary website, if requested
#    if _pp_conf.use_diagnostics_summary:
#        add_to_summary(header[obsparam['object']], filtername, 
#                       len(filenames))

    return None

def add_BiasList(filenames,html_file,imagestretch='linear'):
    """
    add bias list to website
    """

    # update index.html
    html  = '<H2>Bias List</H2>\n'
    html += '%d biases have been detected' % \
            (len(filenames))
            
                ### create frame information table
    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH></TR>\n"
            
    for idx, filename in enumerate(filenames):
        ### fill table
        hdulist = fits.open(filename, ignore_missing_end=True)

        #framefilename = _SP_conf.diagroot + '/' + filename + '.png'
        if not os.path.isdir('diagnostics'):
            os.mkdir('diagnostics')
        if not os.path.isdir('diagnostics/bias'):
            os.mkdir('diagnostics/bias')    
        framefilename = 'diagnostics/bias' + filename + '.png'
        
        
        html += ("<TR><TD>%d</TD><TD><A HREF=\"%s\">%s</A></TD> <TH> <IMG SRC=\"%s\">  </TH></TR>\n") % \
            (idx+1, framefilename, filename,framefilename)   
        ### create frame image
        imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
        imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
        #resize image larger than 1000px on one side

        norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

        plt.close()
        hdulist.close()
        del(imgdat)             

    html += '</TABLE>\n'
    
    append_website(html_file, html, 
                   replace_below="<H2>Bias list</H2>\n")

    return None


def add_BiasSummary(filenames,MasterBias,html_file, imagestretch='linear'):
    """
    add bias processing summary to website
    """


    if not os.path.isdir('diagnostics'):
        os.mkdir('diagnostics') 
    hdulist = fits.open(MasterBias)    
    framefilename = 'diagnostics/' + MasterBias + '.png'
    
    imgdat = hdulist[0].data
    
    
    html  = '<H2>Bias Processing summary</H2>'
    html += '<p> %d bias have been processed </p>' % \
            (len(filenames))
    html += '<p> Master bias file created as <A HREF=\"%s\">%s</A> </p>' % \
            (framefilename, MasterBias)
            
    html +=  ('<p><IMG SRC=\"%s\"></p>')% \
            (framefilename)           
            
    html += '\n'
            
    html += 'Statistics of the master bias:'

    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += ("<TR><TD>Mean</TD><TD>%f</TD></TR><TR><TD>Median</TD><TD>%f</TD></TR><TR><TD>Std</TD><TD>%f</TD></TR>\n") % \
            (np.mean(imgdat),np.median(imgdat),np.std(imgdat))
    html += '</TABLE>\n'



    
    
    ### Create frame image         

    imgdat = hdulist[0].data
    # clip extreme values to prevent crash of imresize
    imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
    imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
    #resize image larger than 1000px on one side
    norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])   

    plt.figure(figsize=(5, 5))

    img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
    # remove axes
    plt.axis('off')
    img.axes.get_xaxis().set_visible(False)
    img.axes.get_yaxis().set_visible(False)

    plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

    plt.close()
    hdulist.close()
    del(imgdat)


    append_website(html_file, html, 
                   replace_below="<H2>Bias Processing summary</H2>\n")                   

    return None


def add_FlatSummary(filenames,MasterFlat,html_fileimagestretch='linear'):
    """
    add bias processing summary to website
    """


    if not os.path.isdir('diagnostics'):
        os.mkdir('diagnostics') 
    hdulist = fits.open(MasterFlat)    
    framefilename = 'diagnostics/' + MasterFlat + '.png'
    
    imgdat = hdulist[0].data
    
    
    html  = '<H2>Flat Processing summary</H2>'
    html += '<p> %d flats have been processed </p>' % \
            (len(filenames))
    html += '<p> Master flat file created as <A HREF=\"%s\">%s</A> </p>' % \
            (framefilename, MasterFlat)
            
    html +=  ('<p><IMG SRC=\"%s\"></p>')% \
            (framefilename)           
            
    html += '\n'
            
    html += 'Statistics of the master flat:'

    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += ("<TR><TD>Mean</TD><TD>%f</TD></TR><TR><TD>Median</TD><TD>%f</TD></TR><TR><TD>Std</TD><TD>%f</TD></TR>\n") % \
            (np.nanmean(imgdat),np.nanmedian(imgdat),np.nanstd(imgdat))
    html += '</TABLE>\n'



    
    
    ### Create frame image         

    imgdat = hdulist[0].data
    # clip extreme values to prevent crash of imresize
    imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
    imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
    #resize image larger than 1000px on one side
    norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])   

    plt.figure(figsize=(5, 5))

    img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
    # remove axes
    plt.axis('off')
    img.axes.get_xaxis().set_visible(False)
    img.axes.get_yaxis().set_visible(False)

    plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

    plt.close()
    hdulist.close()
    del(imgdat)


    append_website(html_file, html, 
                   replace_below="<H2>Flat Processing summary</H2>\n")                   

    return None

def add_FlatList(filenames,html_file,imagestretch='linear'):
    """
    add flat list to website
    """

    # update index.html
    html  = '<H2>Flat List</H2>\n'
    html += '%d flats have been detected' % \
            (len(filenames))
              
                ### create frame information table
    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH></TR>\n"
            
    for idx, filename in enumerate(filenames):
        ### fill table
        hdulist = fits.open(filename, ignore_missing_end=True)

        if not os.path.isdir('diagnostics'):
            os.mkdir('diagnostics') 

        #framefilename = _SP_conf.diagroot + '/' + filename + '.png'
        framefilename = 'diagnostics/' + filename + '.png'
        
        html += ("<TR><TD>%d</TD><TD><A HREF=\"%s\">%s</A></TD> <TH> <IMG SRC=\"%s\">  </TH></TR>\n") % \
            (idx+1, framefilename, filename,framefilename)    
        ### create frame image
        imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
        imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
        #resize image larger than 1000px on one side

        norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

        plt.close()
        hdulist.close()
        del(imgdat)             

    html += '</TABLE>\n'
    
    append_website(html_file, html, 
                   replace_below="<H2>Flat list</H2>\n")

    return None


def Create_Image(filename,imagestretch='linear'):
    if not os.path.isdir('diagnostics'):
        os.mkdir('diagnostics') 
    
    hdulist = fits.open(filename, ignore_missing_end=True)

    framefilename = 'diagnostics/' + filename + '.png'
    
    imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
    imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
    imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
    #resize image larger than 1000px on one side

    norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
    plt.figure(figsize=(5, 5))

    img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
    plt.axis('off')
    img.axes.get_xaxis().set_visible(False)
    img.axes.get_yaxis().set_visible(False)

    plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

    plt.close()
    hdulist.close()
    del(imgdat)                         



def add_PreProc(filenames,BiasName,FlatName,html_file,imagestretch='linear'):
    
    """
    add pre-processing log to website
    """    

    # update index.html
    html  = '<H2>Pre-processing list</H2>\n'
    html += '%d images have been processed' % \
            (len(filenames))


    html += '<p> Master bias used <A HREF=\"%s\">%s</A> </p>' % \
            (BiasName,BiasName)       




    html += '<p> Master flat used <A HREF=\"%s\">%s</A> </p>' % \
            (FlatName,FlatName)         



    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH></TR>\n"

    for idx, filename in enumerate(filenames):
        ### fill table
        hdulist = fits.open(filename + '.fits', ignore_missing_end=True)
        if not os.path.isdir('diagnostics'):
            os.mkdir('diagnostics') 

        #framefilename = _SP_conf.diagroot + '/' + filename + '.png'
        framefilename = 'diagnostics/' + filename + '.png'
        
        html += ("<TR><TD>%d</TD><TD><A HREF=\"%s\">%s</A></TD> <TH> <IMG SRC=\"%s\">  </TH></TR>\n") % \
            (idx+1, framefilename, filename,framefilename)    
        ### create frame image
        imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
        imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
        #resize image larger than 1000px on one side

        norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

        plt.close()
        hdulist.close()
        del(imgdat)  

    html += '</TABLE>\n'
    
    append_website(html_file, html, 
                   replace_below="<H2>Pre-processing list</H2>\n")      



def add_CosmCorr(filenames,html_file,imagestretch='linear'):
    
    """
    add pre-processing log to website
    """    

    # update index.html
    html  = '<H2>Cosmic-correction list</H2>\n'
    html += '%d images have been processed' % \
            (len(filenames))

    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH></TR>\n"

    for idx, filename in enumerate(filenames):
        ### fill table
        hdulist = fits.open(filename, ignore_missing_end=True)
        if not os.path.isdir('diagnostics'):
            os.mkdir('diagnostics') 

        #framefilename = _SP_conf.diagroot + '/' + filename + '.png'
        framefilename = 'diagnostics/' + filename + '.png'
        
        html += ("<TR><TD>%d</TD><TD><A HREF=\"%s\">%s</A></TD> <TH> <IMG SRC=\"%s\">  </TH></TR>\n") % \
            (idx+1, framefilename, filename,framefilename)    
        ### create frame image
        imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
        imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
        #resize image larger than 1000px on one side

        norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

        plt.close()
        hdulist.close()
        del(imgdat)  

    html += '</TABLE>\n'
    
    append_website(html_file, html, 
                   replace_below="<H2>Cosmic-correction list</H2>\n")   



def add_BckgCorr(filenames,html_file,imagestretch='linear'):
    
    """
    add pre-processing log to website
    """    

    # update index.html
    html  = '<H2>Background-correction list</H2>\n'
    html += '%d images have been processed' % \
            (len(filenames))

    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH></TR>\n"

    for idx, filename in enumerate(filenames):
        ### fill table
        hdulist = fits.open(filename, ignore_missing_end=True)
        if not os.path.isdir('diagnostics'):
            os.mkdir('diagnostics') 

        #framefilename = _SP_conf.diagroot + '/' + filename + '.png'
        framefilename = 'diagnostics/' + filename + '.png'
        
        html += ("<TR><TD>%d</TD><TD><A HREF=\"%s\">%s</A></TD> <TH> <IMG SRC=\"%s\">  </TH></TR>\n") % \
            (idx+1, framefilename, filename,framefilename)    
        ### create frame image
        imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
        imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
        #resize image larger than 1000px on one side

        norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

        plt.close()
        hdulist.close()
        del(imgdat)  

    html += '</TABLE>\n'
    
    append_website(html_file, html, 
                   replace_below="<H2>Background-correction list</H2>\n")         
    
    
    
def add_Extract(filenames_im,filenames_spec,html_file,imagestretch='linear'):
    
    """
    add pre-processing log to website
    """    

    # update index.html
    html  = '<H2>Extract list</H2>\n'
    html += '%d images have been processed' % \
            (len(filenames_im))

    html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
    html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH></TR>\n"

    for idx, filename in enumerate(filenames_im):
        ### fill table
        spec = np.loadtxt(filenames_spec[idx])
        plt.figure()
        plt.plot(spec)
        plt.xlabel('pixels')
        plt.ylabel('Relative reflectance')
        specfilename = 'diagnostics/' + filenames_spec[idx] + '.png'
        plt.savefig(specfilename)
        plt.close()
        hdulist = fits.open(filename, ignore_missing_end=True)
        if not os.path.isdir('diagnostics'):
            os.mkdir('diagnostics') 

        #framefilename = _SP_conf.diagroot + '/' + filename + '.png'
        framefilename = 'diagnostics/' + filename + '.png'
        
        html += ("<TR><TD>%d</TD><TD><A HREF=\"%s\">%s</A></TD> <TH> <IMG SRC=\"%s\">  </TH> <TH> <IMG SRC=\"%s\">  </TH></TR>\n") % \
            (idx+1, framefilename, filename,framefilename,specfilename)    
        ### create frame image
        imgdat = hdulist[0].data
        # clip extreme values to prevent crash of imresize
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                            np.percentile(imgdat, 99))
        imgdat = imresize(imgdat, 
                          min(1., 1000./np.max(imgdat.shape)), 
                          interp='nearest')
        #resize image larger than 1000px on one side

        norm = ImageNormalize(imgdat, interval=ZScaleInterval(),
                      stretch={'linear': LinearStretch(),
                               'log': LogStretch()}[imagestretch])
        
        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        plt.savefig(framefilename, format='png', bbox_inches='tight', 
                    pad_inches=0, dpi=200)

        plt.close()
        hdulist.close()
        del(imgdat)  

    html += '</TABLE>\n'
    
    append_website(html_file, html, 
                   replace_below="<H2>Extract list</H2>\n")            
    
    