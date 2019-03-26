#!/usr/bin/env python

""" SP_Logs - Display or save a logs files
    v1.0: 2018-06-13, mdevogele@lowell.edu
"""

import argparse, shlex

import SP_Toolbox as SP
import SP_Deveny_Toolbox as SP_Dev_Tool
import numpy as np
#import SP_diagnostics as diag

from SP_CheckInstrument import CheckInstrument

from astropy.io import fits

import SP_diagnostics as diag

import _SP_conf

def Log(filename,Outfile,verbode,Diagnostic):
    
    


    telescope, obsparam = CheckInstrument([filename[0]])
    print(telescope)
    if Diagnostic:
        diag.create_website('Log.html')
        html  = '<H2>Log</H2>\n'
        html += "<P><TABLE BORDER=\"1\">\n<TR>\n"
        if telescope == 'GMOSS' or telescope == 'GMOSN':
            html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH></TR>\n"
        else:
            html += "<TH>Instrument</TH><TH>Index</TH> <TH>Obs type</TH> <TH>Object</TH> <TH>Exp time [sec]</TH>  <TH>RA</TH> <TH>DEC</TH> <TH>Airmass</TH> <TH>Time</TH>  <TH>Images</TH></TR>\n"

    for idx, elem in enumerate(filename): 
        Splitted = elem.split('_')
        hdulist = fits.open(elem)
        
        if Diagnostic:
            diag.Create_Image(elem)
            framefilename = 'diagnostics/' + elem + '.png'

        
        INST = obsparam['telescope_instrument']
        OBSTYPE = hdulist[0].header[obsparam['obstype']]
        if telescope == 'SOAR':
            if elem.split('_')[3] == 'ACQ':
                OBSTYPE = 'ACQ'
        OBJECT = hdulist[0].header[obsparam['object']]
        if telescope == 'SOAR':
            if OBSTYPE == 'COMP':
                OBJECT = elem.split('_')[4]     
        EXPTIME = hdulist[0].header[obsparam['exptime']]
        ALPHA = hdulist[0].header[obsparam['ra']]
        DEC = hdulist[0].header[obsparam['dec']]
        AM = hdulist[0].header[obsparam['airmass']]
        TIME = hdulist[0].header[obsparam['date_keyword']]
        IDX = Splitted[1]
        GRAT = hdulist[0].header[obsparam['grating']]
        
        if telescope == 'GMOSS' or telescope == 'GMOSN':
            WAV = hdulist[0].header[obsparam['tilt']] 
            print('%13s %6s %10s %6s %15s %7s %3.3f %3.3f \t %5s %25s %15s' % (INST, IDX, OBSTYPE, str(int(WAV)), OBJECT, str(int(EXPTIME)), ALPHA, DEC,str(AM), TIME, GRAT ))
            if Diagnostic:
                html += "<TH>Idx</TH><TH>Filename</TH> <TH>Images</TH> <TH>Images</TH></TR>\n"
        else:
            print('%13s %6s %10s %15s %7s %12s %12s %5s %25s' % (INST, IDX, OBSTYPE, OBJECT, str(int(EXPTIME)), str(ALPHA), str(DEC),str(AM), TIME))
            if Diagnostic:
                html += ("<TR><TD>%s</TD><TD>%s</TD> <TD>%s</TD>  <TD>%s</TD> <TD>%s</TD> <TD>%s</TD> <TD>%s</TD> <TD>%s</TD> <TD>%s</TD> <TH> <IMG SRC=\"%s\">  </TH></TR>\n") % \
            (INST, IDX, OBSTYPE,OBJECT,str(int(EXPTIME)),str(ALPHA),str(DEC), str(AM), TIME,framefilename) 
                
 #       print(INST + '\t' + OBSTYPE + '\t' + OBJECT + '\t'+ str(EXPTIME) + '\t' +  str(ALPHA) + '\t' + str(DEC) + '\t' + str(AM)  + '\t' + TIME )
        
 
    html += '</TABLE>\n'
    diag.append_website('Log.html',html,replace_below="<H2>Log</H2>\n")


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Create a log file')

    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-o',
                        help='Name of the log file',
                        default='Logs.txt')    

    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True) 


    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
    Verbose = args.v
    Outfile = args.o
    filenames = args.images 
    Diagnostic = args.d
    
    
    # call run_the_pipeline only on filenames
    Log(filenames,Outfile,Verbose,Diagnostic)
    pass
