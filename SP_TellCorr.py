#!/usr/bin/env python

""" SP_TellCorr - Perform the wavelength calibration
    v1.0: 2018-04-20, mdevogele@lowell.edu
"""

import SP_Toolbox as SP
import matplotlib.pyplot as plt
import math
import argparse, shlex
import numpy as np

from SP_CheckInstrument import CheckInstrument

def TellCorr(filenames,Std,Verbose,Target,Date,Instrument):

    DetecFlags = {'Instrument':Instrument}  

    for idx,elem in enumerate(filenames):
    
        if Target == False:
            Target = elem.split('.')[0]
   


        Spec = np.loadtxt(elem).transpose()     
        
        if np.size(Spec,0) != 2:
            print('ERROR: Wrong asteroid spectrum file')
            print('ERROR: It is not a two columns text file')
            print('ERROR: SP_TellCorr takes as argument the output of SP_WavCal')
            return
        
        Wav = Spec[0]
        SpecA = Spec[1]
        
        
        idx = (Wav>0.505)*(Wav<0.605)
        index = np.where(idx)
        SpecA = SpecA/np.nanmedian(SpecA[index])
    
        
        Spec = np.loadtxt(Std).transpose()     

        if np.size(Spec,0) != 2:
            print('ERROR: Wrong solar analog spectrum file')
            print('ERROR: It is not a two columns text file')
            print('ERROR: SP_TellCorr takes as argument the output of SP_WavCal')
            return

        WavSA = Spec[0]
        SpecSA = Spec[1]
    
        SpecSA = SpecSA/np.nanmedian(SpecSA[index])
    
        Wavel = [Wav, Wav]
        Spectre = [SpecA, SpecSA]
    
        SpecN, WavN = SP.Shift_Spec(Spectre,Wavel,**DetecFlags)
        
        SpecA = SpecN[0]
        SpecSA = SpecN[1]
    
    
        Spec = SpecA/SpecSA
        WavNew = np.array(WavN[0])
        
        
        if Instrument == 'Deveny' or Instrument == 'DEVENY':
            Tax = SP.Get_Taxonomy(WavNew[10:-10],Spec[10:-10])
            f = open(Target + '_' + Date + '_' + 'DCT' + '_Tax.dat','w')
            for ele in Tax:
                f.write(ele[0] + ',' + str(ele[1]) + '\n')
            f.close()
            plt.figure()
            SP.Plot_Taxonomy(WavNew[10:-10],Spec[10:-10],Target,'170514','DCT')
     #       plt.savefig('spec.jpg',dpi=1200)
            plt.savefig(Target + '_' + Date + '_' + 'DCT' + '_spec.jpg',dpi=1200)
        
            f = open(Target + '.spec','w')
            for wavel, refl in zip(WavNew[10:-10],Spec[10:-10]):
                if not math.isnan(refl):
                    f.write(str(wavel)+ '\t' + str(refl) + '\n')
            f.close()
         
        if 'GMOS' in Instrument:
            f = open(Target + '.spec','w')
            for wavel, refl in zip(WavNew[10:-10],Spec[10:-10]):
                if not math.isnan(refl):
                    f.write(str(wavel)+ '\t' + str(refl) + '\n')
            
        


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract spectrum')

    parser.add_argument('-d', help='',
                        default = '20171019')
    parser.add_argument('-t',
                        help='Target name',
                        default = False)
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-s',
                        help='Standard star spectrum to use',
                        default = 'SpecOut.spec')    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')
    
    parser.add_argument('-i',
                        help='Instrument')  

    args = parser.parse_args()

    Date = args.d
    Target = args.t
    Verbose = args.v
    Std = args.s
    Instrument = args.i
    filenames = args.images  

    
    TellCorr(filenames,Std,Verbose,Target,Date,Instrument)
    pass