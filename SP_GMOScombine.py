#!/usr/bin/env python

""" SP_GMOSCombine - Combine final GMOS spectra
    v1.0: 2018-06-19, mdevogele@lowell.edu
"""

import argparse, shlex
import numpy as np
import SP_Toolbox as SP
import matplotlib.pyplot as plt
import math


def GMOScombine(filenames,Date,Target,Facility):
    
    if Target == False:
        Target = filenames[0].split('_')[0]

    SpecT = []
    SpecA = []
    for elem in filenames:
        
        with open(elem,'r') as f:
            SpecA = repr(f.read().splitlines())
        SpecA = SpecA.replace('"','').replace('[','').replace(']','').replace("'",'').replace("\\",'').replace('t',' ').replace(',','').split()    
    
        SpecA = np.array(SpecA).astype(float)
        print(len(SpecA))
        SpecA = SpecA.reshape(len(SpecA)/2,2)
        
     #   SpecA[SpecA[:,0]>1,1] = np.nan
        SpecA[SpecA[:,1]<0.001,1] = np.nan

        z = np.ma.polyfit(SpecA[:,0],np.nan_to_num(SpecA[:,1]) , 5)
        p = np.poly1d(z)
        SSP = SpecA[:,1] - p(SpecA[:,0])
        
        median = np.nanmedian(SSP)
        std = np.nanstd(SSP)        
        PP = list(SSP[SSP > median+3*std])
        PG = list(SSP[SSP < median-3*std])
        DO = True    
        while DO: 
            DO = False
            if PP or PG:
                median = np.nanmedian(SSP)
                std = np.nanstd(SSP)
                IDX = SSP > median+3*std
                SpecA[IDX,1] = np.nan
                SSP[SSP > median+3*std] = np.nan
                
                median = np.nanmedian(SSP)
                std = np.nanstd(SSP)
                IDX = SSP < median-3*std
                SpecA[IDX,1] = np.nan
                SSP[SSP < median-3*std] = np.nan
                
                median = np.nanmedian(SSP)
                std = np.nanstd(SSP) 
                PP = list(SSP[SSP > median+3*std])
                PG = list(SSP[SSP < median-3*std])
                DO = True

        Wave = np.linspace(0.45,1,200)
        Dw = (1-0.45)/200
    
        Wave_out = []
        Spec_out = []
        Error = []
    
        Wav = SpecA[:,0]
        Spec = SpecA[:,1]
        
        Spec = Spec/np.nanmedian(Spec[800:850])
        for k in Wave:
            Wave_out.append(k+Dw/2)
            Spec_Tamp = []
        
            condition1 = Wav > k
            condition2 = Wav < k + Dw
            condition = condition1*condition2
            Spec_Tamp.append(np.nanmedian(np.extract(condition, Spec)))
            Spec_out.append(Spec_Tamp)
            Error.append(np.nanstd(np.extract(condition, Spec)*~np.isinf(np.extract(condition, Spec)))/np.sqrt(len(np.extract(condition, Spec))-1))        

        SpecT.append(Spec_out)
        
        
    SpecT = np.array(SpecT)
    SpecTT = np.nanmedian(np.array(SpecT).astype(float),axis=0)
    
    idx = (Wave>0.53)*(Wave<0.57)
    SpecTT = SpecTT/np.nanmedian(SpecTT[idx])
    

#    if np.isnan(SpecTT[-1]):
#        SpecTT[-1] = SpecTT[-25]
    

#    Date = '1704012'
    
    Tax = SP.Get_Taxonomy(Wave,SpecTT)
    f = open(Target + '_' + Date + '_' + Facility + '_Tax.dat','w')
    for ele in Tax:
        f.write(ele[0] + ',' + str(ele[1]) + '\n')
    f.close()
    plt.figure()
    SP.Plot_Taxonomy(Wave,SpecTT,Target,Date,Facility)
     #       plt.savefig('spec.jpg',dpi=1200)
    plt.savefig(Target + '_' + Date + '_' + Facility + '_spec.jpg',dpi=1200)
        
    f = open(Target + '.spec','w')
    for wavel, refl in zip(Wave,SpecTT):
        if not math.isnan(refl):
            f.write(str(wavel)+ '\t' + str(refl) + '\n')
    f.close()
    
    
    
    
    
#    SP.Plot_Taxonomy(Wave,SpecTT,Target,Date,'GMOS')
#    plt.savefig(Target + '_' + Date + '_' + 'GMOS' + '_spec.jpg',dpi=1200)
    
#    f = open(OutFile,'w')
#    for item in SpecTT:
#      f.write("%s\n" % item)
#    f.close()
    
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract spectrum')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
#    parser.add_argument('-target', help='primary targetname override',
#                        default=None)
#    parser.add_argument('-m', help='add flats information to diagnostic.htlm file',
#                        choices=['auto','range'])
#    parser.add_argument('-s', help='If there is several series of flat \n || Options: none: Only one serie \n || index: split according to the index of the files \n || target: split according to the target name in fits headers \n || pointing: split according to telescope pointing according to fits headers',
#                        choices=['none','index','target','pointing'],
#                        default = 'None')
#    parser.add_argument('-b',
#                        help='Name of the master bias to use \n || Can use None if no bias are available',
#                        default='MasterBias.fits') 
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    parser.add_argument('-d',
                        help='Date',
                        default = '20171019')
    
    parser.add_argument('-f',
                        help='Facility',
                        default = 'Gemini')
    
    parser.add_argument('-t',
                        help='Target name',
                        default = False)

    args = parser.parse_args()
    Date = args.d
    Target = args.t
    Facility = args.f
#    Method = args.m
#    Series = args.s
#    MasterBias = args.b
    filenames = args.images  
    
    print(filenames)

    
    
    # call run_the_pipeline only on filenames
    GMOScombine(filenames,Date,Target,Facility)
    pass
    