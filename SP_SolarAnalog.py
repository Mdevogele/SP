#!/usr/bin/env python2
# -*- coding: utf-8 -*-


""" SP_SA - Database containing the usual Solar Analog used by MANOS
    v1.0: 2018-03-15, Mdevogele@lowell.edu
"""


import cPickle as pickle
from astroquery.simbad import Simbad



List_of_SA = ['SA 93-101', 'HD  28099', 'SA 98-978', 'SA 102-1081', 'SA 105-56', 'SA 107-684',
              'SA 107-998', 'SA 110-361', 'SA 112-1333', 'SA 113-276', 'SA 115-271']

List_ID_SA = ['HD  11532', 'HD  28099', 'HD 292561', 'BD+00  2717','BD-00  2719','HD 139287',
              'BD+00  3383', 'TYC  447-508-1','BD-00  4074', 'BD-00  4251B','BD-00  4557']

def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    for elem in List_of_SA:
        result_table = Simbad.query_object(elem)
        save_object(result_table, '/Users/maximedevogele/Documents/PythonPackages/SP/SA/' + elem + '.Sim')
        print(result_table['MAIN_ID'][0])
    




#SA= {'SA 93-101' : {'RA': {'D' : 01, 'M' : 53, 'S' : 18.0 }, 'DEC' :{'D' : 00,  'M' : 22, 'S' : 25 }}, 
#     'Hya 64'    : {'RA': {'D' : 04, 'M' : 26, 'S' : 40.1 }, 'DEC' :{'D' : 16,  'M' : 44, 'S' : 49 }},
#     'SA 98-978' : {'RA': {'D' : 06, 'M' : 51, 'S' : 34.0 }, 'DEC' :{'D' : -00, 'M' : 11, 'S' : 33 }},
#     'SA102-1081': {'RA': {'D' : 10, 'M' : 57, 'S' : 04.4 }, 'DEC' :{'D' : -00, 'M' : 13, 'S' : 12 }},
#     'SA105-56'  : {'RA': {'D' : 13, 'M' : 38, 'S' : 42.1 }, 'DEC' :{'D' : -01, 'M' : 14, 'S' : 14 }},
#     'SA107-684' : {'RA': {'D' : 15, 'M' : 37, 'S' : 40.1 }, 'DEC' :{'D' : 16, 'M' : 44, 'S' : 49 }}
#     }


#qd = {'Objects':{'use_me':1,
#       'Instrument':'GMOS-S','CcdBin':'','RoI':'Full',
#       'Disperser':'R400+_%','CentWave':700.0,'AperMask':'2.0arcsec',
#       'Object':'*',
#       'DateObs':''}}