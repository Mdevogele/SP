#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:31:32 2019

@author: maximedevogele
"""

###### Function that does not appear to be necessary anymore. Keeping them here to get them if they actually were.


def Get_AsteroidList():
    
    filenames = _SP_conf.filenames
    
    Asteroid_List = []
    for idx, elem in enumerate(filenames):
        if elem.split('_')[3] == 'Asteroid':
            Asteroid_List.append(elem)
 
    _SP_conf.Asteroid_filenames = Asteroid_List


def Get_FlatList():
    
    filenames = _SP_conf.filenames
    
    Flat_List = []
    for idx, elem in enumerate(filenames):
        if elem.split('_')[3] == 'FLAT':
            Flat_List.append(elem)
 
    _SP_conf.Flat_filenames = Flat_List

def Get_BiasList():
    
    filenames = _SP_conf.filenames
    
    Bias_List = []
    for idx, elem in enumerate(filenames):
        if elem.split('_')[3] == 'BIAS':
            Bias_List.append(elem)
 
    _SP_conf.Bias_filenames = Bias_List
    
    
def Write_Log(string,Proc_Dir):
    now = datetime.datetime.now()
    f = open('./' + Proc_Dir + '/ProcLogFile','a')
    
    ToWrite = now.isoformat() + ':'
    
    ToWrite = ToWrite + '\t' + string + '\n'
    
    f.write(ToWrite)
    f.close()
    
    print(ToWrite)
    
def Wav_Cal(Arc, Disp = 150,**kw):
    
    ''' need to improve the function to work with any bin value and any R values ''' 
    
    if kw.has_key('Instrument'):
        Instrument = kw['Instrument']
    else:
        Instrument = 'GMOS'     


    if Instrument == 'Deveny':
        print('Instrument = Deveny')
        Master = '/Users/maximedevogele/Documents/PythonPackages/SP/Deveny_Arcs_Master'
        with open(Master) as f:
            Master_Arcs = f.read().splitlines()
        Master_Arcs = np.array(Master_Arcs).astype(float)
        Master_Arcs = Master_Arcs/39
        
        hdulist = fits.open(Arc)
        data = hdulist[0].data
        hdulist.close()
        
        
        data = data - np.median(data);
        data = data / np.mean(data);
    
        Lines = np.median(data[250:260,:],axis=0)
        LL = [Lines, Lines,Lines, Lines,Lines, Lines,Lines, Lines]
        LL2 = [Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs,Master_Arcs, Master_Arcs]
        
        result = ird.similarity(np.array(LL2), np.array(LL), numiter=3)
        
        x = np.array(range(2148))
        Wav = (x+result['tvec'][1])*-4.277+11780
        
    if Instrument == 'GMOS':
        if Disp == 150:
            Master = '/Users/maximedevogele/Documents/Gemini/CuAr/Gemini_MasterArc_R150.fits'
            WavCalFile = '/Users/maximedevogele/Documents/Gemini/CuAr/Wav_Cal_MasterArc_R150'
        if Disp == 400:
            Master = '/Users/maximedevogele/Documents/Gemini/CuAr/Gemini_MasterArc_R400.fits'
            WavCalFile = '/Users/maximedevogele/Documents/Gemini/CuAr/Wav_Cal_MasterArc_R400'

        hdulist = fits.open(Master)
        Master_Arcs = hdulist[0].data
        hdulist.close()
        
        WavCal = []
        with open(WavCalFile) as f:
             WavCal = [x.split() for x in f.readlines()]
        WavCal = np.array(WavCal).astype(float)     
        
        hdulist = fits.open(Arc)
        data = hdulist[0].data
        hdulist.close()
        
        result = ird.similarity(Master_Arcs[990:1000,800:950], data[990:1000,800:950], numiter=3)
        WavCalSci = copy.deepcopy(WavCal)
        
        z = np.polyfit((WavCalSci[:,0]-np.mean(WavCalSci[:,0]))/np.std(WavCalSci[:,0]),WavCalSci[:,1], 4)
        p = np.poly1d(z)
        
        Xaxis = np.linspace(0,3132,3132)
        Wav = p((Xaxis-np.mean(WavCalSci[:,0]))/np.std(WavCalSci[:,0]))

    return Wav


def Auto_Extract_Spec(bla):
    
    print(bla)
    hdulist = fits.open(bla)
    data = hdulist[1].data
    dataC = Correct_Amplifier(data)
    dataCC = Correct_Boxes(dataC)
#    dataCC = dataC
    Center = Detect_Spectra(dataCC)
    Start = (685,Center)
    Trace, bkg, MASK1 = Fit_Trace(dataCC,Start,Range = 400, SClip = True)
    Spec1 = Extract_Spectrum(dataCC,Trace,bkg,FWHM=10,Mask = MASK1)
    Wave1 = Extract_Wave(bla,Trace)
    Spec1S = Sig_Clip_Spec(Spec1,n = 7, sig = 4)
    Spec1N = Normalize_Spectrum(Spec1,Wave1,MASK1*Spec1S,wavelength=7000)
    
    return Wave1,Spec1N


def Get_Wavtrans(name):
    
    hdulist = fits.open(name)
    WavTransFile = hdulist[1].header['WAVTRAN']
    
    F = open('./database/fc' + WavTransFile ,'r')
    
    for i in range(14):
        F.readline()
    
    CoeffTamp = []
    Coeff = []
    for i in range(16):
        CoeffTamp = F.readline()
        Coeff.append(float(CoeffTamp))
    
    Coeff = np.array(Coeff)
    
    XS,YS = np.meshgrid(np.linspace(-1,1,1044),np.linspace(1,-1,1566))
    
    
    Wave = np.polynomial.chebyshev.chebval2d(YS,XS,Coeff.reshape((4,4),order='F'))
    
    return Wave


def Correct_Boxes(data):
    
    B1x1 = 904
    B1x2 = 970
    
    B2x1 = 973
    B2x2 = 1028
    
    
    Cons = Detect_Boxes(data)
    
    for List in Cons:
        print(List)
    
        XSize = len(List) + 1
        YSize1 = B1x2-B1x1
        YSize2 = B2x2-B2x1
        
        
        X1 = []
        X2 = []
        Y1 = [] 
        Y2 = []
        
        X1, Y1 = np.meshgrid(range(XSize), range(YSize1), sparse=False, indexing='ij')
        X2, Y2 = np.meshgrid(range(XSize), range(YSize2), sparse=False, indexing='ij')
    
        X1 = X1.flatten()
        Y1 = Y1.flatten()
    
        X2 = X2.flatten()
        Y2 = Y2.flatten()
    
        Z1 = []
        Z2 = []
        
        Z1 = data[List[0]-1:List[-1]+1,B1x1:B1x2]
        Z2 = data[List[0]-1:List[-1]+1,B2x1:B2x2]
        
        ZC1 = Sigma_Clip(Z1, sig = 2.5)
        ZC2 = Sigma_Clip(Z2, sig = 2.5)
        
#        ZCC1 = Sigma_Clip(Z1, sig = 3)
#       ZCC2 = Sigma_Clip(Z2, sig = 3)
        
    
        A1 = np.array([X1*0+1, X1, Y1, X1*Y1]).T
        A2 = np.array([X2*0+1, X2, Y2, X2*Y2]).T
        
        B1 = Z1.flatten()
        B2 = Z2.flatten()
        
        C1 = ZC1.flatten()
        C2 = ZC2.flatten()
        
        Index1 = [i for i, x in enumerate(C1) if x]
        Index2 = [i for i, x in enumerate(C2) if x]
    
    
        coeff1, r, rank, s = np.linalg.lstsq(A1[Index1], B1[Index1])
        coeff2, r, rank, s = np.linalg.lstsq(A2[Index2], B2[Index2])
    
        ZZ1 = []
        ZZ1 = np.empty((XSize,YSize1), dtype=float)
        for i in range(XSize):
            for j in range(YSize1):
                ZZ1[i,j] = coeff1[0]+coeff1[2]*j + coeff1[1]*i+coeff1[3]*i*j
                
        data[List[0]-1:List[-1]+1,B1x1:B1x2] -= ZZ1    
    
        ZZ2 = [] 
        ZZ2 = np.empty((XSize,YSize2), dtype=float)
        for i in range(XSize):
            for j in range(YSize2):
                ZZ2[i,j] = coeff2[0]+coeff2[2]*j + coeff2[1]*i+coeff2[3]*i*j
    
        data[List[0]-1:List[-1]+1,B2x1:B2x2] -= ZZ2
    
    return data

def Detect_Boxes(data):

    B = range(120)
    B = np.array(B)
    B[0]=0
    B[-1]=0
    B[1:-1] = 1 
    
    TT = []

    for i in range(1044):
        c = data[i,902:1022]
        SIG = Sigma_Clip(c)
        
        M1 = data[i,800:850]
        M2 = data[i,1100:1140]
        
#        MS1 = Sigma_Clip(M1)
#        MS2 = Sigma_Clip(M2)
        
        M = (np.median(data[i,800:850])+np.median(data[i,1100:1140]))/2
        c = c/(abs(M)+1)
        TT.append( abs( np.median(np.convolve(c[SIG], B)) ) )
        
    PREUP = np.median(TT)+0.5*np.std(TT)
    PRESUP = TT > PREUP
        
    
    Indices = []
    for i in range(1044):
        if PRESUP[i]:
            Indices.append(i)
    #        aa = savitzky_golay(dataC[i,902:1040],11,2)
    #        dataCC[i,902:1040] = dataC[i,902:1040] - aa
            #dataC[i,902:1040] - TEST[i]    
    

    Cons = group_consecutives(Indices)
    
    Consec = []
    for i in Cons:
        if len(i)>1:
            Consec.append(i)    
    
    
    return Consec



    
    

def Correct_Amplifier(data):
    
#    data[:,1029:1058] = 0 
    
    for col in range(1566):
        c = data[:,col]
        
        SIG = Sigma_Clip(c)
                
        filtered = c*SIG
        
        z = np.polyfit(range(1044), filtered, 7)
        p = np.poly1d(z)
        
        Fitted = p(range(1044))
        
        data[:,col] = c - Fitted
        
    return data


def Display(name):
    iraf.set(stdimage='imtgmos2')
    gmos.gdisplay(name)
    
def Create_Query(IND_KEY = [7,3,4,5,6],IND_COND = []):
        
        
    SQL_Query =  'SELECT '
    for ind in IND_KEY:
        SQL_Query = SQL_Query + KEY_LIST[ind] +','
    SQL_Query = SQL_Query[0:-1] 
    SQL_Query = SQL_Query + ' FROM obslog '
    
    
    if IND_COND:
        SQL_Query = SQL_Query + 'WHERE '
        for ind in IND_COND:
            SQL_Query = SQL_Query + KEY_LIST[ind] +'=? '
            SQL_Query = SQL_Query + 'AND ' 
        SQL_Query = SQL_Query[:-4]
        
    return SQL_Query


def Create_database(raw = 'raw', Database = 'obsLog.sqlite3'):
    
    print("=== Creating Observer log file===")    
    
    os.system('python /Users/maximedevogele/Documents/Gemini/obslog2.py' + ' ' + Database + ' ' + raw + '/')    



def Get_From_database(**kw):
    
    '''
        Allow to display the information contain in the file database
        
        Arguments:      Type                default values
            dbFile      String              './raw/obsLog.sqlite3'  File containing the database
            IND_KEY     List of integers    [7,3,4,5,6]             Keys to display
            IND_COND    List of integers    []                      Keys to put conditions
            VALUE_COND  Tupple of Strings   ()                      Conditions on the keys
            
        Keys:
            [0]     obs_id
            [1]     use_me
            [2]     N_ext
            [3]     File
            [4]     DateObs
            [5]     TimeObs
            [6]     Instrument
            [7]     Object
            [8]     RA
            [9]     Dec
            [10]    ObsType
            [11]    ObsClass
            [12]    CcdBin
            [13]    RoI
            [14]    NodMode
            [15]    NS_Shift
            [16]    NS_Cycles
            [17]    DTA_Xoffset
            [18]    Filter1
            [19]    Filter2
            [20]    Disperser
            [21]    AperMAsk
            [22]    MaskType
            [23]    Rotator
            [24]    CentWave
            [24]    T_exp
            [26]    Airmass
            [27]    Quality
            [28]    Xoffset            
            [29]    Yoffset

            
        Examples:   
            Display_database() 
                Uses default values and display keys 
                7 the Object,3 the file name,4 the date of observation,
                5 the time of observation, and 6 the Instrument used
            
            Display_database(Ind_Keys = [0,1,2,3])  
                displays only the keys 1,2,3,4
            
            Display_database(Ind_Keys = [0,1,2,3], IND_COND = [7], VALUE_COND = (2017 QB,) : 
                displays only the keys 1,2,3,4 and only entries where the object is 2017 QB
            
    '''
    
    
    
    
    # Define the default value for the query to the database
    defaults = dict(dbFile = './raw/obsLog.sqlite3',
                    IND_KEY = [7,3,4,5,6],
                    IND_COND = [],
                    VALUE_COND = ()
                    )   
     
    # If parameter has not been provided by the user,
    # assign the defaults values to the kw dictionary
    for key in defaults.keys():
        if not kw.has_key(key):
            kw[key] = defaults[key]

    dbFile = kw.get('dbFile')    
    IND_KEY = kw.get('IND_KEY')
    IND_COND = kw.get('IND_COND')
    VALUE_COND = kw.get('VALUE_COND')
    
    
    SQL_Query = Create_Query(IND_KEY = IND_KEY,IND_COND = IND_COND)
    

    db = sqlite3.connect(dbFile)
    c = db.cursor()
    c.execute(SQL_Query, VALUE_COND)
    all_rows = c.fetchall()
    
    PRINT_STRING = ''
    for number in range(len(IND_KEY)):
        PRINT_STRING = PRINT_STRING + '{' + str(number) + '} '
        
    #Header = [KEY_LIST[i] for i in IND_KEY]

    return all_rows
    
    db.close()


def Display_database(**kw):
    
    '''
        Allow to display the information contain in the file database
        
        Arguments:      Type                default values
            dbFile      String              './raw/obsLog.sqlite3'  File containing the database
            IND_KEY     List of integers    [7,3,4,5,6]             Keys to display
            IND_COND    List of integers    []                      Keys to put conditions
            VALUE_COND  Tupple of Strings   ()                      Conditions on the keys
            
        Keys:
            [0]     obs_id
            [1]     use_me
            [2]     N_ext
            [3]     File
            [4]     DateObs
            [5]     TimeObs
            [6]     Instrument
            [7]     Object
            [8]     RA
            [9]     Dec
            [10]    ObsType
            [11]    ObsClass
            [12]    CcdBin
            [13]    RoI
            [14]    NodMode
            [15]    NS_Shift
            [16]    NS_Cycles
            [17]    DTA_Xoffset
            [18]    Filter1
            [19]    Filter2
            [20]    Disperser
            [21]    AperMAsk
            [22]    MaskType
            [23]    Rotator
            [24]    CentWave
            [24]    T_exp
            [26]    Airmass
            [27]    Quality
            
        Examples:   
            Display_database() 
                Uses default values and display keys 
                7 the Object,3 the file name,4 the date of observation,
                5 the time of observation, and 6 the Instrument used
            
            Display_database(Ind_Keys = [0,1,2,3])  
                displays only the keys 1,2,3,4
            
            Display_database(Ind_Keys = [0,1,2,3], IND_COND = [7], VALUE_COND = (2017 QB,) : 
                displays only the keys 1,2,3,4 and only entries where the object is 2017 QB
            
    '''
    
    
    
    
    # Define the default value for the query to the database
    defaults = dict(dbFile = './raw/obsLog.sqlite3',
                    IND_KEY = [7,3,4,5,6],
                    IND_COND = [],
                    VALUE_COND = ()
                    )   
     
    # If parameter has not been provided by the user,
    # assign the defaults values to the kw dictionary
    for key in defaults.keys():
        if not kw.has_key(key):
            kw[key] = defaults[key]

    dbFile = kw.get('dbFile')    
    IND_KEY = kw.get('IND_KEY')
    IND_COND = kw.get('IND_COND')
    VALUE_COND = kw.get('VALUE_COND')
    
    
    SQL_Query = Create_Query(IND_KEY = IND_KEY,IND_COND = IND_COND)
    

    db = sqlite3.connect(dbFile)
    c = db.cursor()
    c.execute(SQL_Query, VALUE_COND)
    all_rows = c.fetchall()
    
    PRINT_STRING = ''
    for number in range(len(IND_KEY)):
        PRINT_STRING = PRINT_STRING + '{' + str(number) + '} '
        
    Header = [KEY_LIST[i] for i in IND_KEY]
    print(PRINT_STRING.format(*Header))
    for row in all_rows:
        print(PRINT_STRING.format(*row))

    db.close()


def Get_ObjectList(dbFile = './raw/obsLog.sqlite3'):
    
    qr = {'use_me':1,
           'Object':'*'}
    
    ObjectList = list(set(fs.objectListQuery(dbFile,fs.createQuery('object', qr, Date = False),qr)))

    print(ObjectList)
    
    return ObjectList

def Get_DateList(dbFile = './raw/obsLog.sqlite3'):
    
    qr = {'use_me':1,'Object':'*',
           'DateObs':''}
    
    DateList = list(set(fs.dateListQuery(dbFile,fs.createQuery('dateobs', qr, Date = False),qr)))
    
    print('The targets were observed on the: \n')
    for i in DateList: print(str(i))
    print('\n')
    
    return DateList
    
    
def Get_BinList(dbFile = './raw/obsLog.sqlite3'):
    
    qr = {'use_me':1,'Object':'*',
           'DateObs':'','CentWave':''}
    
    BinList = list(set(fs.ccdbinListQuery(dbFile,fs.createQuery('ccdbin', qr, Date = False),qr)))
   
    print(BinList)

    return BinList


def Get_ObjectsFileList(dbFile = './raw/obsLog.sqlite3'):
    
    ObjList = Get_ObjectList()
    
    ObjFileList={}
    for Object in ObjList:
        qr = {'use_me':1,
              'Object':Object,
              'DateObs':''}
        ObjectFileList = list(set(fs.fileListQuery(dbFile,fs.createQuery('sciSpecSimp', qr, Date = False),qr)))
        ObjFileList[Object] = ObjectFileList
    
    print(ObjFileList)

    return ObjFileList
    

def Check(filename):
    
    instruments =[] 
    for elem in filename:
        if not os.path.isfile(elem):
            raise IOError('No such file: ' + str(elem))
            
        try:
            hdulist = fits.open(elem, ignore_missing_end=True)
        except IOError:
            print('ERROR: cannot open file %s' % elem)
            continue

        header = hdulist[0].header
        for key in _SP_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])
                break

    if len(np.unique(instruments)) == 0:
        raise KeyError('cannot identify telescope/instrument; please update' + \
                       '_SP_conf.instrument_keys accordingly')

    if len(np.unique(instruments)) > 1:
        raise Warning('More than one instrument identified')    
    
    
    if len(filename) == 0:
        raise IOError('cannot find any data...')
        
def bckg_sub(files,out):

    hdulist = fits.open(files)
    image = hdulist[0].data
    
    Mask = np.zeros((2088,3132))
    
    Mask[:,:] =  0
    
    Mask[0:600,:] = 1
    Mask[CHIP_GAP2_B2[0][0]:CHIP_GAP2_B2[0][1],:] = 1
    Mask[CHIP_GAP2_B2[1][0]:CHIP_GAP2_B2[1][1],:] = 1
    Mask[800:1200,:] = 1
    Mask[1400:,:] = 1
    
    Mask = Mask.astype(bool)
    
    image2 = ma.masked_array(image,mask = Mask)
    
    
    X = range(2088)
    X = np.array(X)
    
    for i in range(3132):
        image3 = Sigma_Clip(image2[:,i])
        XX = ma.masked_array(X,image3.mask)
        z = np.ma.polyfit(XX,image3 , 2)
        p = np.poly1d(z)
        image[:,i] = image[:,i] - p(X)
        
    hdulist.writeto(out)
    

def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError("Cannot covert {} to a bool".format(s))
