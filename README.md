# SP: Spectroscopic Pipeline

2020-02-06 [QY]: creating the readme.

## Setup

### Add to .bashrc:

```
export PATH="$HOME/Library/SP:$PATH"
export SPECPIPEDIR=$HOME/Library/SP
```

### Set files in the directory to executable:

```
chmod 700 *.py
```

### Install dependencies

```
conda install -c conda-forge imreg_dft pillow astropy matplotlib astroquery future itertools simplejson
```

### Requirement

scipy <=1.2.0 as 1.3.0+ does not have toimage(). Python 3.5 seems to work.

## Run-through: DeVeny

1. Renaming the data according to the FITS header: `SP_Prepare.py *fits` [note: SP will automatically make a copy of the raw data, as it will rename/rewrite all the data.]
2. Print out a summary of all data: `SP_Log.py *fits`
3. Creating master bias: `SP_Bias.py *BIAS* -o MasterBias.fits`
4. Creating master flat: `SP_Flat.py *FLAT* -b MasterBias.fits -o MasterFlat.fits -m none` [note: be aware of the -m parameter, if there is only one set of flats, “none” is sufficient.]
5. Bias subtraction and flat division: `SP_Preproc.py *MyObject* -b MasterBias.fits -f MasterFlat_1.fits` [note: do the same thing for other targets and standard stars]
6. Background subtraction: `SP_BckgSub.py *Procc.fits`
7. Spectrum extraction: `SP_Extract.py *Bckg.fits`
8. Combine target: `SP_Combine.py *MyObject*_Extracted.txt -o MyObject_Comb.txt`
9. Combine standard stars: `SP_Combine.py *Standard*_Extracted.txt -o SA_Comb.txt`
10. Wavelength calibration (target): `SP_WavCal.py MyObject_Comb.txt -a *ARCS* -m template -o MyObject.sp`
11. Wavelength calibration (standard star): `SP_WavCal.py Standard_Comb.txt -a *ARCS* -m template -o Standard.sp`
12. Telluric correction: `SP_TellCorr.py MyObject.sp -s Standard.sp -i Deveny`
