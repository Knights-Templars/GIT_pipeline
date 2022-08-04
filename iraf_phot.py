import os

# Image directory

#================================================================================#

data_dir = '/home/anirban.dutta/SN2022erq_Reduction/SA110/For_Phot/SN_ALIGNED/'

os.chdir(data_dir)

common_text = 'wcs_*.fits'

#================================================================================#

import re
import glob
import shutil
import datetime
import warnings
import numpy as np
import pandas as pd
from pyraf import iraf
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

#================================================================================#

config_sex = '/home/anirban.dutta/astromatic/config_sextractor.sex'
param_sex = '/home/anirban.dutta/astromatic/default.param'

#================================================================================#

# Stop showing warnings. Not a good practice !
warnings.filterwarnings("ignore")

#================================================================================#

OBJECT_NAME = 'sa110_340'
OBJECT_RA = '18:41:28.44'
OBJECT_DEC = '+00:15:23.0'
TYPE = 'SNIa'
REDSHIFT = 0.01
DISCOVERY_DATE = '2021-08-23 08:55:40.800'
HOST_GALAXY = 'Unknown'
OBJECT = 'sa110_340'

#================================================================================#

# Observatory Details
#--------------------------------------------------------------------------------#
observatory = 'iao'
obs_lat = '32:46:46'
obs_long = '78:57:51'
latitude = 32.7794 * u.deg
longitude = 78.9642 * u.deg
altitude = 4500 * u.m
tz = +5.5 * u.hour
name = 'Indian Astronomical Observatory, Hanle'
#--------------------------------------------------------------------------------#
# New CCD Specification

read_noise_new = 12.0  # 5.75
ccd_gain_new = 1.04   # 0.28
data_max_new = 55000  # 700000

# Old CCD Specification

read_noise_old = 4.87
ccd_gain_old = 1.22
data_max_old = 55000
#--------------------------------------------------------------------------------#


# Important Image Header Keywords

RA_key = 'RA'
DEC_key = 'DEC'
UT_key = 'UT'
DATE_key = 'DATE-OBS'
FILTER_key = 'FILTER'
AIRMASS_key = 'AIRMASS'
EXPTIME_key = 'EXPTIME'
#--------------------------------------------------------------------------------#

# load IRAF packages
# _doprint = 0 does not print all the subpackages under eack package
#--------------------------------------------------------------------------------#
iraf.noao(_doprint = 0)
iraf.imred(_doprint = 0)
iraf.ccdred(_doprint = 0)
iraf.digiphot(_doprint = 0)
iraf.daophot(_doprint = 0)
iraf.ptools(_doprint = 0)
iraf.ccdred.instrument = 'ccddb$kpno/direct.dat'

#================================================================================#

def remove_file(file_name):
    '''
    file_name: name of the file to remove
    '''
    os.remove(file_name)

def remove_similar_files(common_text):
    
    '''
    common_text: A string (e.g. *.fits, *.list) used for removing
    similar kind of files
    '''
    
    for file_name in glob.glob(common_text):
        remove_file(file_name)
        
def group_similar_files(text_list, common_text, exceptions = ''):
    
    '''
    text_list: A text file used to store list of files
    common_text: A string (e.g. *.fits, *.list) used for grouping similar
    kinds of files
    exceptions: string of file name to exclude in grouping
    
    returns: list of grouped files
    '''
    
     
    list_files = glob.glob(common_text)
    if exceptions != '':
        list_exceptions = exceptions.split(',')
        for text in list_exceptions:
            list_files = filter(lambda x: not re.search(text, x), list_files)
        
    list_files.sort()
    if len(text_list) != 0:
        with open(text_list, 'w') as f:
            for file_name in list_files:
                f.write(file_name+'\n')
                
    return list_files         


def text_list_to_python_list(text_file):
    
    '''
    text_file: A text file from which a python list will be made
    
    returns: A python list
    '''
    
    if os.path.isfile(text_file):
        with open(text_file, 'r+') as f:
            python_list = f.read().split()
            
    return python_list

def python_list_to_text_list(python_list, text_file):
    
    '''
    python_list: A python list
    text_file: A file for saving a python list
    '''
    
    with open(text_file, 'w') as f:
        for element in python_list:
            f.write(str(element) + '\n') 
        f.close()
        
#================================================================================#

def imexam_fwhm(text_list, coord_file, log_imexam = 'log_imexam'):
    
    '''
    Run the Imexam task in iraf
    text_list: a list of images on which imexam will be run
    coord_file: A file of stars coordinate on which imexam will be run
    log_imexam: A file to store the imexam results
    
    '''
    
    list_files = text_list_to_python_list(text_list)
    
    task = iraf.images.tv.imexam
    task.unlearn()
    
    for file_name in list_files:
        task.ncoutpu = 101
        task.nloutpu = 101
        task.image = ''
        task.logfile = log_imexam
        task.keeplog = 'yes'
        task.defkey = 'a'
        task.autored = 'yes'
        task.allfram = 'yes'
        task.nframes = 0
        task.ncstat = 5
        task.nlstat = 5
        task.graphcu = ''
        task.imagecu = coord_file
        task.wcs = 'logical'
        task.xformat = ''
        task.yformat = ''
        task.graphic = 'stdgraph'
        task.use_dis = 'no'
        task.mode = 'ql'
        
        task(input = file_name, frame = 1)
        
        
def datapars(fwhm_value, data_max=data_max_new, read_noise=read_noise_new, ccd_gain=ccd_gain_new, 
             exposure=EXPTIME_key, airmass=AIRMASS_key, filter_=FILTER_key):
    '''
    Edit data parameters for photometry
    fwhm_value: median fwhm of the stars in the field
    data_max: saturation value of the CCD
    read_noise: read noise (e-)of the CCD
    ccd_gain: gain(e-/ADU) of the CCD
    exposure: Exposure time of the image (header keyword)
    airmass: Airmass at which observation was done (header keyword)
    filter: Bandpass of observation (header keyword)
    ut: Universal Time (header keyword)
    '''
    
    task = iraf.noao.digiphot.daophot.datapars
    task.unlearn()
    
    task.scale = 1.0
    task.fwhmpsf = float(fwhm_value)
    task.emissio = 'yes'
    task.sigma = 'INDEF'
    task.datamin = 'INDEF'
    task.datamax = data_max
    task.noise = 'poisson'
    task.ccdread = ''
    task.gain = ''
    task.readnoi = read_noise
    task.epadu = ccd_gain
    task.exposur = exposure
    task.airmass = airmass
    task.filter = filter_
    task.obstime = exposure
    task.itime = 1.
    task.xairmass = 'INDEF'
    task.ifilter = 'INDEF'
    task.otime = 'INDEF'
    task.mode = 'ql'
    
    
def centerpars(center='centroid'):
    
    '''
    Edit Centering parameters for photometry
    centering algorithm: centroid, gauss, none, ofilter
    '''
    
    task = iraf.noao.digiphot.daophot.centerpars
    task.unlearn()
    
    task.calgori = center  # centering algorithm
    task.cbox = 3
    task.cthresh = 0
    task.minsnra = 1.
    task.cmaxite = 10
    task.maxshif = 0.5
    task.clean = 'no'
    task.rclean = 1.
    task.rclip = 2.
    task.kclean = 3.
    task.mkcente = 'no'
    task.mode = 'ql'
    

def fitskypars(fwhm_value, n=3, mode='mode'):
    
    '''
    Edit Sky parameters for photometry
    fwhm_value: median fwhm of the stars in the field 
    n: Number of fwhm in pixels the sky annulus will be
    mode: sky algorithm | median, mode, center, gauss
    '''
    task = iraf.noao.digiphot.daophot.fitskypars
    task.unlearn()
    
    task.salgori = mode
    task.annulus = n*float(fwhm_value)
    task.dannulus = 3
    task.skyvalu = 0.
    task.smaxite = 10
    task.sloclip = 0.
    task.shiclip = 0.
    task.snrejec = 50
    task.sloreje = 3
    task.shireje = 3
    task.khist = 3
    task.binsize = 0.1
    task.smooth = 'no'
    task.rgrow = 0.
    task.mksky = 'no'
    task.mode = 'ql'
    

def photpars(aperture_values):
    
    '''
    Edit the photometry parameters.
    aperture_values: Values of the aperture at which to perform
                    photmetry.
    Can be single valued, list of aperture values
    e.g. 5, [5,10,15,20], 5:15:1
    
    '''
    
    task = iraf.noao.digiphot.daophot.photpars
    task.unlearn()
    
    task.weighti = 'constant'
    task.apertur = aperture_values
    task.zmag = 25
    task.mkapert = 'no'
    task.mode = 'ql'
    
    
def daopars(fwhm_value, psf_radius, recenter='yes'):
    
    '''
    Edit the daophot fitting parameters for photometry
    fwhm_value: median fwhm of the stars in the field
    psf_radius: A constant multiplicative term
    
    '''
    
    psf_aperture = float(psf_radius)*float(fwhm_value)
    
    task = iraf.noao.digiphot.daophot.daopars
    task.unlearn()
    
    task.functio = 'moffat25'   # gauss, moffat15, moffat25, lorentz, penny1, penny2, auto
    task.varorde = 2
    task.nclean = 3
    task.fitsky = 'yes'
    task.recente = 'yes'
    task.matchra = float(fwhm_value)
    task.psfrad = psf_aperture
    task.fitrad = 1.5*float(fwhm_value)
    task.sannulu = 3*float(fwhm_value)
    task.wsannul = 3
    
    
def phot(file_name, coord_file):
    
    '''
    file_name: name of the file on which to run photometry
    coord_file: A file of stars coordinate on which photometry will be run
    output: image.mag.?
    '''
    
    task = iraf.noao.digiphot.daophot.phot
    task.unlearn()
    
    task.interac = 'no'
    task.radplot = 'no'
    task.verbose = 'no'
    task.verify = 'no'
    task.update = 'no'
    
    task(image = file_name, coords = coord_file, output = 'default')
    

def pstselect(file_name, magfile_name, fwhm_value, psf_radius, data_max):
    
    '''
    file_name: name of the file on which to run the task pstselect
    magfile_name: photometry file from phot task
    fwhm_value: median fwhm of the stars in the field
    psf_radius: A constant multiplicative term
    
    image: Image for which to build psf
    photfile: Input mag file | image.mag.?
    pstfile: output psf star list | image.pst.?
    
    '''
    
    datapars(fwhm_value)
    daopars(fwhm_value, psf_radius)
    
    task = iraf.noao.digiphot.daophot.pstselect
    task.unlearn()
    
    task.interac = 'no'
    task.verify = 'no'
    task.verbose = 'no'
    task.update = 'no'
    
    task(image = file_name, photfile = magfile_name, pstfile = 'default', 
         maxnpsf = 25)
    
    
    
def psf(file_name, magfile_name, pstfile_name, fwhm_value, psf_radius):
    
    '''
    file_name: name of the file on which to run the task psf
    magfile_name: photometry file from phot task
    pstfile_name: pstselect file from the task pstselect | image.pst.?
    fwhm_value: median fwhm of the stars in the field
    psf_radius: A constant multiplicative term
    
    image: Input image for which to build psf
    photfile: Input photomery file | image.mag.?
    pstfile: Input psf star list | image.pst.?
    psfimage: Output psf image | image.psf.?
    opstfile: Output psf star list | image.pst.?
    groupfile: Output psf star group file | image.psg.?
    '''
    
    
    datapars(fwhm_value)
    daopars(fwhm_value, psf_radius)
    
    task = iraf.noao.digiphot.daophot.psf
    task.unlearn()
    
    task.plotfil = ''
    task.matchby = 'yes'
    task.interac= 'no'
    task.showplo = 'no'
    task.verify = 'no'
    task.update = 'no'
    task.verbose = 'no'
    
    task(image = file_name, photfile = magfile_name, pstfile = pstfile_name, 
         psfimage = 'default', opstfile = 'default', groupfil = 'default')
    
    
def allstar(file_name, magfile_name, psffile_name, fwhm_value, datamax, 
            psf_radius, recenter='yes'):
    
    '''
    file_name: name of the file on which to run allstar
    magfile_name: photometry file from phot task
    psffile_name: pstselect file from the task psf | image.psf.?
    fwhm_value: median fwhm of the stars in the field
    data_max: saturation value of the CCD
    psf_radius: A constant multiplicative term
    
    image: Image corresponding to photometry
    photfile: Input photometry file | image.mag.?
    psfimage: Psf image | image.psf.?
    allstarf: Output photometry file | image.als.?
    rejfilef: Output rejection file | image.arj.?
    subimage: Subtracted image | image.sub.?
    '''
    
    datapars(fwhm_value)
    daopars(fwhm_value, psf_radius, recenter)
    
    
    task = iraf.noao.digiphot.daophot.allstar
    task.unlearn()
    
    task.verbose = 'no'
    task.verify = 'no'
    task.update = 'no'
    
    task(image = file_name, photfile = magfile_name, psfimage = psffile_name, 
         allstarf = 'default', rejfile = 'default', subimage = 'default')

def sort_als(file_name):

    task = iraf.noao.digiphot.ptools.psort
    task.unlearn()

    task(infiles=file_name, field='ID', ascend='yes')


def run_sextractor(ctext):

    for text in ['*.list']:
        remove_similar_files(common_text=text)

    file_list='ImageList.list'
    list_files=group_similar_files(file_list, common_text=ctext)

    for file_name in list_files:
        command="sex %s -c %s -CATALOG_NAME %s -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME %s -MAG_ZEROPOINT 25.0" % (file_name, config_sex, file_name+'.ldac', param_sex)
        os.system(command)
        print('Executing command: %s\n' % command)        

def get_table_from_ldac(filename, frame=1):
    
    """
    Load an astropy table from a fits_ldac by frame 
    
    """
    from astropy.table import Table
    if frame>0:
        frame = frame*2
    tbl = Table.read(filename, hdu=frame)
    return tbl          


def calculate_fwhm(textlist_files, coord_file='stars.coo', 
                   log_imexam='log_imexam'):
    
    '''
    textlist_files: text file of a list of images on which to run imexam
    coord_file: A file of stars coordinate on which photometry will be run
    log_imexam: A file to store the imexam results
    '''
    
    list_files = text_list_to_python_list(textlist_files)
    imexam_fwhm(textlist_files, coord_file, log_imexam)
    coord_df = pd.read_csv(coord_file, sep='\s+', comment='#', header=None)
    rows, columns = coord_df.shape
    # This is applicable for imexam log text file
    n = (2+rows)*len(list_files)
    columns = ['COL', 'LINE', 'X', 'Y', 'R', 'MAG', 'FLUX', 'SKY', 'PEAK',
               'E', 'PA', 'BETA', 'ENCLOSED', 'MOFFAT', 'DIRECT']
    imexam_df = [pd.read_csv(log_imexam, sep='\s+', names=columns, skiprows=i, 
                nrows=rows, header=None, comment='#') for i in range(2, n, rows+2)]
    list_fwhm = [df['MOFFAT'].values for df in imexam_df]
    list_median_fwhm = []
    for fwhm_values in list_fwhm:
        fwhm_list = [float(value) for value in fwhm_values if value != 'INDEF']
        print (type(fwhm_list))
        median = np.median(fwhm_list)
        list_median_fwhm.append(round(median, 1))
        
    return list_median_fwhm


def calculate_airmass(textlist_files):
    
    list_airmass = []
    
    list_files = text_list_to_python_list(textlist_files)
    hct = EarthLocation.from_geodetic(lat=latitude, lon=longitude, height=altitude)
    for file_name in list_files:
        hdu = fits.open(file_name, mode='update')
        header = hdu[0].header

        if 'TM_START' in header.keys():
            date_obs = header['DATE-OBS']
            time_start = header['TM_START']
            
            ra_ = OBJECT_RA
            dec_ = OBJECT_DEC
            c = SkyCoord(ra_, dec_, unit=(u.hourangle, u.deg))
            ra = c.ra.deg
            dec = c.dec.deg
            time_utc = str(datetime.timedelta(seconds=int(time_start)))
            datetime_utc = date_obs+' '+time_utc
            time = Time(datetime_utc) 
        else:
            date_time = header['DATE-OBS'].split('T')
            time_obj = date_time[0]+' '+date_time[1] 
            time = Time(time_obj)
            
            ra = OBJECT_RA
            dec = OBJECT_DEC
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            ra = c.ra.deg
            dec = c.dec.deg
        

        coord = SkyCoord(ra, dec, unit='deg')
        altaz_ = coord.transform_to(AltAz(obstime=time, location=hct))
        airmass = altaz_.secz.value
        print('The image %s has been observed at an airmass of %f'%(file_name, airmass))
        
        list_keywords = ['RA', 'DEC', 'AIRMASS']
        dict_header = {'RA':ra, 'DEC':dec,
                      'AIRMASS': airmass}
        
        for key in list_keywords:
            if key in header.keys():
                header.remove(key, remove_all=True)
            header.append(card=(key, dict_header[key]))
            
        hdu.flush()
        hdu.close()
        
        list_airmass.append(airmass)
        
    return list_airmass


def use_aperture(fwhm, n=4, cog=False, single_ap=True):
    
	"""
    Calculates apertures to be calculated in terms of 'Pixels' from a string supplying apertures
    in terms of FWHM value of the image.
    Args:
       
        fwhm        : FWHM of the image to which photometry is being done
        n           : No. of fwhm for the bigger aperture
    Returns:
        aper_values : String containing apertures to be used for photometry
    """
	
	smaller_ap = fwhm
	larger_ap = n*float(fwhm)
	
	if not single_ap:
		if cog:
			aperture = str(smaller_ap)+':'+str(larger_ap)+':'+str(1)
		else:
			aperture = str(smaller_ap)+':'+str(larger_ap)
	else:
		aperture = fwhm
	
	return aperture

def aper_phot(file_name, fwhm, coord_file, cog, single_ap, 
		data_max, center='centroid'):
    """
    Performs aperture photometry (PHOT task) on the files in the list 'list_files'.
    Selects candidate stars from the coordinate file 'coord_file'.
    Args:
        textlist_files  : List of all FITS files on which aperture photometry 
        		   is to be performed.
        textlist_fwhm   : List of Mean FWHM values of all the FITS files
        coord_file      : Name of the coordinate file containing candidate star
        phot_radius     : String containing the apertures at which photometry 
        		   is to be done("1,4").
        data_max        : Maximum good pixel value
    Returns:
        None
    """
    #list_files = text_list_to_python_list(textlist_files)

    #for index in range(0, len(list_files)):
    aperture_values = use_aperture(fwhm, cog=cog, 
        				single_ap=single_ap)

    datapars(fwhm, data_max)
    centerpars(center)
    fitskypars(fwhm)
    photpars(aperture_values)
    phot(file_name=file_name, coord_file=coord_file)
    
    display_text("Aperture Photometry Is Completed For Aperture Values (x FWHM): {0}".format(aperture_values))
    	
def psf_phot(file_name, fwhm, mag_suffix='.mag.1', mag_apply = '.mag.1',psf_radius= 2.0, 
	     psf_apply=2.0, data_max=data_max_new, recenter='yes'):
    """
    Performs PSF Photometry on the text list 'list_files'.
    Args:
        textlist_files  : List of all FITS files on which PSF photometry is to 
        		   be performed.
        textlist_fwhm   : List of Mean FWHM values of all the FITS files
        mag_suffix      : Suffix of the mag files from PHOT task to be used for 
        		   selecting candidate stars
        psf_radius      : PSF fit radius in units of FWHM
        data_max        : Maximum good pixel value
    Returns:
        None
    """
    global run_count

    #list_files = text_list_to_python_list(textlist_files)

    #for index in range(0, len(list_files)):
    file_mag = file_name + mag_suffix
    file_mag2 = file_name + mag_apply
    pstselect(file_name, file_mag, fwhm, psf_radius, 
                 data_max)

    file_pst = file_name + '.pst.' + str(run_count * 2 - 1)
    psf(file_name, file_mag, file_pst, fwhm, psf_radius)

    file_psf = file_name + '.psf.' + str(run_count) + '.fits'
    allstar(file_name, file_mag2, file_psf, fwhm, data_max, 
                psf_apply, recenter)

    #run_count += 1
    display_text("PSF Photometry Is Completed For PSF Radius = {0} * FWHM".format(psf_radius))
  
def display_text(text):

    print('#'+'-'*(10+len(text))+'#')
    print('#'+('-'*5) + str(text) +('-'*5)+'#')
    print('#'+'-'*(10+len(text))+'#')   


# Plot SN

def read_data(image_name):

    '''
    image_name: name of the image

    Returns: data_array, image_header
    '''
    
    image = fits.open(image_name)
    image_data = image[0].data
    image_header = image[0].header
    
    return image_data, image_header




#================================================================================#

remove_resfile = True

if remove_resfile:
    for text in ['*.als*', '*.arj*','*.psf*', '*.psg*', '*.pst*', 
                 '*.sub*', 'OUTPUT*', 'output*', '*.mag*', 'list_*', '*.coo', '*.list']:
        remove_similar_files(common_text=text)

#================================================================================#

text_list = 'list_files'
textlist_fwhm = 'list_fwhm'
list_files = group_similar_files(text_list, common_text=common_text)

print ("The list of files on which photometry will be performed")
print (list_files)

list_airmass = calculate_airmass(textlist_files=text_list)

list_fwhm =[]

run_sextractor(ctext=common_text)

for file_name in list_files:

    if os.path.exists("stars.coo"):
        os.remove("stars.coo")

    sourceTable = get_table_from_ldac(file_name+'.ldac')

    cleanSources = sourceTable[(sourceTable['FLAGS']==0) & (sourceTable['FWHM_WORLD'] < 2) & 
                           (sourceTable['XWIN_IMAGE']<3500) & (sourceTable['XWIN_IMAGE']>500) &
                           (sourceTable['YWIN_IMAGE']<3500) &(sourceTable['YWIN_IMAGE']>500) &
                          (sourceTable['ELLIPTICITY'] < 0.4)]
                          

    fwhm = np.median(sigma_clip(cleanSources['FWHM_IMAGE'])).tolist()
    print ('List of Fwhm for all the images')
    print ('#------------------------------#')
    print (fwhm)
    list_fwhm.append(fwhm)     

    print ("Performing photometry for image %s"%file_name)
    coord_file_stars = 'stars.coo'    

    x = cleanSources['XWIN_IMAGE']
    y = cleanSources['YWIN_IMAGE']

    col_format = "{0:>4f}{1:15f}\n"
    with open("stars.coo", 'w') as f:
        for x in zip(x, y):
            f.write(col_format.format(*x))
    f.close()  



    aper_phot(file_name, fwhm, coord_file_stars, cog = True, single_ap = False, 
	        data_max=data_max_new)                    


    run_count =1

    psf_phot(file_name, fwhm, data_max=data_max_new)

'''
for index, file_name in enumerate(list_files):

    if os.path.exists("sn.coo"):
        os.remove("sn.coo")

    data, header = read_data(file_name)


    corr_deg = SkyCoord('17 56 02.513 +18 21 14.07', frame='icrs', unit=(u.hourangle, u.deg))

    ra_deg = corr_deg.ra.value
    dec_deg = corr_deg.dec.value

    print (ra_deg, dec_deg)


    w = WCS(header)

    x_sn, y_sn = w.all_world2pix(ra_deg, dec_deg, 1)

    print (x_sn, y_sn)


    coordinate_file ="sn.coo"

    if os.path.exists(coordinate_file):
        os.remove(coordinate_file)  

    col_format = "{0:>4f}{1:15f}\n"
    with open("sn.coo", 'w') as f:
    #for x in zip(x_sn, y_sn):
        f.write(col_format.format(x_sn, y_sn))
    f.close()

    sort_als(file_name+'.als.1')


    coord_file_sn = 'sn.coo'

    print ("Fwhm at which photometry will be performed", list_fwhm[index])

    aper_phot(file_name, list_fwhm[index], coord_file_sn, cog = False, single_ap = False, 
	            data_max=data_max_new)    

    psf_phot(file_name, list_fwhm[index], mag_suffix='.mag.1', mag_apply='.mag.2',
	         data_max=data_max_new)    
'''
