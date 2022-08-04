# A Python Code for performing Alignment of Astronomical Images having WCS 
# Author - Anirban Dutta
# Version - 2.0
# Date: 20/10/2019

#---------------------------------------------------------------------------------------------------------------------#
# Import the necessary modules and packages

import os
import numpy as np
import sys
import shutil
import glob
from astropy.io import fits
import subprocess
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import warnings
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.io import ascii
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import re


#from pyraf import iraf



#---------------------------------------------------------------------------------------------------------------------#
# FEW PARAMETERS WHICH THE USER CAN CHANGE

Astrometry = False
make_catalog = True
Precise_align = True
coarse_align = False
align_check = False
TELESCOPE = input("Enter the Telescope name(HCT/GROWTH) for which you want to run the code [HCT/GROWTH]:")

# Details regarding the OBJECT under study #

OBJECT='OBJECT'
OBJECT_NAME='sa110_340'
Right_Ascension='RA'
Declination='DEC'
RA='18:41:28.44' 
DEC='+00:15:23.0'

# Change the Working Directory

working_directory = '/home/anirban.dutta/SN2022erq_Reduction/SA110/For_Phot/'
os.chdir(working_directory)
print (working_directory)

#---------------------------------------------------------------------------------------------------------------------#

# Test whether SCAMP, SWARP and SEXTRACTOR are installed properly

dependencies=['swarp', 'scamp', 'sex']

def test_dependencies(dep):
    try:
       subprocess.Popen(dep, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
       print("%s is installed properly. OK" % dep)
       return 1
    except ImportError:
        print("==%s IS NOT INSTALLED PROPERLY" % dep)
        return 0
        
#i=0
#for dep in dependencies:
#    i+=test_dependencies(dep)
#print("\n%i out of %i dependencies installed properly." % (i, len(dependencies)))
#if i!=len(dependencies):
#    print('**********Please Install the programs before continuing**********')
#else:
#    print('**********You are ready to continue**********')
    
#---------------------------------------------------------------------------------------------------------------------#

# Path for config files
config_sex = '/home/anirban.dutta/astromatic/config_sextractor.sex'
param_sex = '/home/anirban.dutta/astromatic/default.param'
config_scamp = '/home/anirban.dutta/astromatic/config.scamp'
config_swarp = '/home/anirban.dutta/astromatic/config.swarp'


#---------------------------------------------------------------------------------------------------------------------#
# Telescope and CCD Specifications:

Telescope_Place= 'Hanle'
Telescope_1 = 'HCT'
CCD_name_1='RTS2'
read_noise_1=4.87 # electrons
gain_1=1.22 # electron/ADU
gain_1_new=0.28 # electron/ADU
read_noise_1_new=5.75 # electrons
data_max_new = 700000
data_max_1=55000 
pxscale_HCT = 0.296 # arcsec/pixel
scale_low_1 = 0.1
scale_high_1 = 0.5
Telescope_2 = 'GROWTH'
CCD_name_2='RTS2'

read_noise_2=12.0 # electrons
gain_2=1.04  # electrons/ADU
data_max_2=55000


pxscale_GROWTH = 0.676 # arcsec/pixel 
scale_low_2 = 0.5
scale_high_2 = 1.0
pxscale_PS1 = 0.25 # arcsec/pixel
gain_PS1 = 1.01

#---------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------#

cwd=os.getcwd()
DIR_aligned=cwd+"/SN_ALIGNED/"

#---------------------------------------------------------------------------------------------------------------------#
# function for removing files

def remove_file(file_name):
    try:
        os.remove(file_name)
    except OSError:
        pass
        
# function for removing files having similar names

def remove_similar_files(common_text):
    for residual_file in glob.glob(common_text):
        remove_file(residual_file)

def group_similar_files(text_list, common_text, exceptions=''):
    list_files=glob.glob(common_text)
    if exceptions !='':
        list_files=list(filter(lambda z: not re.search(exceptions, z), list_files))

    list_files.sort()
    if len(text_list) !=0:
        with open(text_list, 'w') as f:
            for file_name in list_files:
                f.write(file_name+'\n')
                
    return list_files    
    


def text_list_to_python_list(text_list):
    if os.path.exists(text_list):
        with open(text_list, 'r+') as f:
            python_list=f.read.split()
            
    return python_list
    
def python_list_to_text_list(python_list, text_list):
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(str(element)+'\n')
            
            
            
def list_lists_to_list(list_lists, text_list):
    list_name=[]
    for file_name in list_lists:
        with open(file_name, 'r') as f:
            file_list=f.read().split()
            for element in file_list:
                list_name.append(element)
    python_list_to_text_list(list_name, text_list)
    
    return list_name
    
class color:

   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

        
    
    
def display_text(text):

    print('#'+'-'*(10+len(text))+'#')
    print('#'+('-'*5)+str(text)+('-'*5)+'#')
    print('#'+'-'*(10+len(text))+'#')        
    
#---------------------------------------------------------------------------------------------------------------------#


for text in ['*.xyls', '*.axy', '*.corr', '*.match', '*.new', '*.wcs', '*.solved', '*.rdls', '*.png', '*.ps', '*.coo', '*.list', "*_resamp.fits", 
             "*_resamp.weight.fits", "log*", "list_*", '*.txt', '*.ldac', 'awcs_*.fits']:

    remove_similar_files(common_text=text)

#iraf.noao(_doprint=0)
#iraf.images(_doprint=0)

# function to edit the header

def copy_header(file_name):
        
    headerlist = fits.open(file_name, mode = 'update')        
    file_header = headerlist[0].header
    date_obs = file_header['DATE-OBS']
    exposure_time = file_header['EXPTIME']
    Filter = file_header['FILTER']
    
    object_ra = RA
    object_dec = DEC
    Object = OBJECT_NAME
        
    #print object_ra
    #print object_dec
    #print Object      
        
    if re.search(':', object_ra):
        c = SkyCoord(object_ra, object_dec, unit=(u.hourangle, u.deg))
        #print c
        ra_deg = c.ra.degree
        dec_deg = c.dec.degree
        RA_degrees = round(ra_deg, 5)
        DEC_degrees = round(dec_deg, 5)
    else:
        print("The RA and DEC are already in degrees") 
        
        RA_degrees == object_ra 
        DEC_degrees == object_dec     
        
                
    list_keywords = ['DATE-OBS', 'EXPTIME', 'FILTER', 'RA', 'DEC', 'OBJECT', 'RA_DEG', 'DEC_DEG'] 
    dict_header = {'DATE-OBS': date_obs, 'EXPTIME': exposure_time, 'FILTER': Filter, 'RA': object_ra, 'DEC': object_dec, \
                 'OBJECT': Object, 'RA_DEG': RA_degrees, 'DEC_DEG': DEC_degrees}
    comment = {'DATE-OBS': 'The date of observation', 'EXPTIME': 'Time of exposure of the frame', 'FILTER': 'Filter', \
                'RA': 'Right Ascension of the target', 'DEC': 'Declination of the target', 'OBJECT': 'Name of the Object under study', \
                'RA_DEG': 'Right Ascension in degrees', 'DEC_DEG': 'Declination in degrees'}
    
    for keyword in list_keywords:
        if keyword in file_header.keys():
            file_header.remove(keyword, remove_all = True)
        file_header.append(card = (keyword, dict_header[keyword], comment[keyword])) 
        
    headerlist.flush()
    headerlist.close()                 
        
                                          
        
# -------------------------------------------------------------------------------------------------------#

# Ignoring warnings in output 

warnings.filterwarnings("ignore")

# -------------------------------------------------------------------------------------------------------#

"""
FUNCTIONS TO CONVERT FITS FILES OR ASTROPY TABLES TO FITS_LDAC FILES AND
VICE VERSA.

"""

def convert_hdu_to_ldac(hdu):
    
    """
    Convert an hdu table to a fits_ldac table
   
    """

    tblhdr = np.array([hdu.header.tostring(',')])
    col1 = fits.Column(name='Field Header Card', array=tblhdr, format='13200A')
    cols = fits.ColDefs([col1])
    tbl1 = fits.BinTableHDU.from_columns(cols)
    tbl1.header['TDIM1'] = '(80, {0})'.format(len(hdu.header))
    tbl1.header['EXTNAME'] = 'LDAC_IMHEAD'
    tbl2 = fits.BinTableHDU(hdu.data)
    tbl2.header['EXTNAME'] = 'LDAC_OBJECTS'
    return (tbl1, tbl2)

def convert_table_to_ldac(tbl):
    
    """
    Convert an astropy table to a fits_ldac
    
    """
    import tempfile
    f = tempfile.NamedTemporaryFile(suffix='.fits', mode='rb+')
    tbl.write(f, format='fits')
    f.seek(0)
    hdulist = fits.open(f, mode='update')
    tbl1, tbl2 = convert_hdu_to_ldac(hdulist[1])
    new_hdulist = [hdulist[0], tbl1, tbl2]
    new_hdulist = fits.HDUList(new_hdulist)
    return new_hdulist

def save_table_as_ldac(tbl, filename, **kwargs):
    
    """
    Save a table as a fits LDAC file
    
    """
    hdulist = convert_table_to_ldac(tbl)
    hdulist.writeto(filename, **kwargs)



def get_table_from_ldac(filename, frame=1):
    
    """
    Load an astropy table from a fits_ldac by frame 
    
    """
    from astropy.table import Table
    if frame>0:
        frame = frame*2
    tbl = Table.read(filename, hdu=frame)
    return tbl

# -------------------------------------------------------------------------------------------------------#

            
# function to run Astrometry.net and find initial estimate of WCS of the images
  
def run_astrometry(ctext, scale_low, scale_high):

    if os.path.exists('astrometryfiles'):
        os.remove('astrometryfiles')

    list_unsolved = []		
            
    file_list='astrometryfiles'
    list_files=group_similar_files(file_list, common_text=ctext, exceptions = 'wcs_ps1_reference')
    for filename in list_files:
        hdul=fits.open(filename)
        RA=hdul[0].header['RA_DEG']
        DEC=hdul[0].header['DEC_DEG']
        OBJECT=hdul[0].header['OBJECT']
        data_image=hdul[0].data
        mean, median, std_dev=sigma_clipped_stats(data_image, sigma=3.0)
        nsigma=20*std_dev
        downsample = [1,2,3,4,5]

        astrometry_command='solve-field'+" "+"--ra "+str(RA)+","+" --dec "+str(DEC)+"," \
                        +" --radius 1.0"+","+" "+"--cpulimit 30"+" --downsample "+str(downsample[1])+ \
                        " "+"--tweak-order 0"+"," \
                        +" --overwrite"+" "+"--resort"+" --new-fits"+" "+"wcs_"+filename+" "+filename
        print("----------Astrometry.net-->Solve-Field is running----------")
        print(astrometry_command)
        os.system(astrometry_command) 
        if os.path.exists('wcs_'+filename):
            print('The file named %s has been solved' % filename)
        else:
            astrometry_command='solve-field'+" "+"--ra "+str(RA)+","+" --dec "+str(DEC)+"," \
                        +" --radius 1.0"+","+" --crpix-center"+" "+"--cpulimit 20"+" --downsample "+ str(downsample[1])+" "+"--tweak-order 0"+"," \
                        +" --overwrite"+" "+"--resort"+" --new-fits"+" "+"wcs_"+filename+" "+filename
            print("----------Astrometry.net-->Solve-Field is running----------")
            print(astrometry_command)
            os.system(astrometry_command) 
            
        if os.path.exists('wcs_'+filename):
            print('The file named %s has been solved' % filename)
        else:
            print("The file named %s has not been solved in the first attempt" % filename)
            astrometry_command='solve-field'+" "+"--ra "+str(RA)+","+" --dec "+str(DEC)+"," \
                        +" --radius 1.0"+","+" --crpix-center"+" "+"--cpulimit 30"+" --downsample "+str(downsample[2])+" "+"--tweak-order 0"+"," \
                        +" --overwrite"+" "+"--resort"+" --new-fits"+" "+"wcs_"+filename+" "+filename
            print("----------Astrometry.net-->Solve-Field is running----------")
            print(astrometry_command)
            os.system(astrometry_command)
            
        if not os.path.exists('wcs_'+filename):
            print("The file named %s has not been solved in the second attempt" % filename)
            astrometry_command='solve-field'+" "+"--ra "+str(RA)+","+" --dec "+str(DEC)+"," \
                        +" --radius 1.0"+","+" --crpix-center"+" "+"--scale-low "+str(scale_low)+" "+"--sigma "+str(std_dev)+"," \
                        +" "+"--scale-high "+str(scale_high)+" "+"--cpulimit 30"+" --downsample "+str(downsample[3])+" "+"--tweak-order 0"+"," \
                        +" --overwrite"+" "+"--resort"+" --new-fits"+" "+"wcs_"+filename+" "+filename
            print("----------Astrometry.net-->Solve-Field is running----------")
            print(astrometry_command)
            os.system(astrometry_command)
        
        if not os.path.exists('wcs_'+filename):
            print('The file %s has not been solved. Please check again' % filename)
            list_unsolved.append(filename) 
        if list_unsolved == []:
            print('Congratulations all your files are solved by Astrometry.net')
        else:
            print('The list of your unsolved files. Have a look at them again. May be you can change some parameters like Downsampling, Tweak Order!!',                          list_unsolved)
            python_list_to_text_list(list_unsolved, 'list_unsolved')
                                                                            

        
        
        
'''
Function to query gaia catalog with a given RA, DEC and a search radius
'''


def make_gaia_catalog(ra, dec, radius_deg, catalog_min_mag, catalog_max_mag, catname):
    job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.panstarrs1_best_neighbour AS pbest, \
        gaiadr2.panstarrs1_original_valid AS ps1 WHERE g.source_id = pbest.source_id AND pbest.original_ext_source_id = \
            ps1.obj_id AND CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND ps1.r_mean_psf_mag > %.2f AND \
                ps1.r_mean_psf_mag < %.2f AND pmra IS NOT NULL AND pmdec IS NOT NULL AND abs(pmdec) > 0 AND abs(pmdec) < 40 AND \
                    abs(pmra)>0 AND abs(pmra) < 40 AND ps1.n_detections > 6 AND pbest.number_of_mates=0 AND \
                        pbest.number_of_neighbours=1;"%(ra, dec, radius_deg, catalog_min_mag, catalog_max_mag), dump_to_file = False)
    
    p = job.get_results()
    
    # convert RA and DEC errors from mas(milli arc second) to degrees

    p['ra_errdeg'] = p['ra_error'] / 3.6e6
    p['dec_errdeg'] = p['dec_error'] / 3.6e6
    p['FLAGS']=0
    
    
    # list_cols = ['astrometric_n_obs_al' , 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al',
    # 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'astrometric_weight_al', 'astrometric_pseudo_colour',
    # 'astrometric_pseudo_colour_error', 'mean_varpi_factor_al', 'astrometric_matched_observations', 'visibility_periods_used', 'astrometric_sigma5d_max', 'frame_rotator_object_type', 
    # 'matched_observations', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', 
    # 'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 
    # 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 
    # 'rv_nb_transits', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'l', 'b', 'ecl_lon', 'ecl_lat', 'priam_flags', 'teff_val', 'teff_percentile_lower', 
    # 'teff_percentile_upper', 'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper', 'e_bp_min_rp_val', 'e_bp_min_rp_percentile_lower', 'e_bp_min_rp_percentile_upper', 
    # 'flame_flags', 'radius_val', 'radius_percentile_lower', 'radius_percentile_upper', 'lum_val', 'lum_percentile_lower', 'lum_percentile_upper', 'gaia_astrometric_params', 
    # 'obj_name', 'obj_id', 'ra_2', 'dec_2', 'ra_error_2', 'dec_error_2', 'epoch_mean', 'zone_id', 'obj_info_flag', 'quality_flag', 'designation', 'phot_variable_flag']
    
    p.remove_columns(['astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al', \
    'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'astrometric_weight_al', 'astrometric_pseudo_colour', \
    'astrometric_pseudo_colour_error', 'mean_varpi_factor_al', 'astrometric_matched_observations', 'visibility_periods_used', 'astrometric_sigma5d_max', 'frame_rotator_object_type', \
    'matched_observations', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', \
    'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', \
    'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', \
    'rv_nb_transits', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'l', 'b', 'ecl_lon', 'ecl_lat', 'priam_flags', 'teff_val', 'teff_percentile_lower', \
    'teff_percentile_upper', 'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper', 'e_bp_min_rp_val', 'e_bp_min_rp_percentile_lower', 'e_bp_min_rp_percentile_upper', \
    'flame_flags', 'radius_val', 'radius_percentile_lower', 'radius_percentile_upper', 'lum_val', 'lum_percentile_lower', 'lum_percentile_upper', 'gaia_astrometric_params', \
    'obj_name', 'obj_id', 'ra_2', 'dec_2', 'ra_error_2', 'dec_error_2', 'epoch_mean', 'zone_id', 'obj_info_flag', 'quality_flag', 'DESIGNATION', 'phot_variable_flag', \
    'datalink_url', 'original_ext_source_id'])
    
    if os.path.exists(catname+'.txt'):
        os.remove(catname+'.txt')
    ascii.write(p, catname+'.txt')
    
    if os.path.exists(catname + '.ldac'):
        os.remove(catname + '.ldac')
    save_table_as_ldac(p, catname + '.ldac')

'''
Run Sextractor on the list of images and create thier catalogs
'''

def run_sextractor(ctext):

    for text in ['*.list']:
        remove_similar_files(common_text=text)

    file_list='ImageList.list'
    list_files=group_similar_files(file_list, common_text=ctext)

    for file_name in list_files:
        command="sex %s -c %s -CATALOG_NAME %s -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME %s -MAG_ZEROPOINT 25.0" % (file_name, config_sex, file_name+'.ldac', param_sex)
        os.system(command)
        print('Executing command: %s\n' % command)
        
        
def file_rename(ctext):
    
    if os.path.exists('scampheadfiles'):
        os.remove('scampheadfiles')

    file_list='scampheadfiles'
    list_files=group_similar_files(file_list, common_text=ctext)
    
    for filename in list_files:
        name=filename[:-10]
        os.rename(filename, name+'.head')
        

def run_scamp(ctext,catname):

    if os.path.exists('scampfilelist'):
         os.remove('scampfilelist')
        
                 
    scamp_list='scampfilelist'
    file_list=group_similar_files(scamp_list, common_text=ctext)
    scamp_command="scamp -c %s @%s -ASTREFCAT_NAME %s" %(config_scamp, scamp_list, catname)
    #scamp_command="scamp -c %s @%s" %(config_scamp, scamp_list)
    os.system(scamp_command)
    print('Executing command: %s\n' % scamp_command)   
    
    file_rename(ctext='*.head')    
        
        
def run_swarp(pxscale, ra, dec, gain, file_name, bkg_sub = 'N'):
    
    output_name = 'a'+file_name
    if os.path.exists(output_name):
        os.remove(output_name)
    swarp_command = 'swarp'+" "+'-c'+" "+str(config_swarp)+" "+"-IMAGEOUT_NAME"+" "+str(output_name)+" "+ \
    '-RESAMPLE'+" "+"Y"+" "+'-SUBTRACT_BACK'+" "+bkg_sub+" "+"-PIXELSCALE_TYPE"+" "+"MANUAL"+ " -PIXEL_SCALE"+" "+ \
    str(pxscale)+" "+"-CENTER"+" "+str(ra)+","+str(dec)+" "+"-GAIN_DEFAULT "+ str(gain)+" "+str(file_name)

    print("----------SwArP is running----------")
    print(swarp_command)
    os.system(swarp_command)
    
    
def rename_reference_wcs(ctext):
    
    if os.path.exists('referencelist'):
        os.remove('referencelist')
        
    reference_files='referencelist'
    list_reference=group_similar_files(reference_files, common_text=ctext)
    
    for file_name in list_reference:
        os.rename(file_name, 'wcs_'+file_name)
        
def rename_reference(ctext):
    
    if os.path.exists('referencelist'):
        os.remove('referencelist')
        
        
    reference_files='referencelist'
    list_reference=group_similar_files(reference_files, common_text=ctext)
    
    for file_name in list_reference:
        os.rename(file_name, file_name[4:])

def panstarrs_query(ra_deg, dec_deg, rad_deg, maxmag=18,
                    maxsources=5000):
    """
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['objID', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000',
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag',
                             'zmag', 'e_zmag',
                             'ymag', 'e_ymag'],
                    column_filters={"gmag":
                                    ("<%f" % maxmag)},
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="II/349/ps1")[0]   


def read_data(image_name):

    '''
    image_name: name of the image

    Returns: data_array, image_header
    '''
    
    image = fits.open(image_name)
    image_data = image[0].data
    image_header = image[0].header
    
    return image_data, image_header


def get_ccd_box(data, header):
    
    '''
    data: image data
    header: image header

    Returns: get minimum ra and dec of the image
    '''
    corners = [[0.], [0.5], [1.]] * np.array(data.T.shape)[::-1]
    print (corners)
    w = WCS(header)

    (ra_min, dec_min), (ra_ctr, dec_ctr), (ra_max, dec_max) = w.all_pix2world(corners, 0.)
    
    #print (ra_min, dec_min, ra_max, dec_max, dec_max - dec_min)
    
    #if ra_min > ra_max:
    #    ra_min, ra_max = ra_max, ra_min
    #if dec_min > dec_max:
    #    dec_min, dec_max = dec_max, dec_min
    #max_size_dec = 0.199
    #if dec_max - dec_min > max_size_dec:
    #    dec_min = dec_ctr - max_size_dec / 2.
    #    dec_max = dec_ctr + max_size_dec / 2.    

    return (ra_ctr, dec_ctr)                                   
        
def check_align(ctext):

    run_sextractor(ctext)
    resampled_catalogs = 'catalog.list'
    list_catalogs = group_similar_files(resampled_catalogs, common_text = 'awcs_*.ldac')
    print (list_catalogs)

    for file_name in list_catalogs:
        sourceTable = get_table_from_ldac(file_name)
        data, header = read_data(file_name[:-5])

        #date_obs = header['DATE-OBS']
        #filt = header['FILTER']

        center_coords = get_ccd_box(data, header)

        print ("The central coordinates are:", center_coords)

    # Query PS1 catalog around the center with a 1.0 deg radius
        ps1_catalog = panstarrs_query(center_coords[0], center_coords[1], 1.0)

    # Read the image coordinates in pixels

        w = WCS(header)
        image_coords = w.all_world2pix(ps1_catalog['RAJ2000'], ps1_catalog['DEJ2000'], 1)

        good_cat_stars = ps1_catalog[np.where((image_coords[0] > 500) & 
                 (image_coords[0] < 3500) & (image_coords[1] > 500) &
                 (image_coords[1] < 3500))]

        cleanSources = sourceTable[(sourceTable['FLAGS']==0) & (sourceTable['FWHM_WORLD'] < 2) & 
                           (sourceTable['XWIN_IMAGE']<3500) & (sourceTable['XWIN_IMAGE']>500) &
                           (sourceTable['YWIN_IMAGE']<3500) &(sourceTable['YWIN_IMAGE']>500) &
                          (sourceTable['ELLIPTICITY'] < 0.2)]

        sourceCatCoords = SkyCoord(ra=cleanSources['X_WORLD'], dec=cleanSources['Y_WORLD'], 
                               frame='icrs', unit='degree')
        ps1CatCoords = SkyCoord(ra=good_cat_stars['RAJ2000'], dec=good_cat_stars['DEJ2000'], 
                            frame='icrs', unit='degree')


         #Now cross match sources

        #Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
        photoDistThresh = 1.0
        idx_image, idx_ps1, d2d, d3d = ps1CatCoords.search_around_sky(sourceCatCoords, 
                                                                  photoDistThresh*u.arcsec)
        #idx_image are indexes into sourceCatCoords for the matched sources, while idx_ps1 are indexes into ps1CatCoords for the matched sources

        print('Found %d good cross-matches within %f pixels'%(len(idx_image), photoDistThresh/0.696))            

        astrometry_offset_ra = np.array(good_cat_stars['RAJ2000'][idx_ps1] - cleanSources['X_WORLD'][idx_image])
        astrometry_offset_dec = np.array(good_cat_stars['DEJ2000'][idx_ps1] - cleanSources['Y_WORLD'][idx_image])
#Compute sigma clipped statistics
        offset_ra_mean, offset_ra_med, offset_ra_std = sigma_clipped_stats(astrometry_offset_ra)
        offset_dec_mean, offset_dec_med, offset_dec_std = sigma_clipped_stats(astrometry_offset_dec)

        #print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))

        print ("Mean offset in ra:", offset_ra_mean)
        print ("Mean offset in dec:", offset_dec_mean)

#rename_reference_wcs(ctext='reference_*.fits')     

display_text('ASTROMETRY IS BEING PERFORMED ON THE IMAGES USING ASTROMETRY.NET')
if Astrometry:
    for text in ['wcs_*.fits']:
        remove_similar_files(text)

    list_astrometry = group_similar_files('', 'cfb*.fits')
    #print list_astrometry
    for file_name in list_astrometry:
        copy_header(file_name)
    #for text in ['wcs_*']:
    #	remove_similar_files(common_text = text)
    if TELESCOPE == 'GIT':
        group_similar_files('list_object', 'cfb*.fits')
        run_astrometry(ctext = '*.fits', scale_low = scale_low_2, scale_high = scale_high_2)
    else:
        group_similar_files('list_object', 'cfb*.fits')
        run_astrometry(ctext = '*.fits', scale_low = scale_low_2, scale_high = scale_high_2)        

list_wcsfiles=group_similar_files('list_wcsobject', 'wcs_*.fits')

img_header=fits.open(list_wcsfiles[0])

data=img_header[0].data
header=img_header[0].header
w=WCS(header)
[ra,dec]=w.all_pix2world(data.shape[0]/2, data.shape[1]/2,1)
print('The Right Ascension of the Center of the field is:', ra)
print('The Declination of the Center of the field is:', dec)

display_text('GENERATING GAIA CATALOG')
print('EXTERNAL GAIA CATALOG IS BEING GENERATED AROUND THE RA %f and DEC %f with a search radius of %d degree' %(ra, dec, 1))
make_gaia_catalog(ra, dec, 0.5, 10, 20, catname='gaiacatalog')

display_text('The GAIA catalog has been generated')
display_text('RUNNING SEXTRACTOR TO GENERATE CATALOG FILES FOR EACH IMAGE')

if make_catalog:
    #for text in ['*.ldac']:
    #	remove_similar_files(text)
    run_sextractor(ctext='wcs_*.fits')


if coarse_align:
    for text in ['awcs*.fits']:
        remove_similar_files(text)
    for file_name in list_wcsfiles:
        if not re.search('wcs_reference', file_name):
            run_swarp(pxscale_HCT, ra, dec, gain_1_new, file_rename)




if Precise_align:
    for text in ['*.head']:
        remove_similar_files(text)
    display_text('RUNNING SCAMP TO SOLVE FOR ASTROMETRY')
    run_scamp(ctext='wcs_*.ldac', catname='gaiacatalog.ldac')	


    display_text('RUNNING SWARP TO ALIGN THE IMAGES')
    for text in ['wcs*.fits']:
        remove_similar_files(text)
    for file_name in list_wcsfiles:
        if not re.search('wcs_reference', file_name):
            run_swarp(pxscale_GROWTH, ra, dec, gain_2, file_name)


    print('Precise align is enabled')

'''
    if TELESCOPE == 'HCT':
        
        
        
        display_text('RUNNING SWARP TO ALIGN THE IMAGES')
        for file_name in list_wcsfiles:
            if not re.search('wcs_reference', file_name):
                run_swarp(pxscale_HCT, ra, dec, gain_1_new, file_name)
            else:
                run_swarp(pxscale_HCT, ra, dec, gain_1_new, file_name)                    
        #elif TELESCOPE =='GROWTH':
       # 	display_text('RUNNING SCAMP TO SOLVE FOR ASTROMETRY')
       # 	run_scamp(ctext='wcs_*.ldac', catname='gaiacatalog.ldac')
       # 	display_text('RUNNING SWARP TO ALIGN THE IMAGES')
       # 	for file_name in list_wcsfiles:
       # 		if not re.search('wcs_ps1_reference', file_name):
       # 			run_swarp(pxscale_GROWTH, ra, dec, gain_2, file_name, bkg_sub = 'Y')
       #     	    	else:               
#		               run_swarp(pxscale_GROWTH, ra, dec, gain_PS1, file_name, bkg_sub = 'N') 

'''               
if align_check:
    check_align(ctext = 'awcs_*.fits') 


list_resampfiles = group_similar_files('', common_text = 'acfb*.fits')
if list_resampfiles != []:
    for file_name in list_resampfiles:
        if not re.search('wcs_reference', file_name):
            print(file_name)
            #copy_header(file_name)          


display_text('ALL THE FILES ARE ALIGNED AND WILL BE COPIED TO SN_ALIGNED')
if os.path.exists(DIR_aligned):
    shutil.rmtree(DIR_aligned)
    
os.mkdir(DIR_aligned)

for file_name in group_similar_files('', common_text='acfb*.fits'):
    if os.path.isfile(DIR_aligned+file_name):
        os.remove(DIR_aligned+file_name)
    shutil.copy(file_name, DIR_aligned+file_name)
    
#for text in ['*.xyls', '*.axy', '*.corr', '*.match', '*.new', '*.wcs', '*.solved', '*.rdls', '*.png', '*_resamp.weight.fits', '*.ps']:
#    remove_similar_files(common_text=text)

display_text('Congratulations !! All the files have been aligned')                  
     
        
        
        
        

    
