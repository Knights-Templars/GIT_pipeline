import os

# Image directory

#================================================================================#

data_dir = '/home/anirban.dutta/SN2022erq_Reduction/SA110/For_Phot/SN_ALIGNED/'

os.chdir(data_dir)

sn_obj = 'sa110_340'
common_text = 'wcs_*.fits'

#================================================================================#

# Import modules and libraries
#================================================================================#

import os
import re
import glob
import warnings
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats

#================================================================================#


config_sex = '/home/anirban/astromatic/config_sextractor.sex'
param_sex = '/home/anirban/astromatic/default.param'

#================================================================================#

# Stop showing warnings. Not a good practice !
warnings.filterwarnings("ignore")

#================================================================================#

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
        list_exceptions = exceptions.split(',')
        for text in list_exceptions:
            list_files=list(filter(lambda z: not re.search(text, z), list_files))
            
    list_files.sort()
    if len(text_list) !=0:
        with open(text_list, 'w') as f:
            for file_name in list_files:
                f.write(file_name+'\n')
                
    return list_files    
    
    
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

def sdss_query(ra_deg, dec_deg, rad_deg, maxmag=18,
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
    vquery = Vizier(columns=['_RAJ2000', '_DEJ2000',
                             'umag', 'e_umag'
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag',
                             'zmag', 'e_zmag'],
                    column_filters={"gmag":
                                    ("<%f" % maxmag)},
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="V/147")[0]                                


def plot_catalogs_stars(catalog, data, header):
    
    '''
    Plot stars in a catalog.
    
    catalog: catalog containing the details of the stars.
    data: image data.
    header: image header.
    
    '''
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca()
    mean, median, std = sigma_clipped_stats(data)
    w = WCS(header)
    plt.imshow(data, vmin=median-3*std, vmax=median+3*std)
    ps1_imcoords = w.all_world2pix(catalog['RAJ2000'], catalog['DEJ2000'], 0)
    
    #print (ps1_imcoords)

    circles = [plt.Circle((ps1_imcoords[0][i], ps1_imcoords[1][i]), 
                      radius = 10, edgecolor='r', facecolor='None') for i in range(len(ps1_imcoords[0]))]
    for c in circles:
        ax.add_artist(c)
    
    plt.show()

#================================================================================#


        
text_list = 'list_files'
list_files = group_similar_files(text_list, common_text=common_text, exceptions='psf,sub')

for index, file_name in enumerate(list_files):

    data, header = read_data(file_name)

    date_obs = header['DATE-OBS']
    filt = header['FILTER']

    print ("Working on %s in filter %s"%(file_name, filt))

    center_coords = get_ccd_box(data, header)

    print ("The central coordinates are:", center_coords)

    # Query PS1 catalog around the center with a 1.0 deg radius
    ps1_catalog = sdss_query(center_coords[0], center_coords[1], 1.0)

    # Read the image coordinates in pixels

    w = WCS(header)
    image_coords = w.all_world2pix(ps1_catalog['_RAJ2000'], ps1_catalog['_DEJ2000'], 1)

    # Select catalog stars putting stars at the edge as outliers 

    good_cat_stars = ps1_catalog[np.where((image_coords[0] > 500) & 
                 (image_coords[0] < 3500) & (image_coords[1] > 500) &
                 (image_coords[1] < 3500))]

    sourceTable = get_table_from_ldac(file_name+'.ldac')                 

    cleanSources = sourceTable[(sourceTable['FLAGS']==0) & (sourceTable['FWHM_WORLD'] < 2) & 
                           (sourceTable['XWIN_IMAGE']<3500) & (sourceTable['XWIN_IMAGE']>500) &
                           (sourceTable['YWIN_IMAGE']<3500) &(sourceTable['YWIN_IMAGE']>500) &
                          (sourceTable['ELLIPTICITY'] < 0.4)]

    # Match the clean sextractor catalog with the ps1 catalog within 3 pixels
    #================================================================================#

    sourceCatCoords = SkyCoord(ra=cleanSources['X_WORLD'], dec=cleanSources['Y_WORLD'], 
                               frame='icrs', unit='degree')
    ps1CatCoords = SkyCoord(ra=good_cat_stars['_RAJ2000'], dec=good_cat_stars['_DEJ2000'], 
                            frame='icrs', unit='degree')

    #Now cross match sources

    #Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
    photoDistThresh = 3.0
    idx_image, idx_ps1, d2d, d3d = ps1CatCoords.search_around_sky(sourceCatCoords, 
                                                                  photoDistThresh*u.arcsec)
    #idx_image are indexes into sourceCatCoords for the matched sources, while idx_ps1 are indexes into ps1CatCoords for the matched sources

    print('Found %d good cross-matches within %f pixels'%(len(idx_image), photoDistThresh/0.696))                          

    psffile_name = file_name + '.als.1'

    psf_table = ascii.read(psffile_name)

    psf_ra, psf_dec = w.all_pix2world(psf_table['XCENTER'], psf_table['YCENTER'], 1)

    psfsourceCatCoords = SkyCoord(ra=psf_ra, dec=psf_dec, frame='icrs', unit='degree')

    photoDistThresh = 5.0
    idx_psfimage, idx_psfps1, d2d, d3d = ps1CatCoords.search_around_sky(psfsourceCatCoords, photoDistThresh*u.arcsec)

    print('Found %d good cross-matches'%len(idx_psfps1))


    psfoffsets = np.array(good_cat_stars[filt+'mag'][idx_psfps1] - psf_table['MAG'][idx_psfimage])
#Compute sigma clipped statistics
    zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
    print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))


    # sn_magfile = file_name + '.als.2'
    # stars_magfile = file_name + '.als.1'

    # sn_magtable = ascii.read(sn_magfile)
    # stars_magtable = ascii.read(stars_magfile)

    # psf_mag = sn_magtable['MAG']
    # psf_err = sn_magtable['MERR']

    # psf_mag_stars = stars_magtable['MAG']
    # psf_mag_stars_err = stars_magtable['MERR']

    # zp = zero_psfmean
    # zp_err = zero_psfstd

    # sn_mag = psf_mag + zp
    # stars_mag = psf_mag_stars + zp

    # sn_err = np.sqrt(psf_err**2 + zp_err**2)
    # stars_mag_err = np.sqrt(psf_mag_stars_err**2 + zp_err**2)

    # file_name_save = 'mag_'+filt+'.txt'

    # if os.path.exists(file_name_save):
    #     os.remove(file_name_save)

    # f = open(file_name_save, 'w')
    # f.write("{0:>4s}{1:>20s}{2:>15s}{3:>15s}\n\n".format('MAG', 'MERR', 'RA','DEC'))

    # for l in zip(stars_mag, stars_mag_err, psf_ra, psf_dec):
    #         f.write("{0:>4f}{1:>20f}{2:>15f}{3:>15f}\n\n".format(*l))
        
        
    # f.close()

    # print('%s   %s  %s  %s  %.2f  %.2f '%(file_name, sn_obj, date_obs, filt, sn_mag, sn_err))