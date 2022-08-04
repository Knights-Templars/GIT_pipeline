
#--------------------------------------------------------------------------------#

# This script will be used for creating a master-bias frame

# preprocess_bs_git_v1.py 

#--------------------------------------------------------------------------------#

# User can change this

cwd = '/home/anirban.dutta/SN2022erq_Reduction/bias/'
OBJECT_NAME='SN2022erq'
RA_object = '18:33:25.360'
DEC_object = '+44:05:11.65'
Standard_Name = 'SA'

#--------------------------------------------------------------------------------#


# Starting of Code
# Import required libraries

import os
import re
import sys
import glob
import shutil
import numpy as np
from pyraf import iraf
from astropy.io import fits

#--------------------------------------------------------------------------------#

Software_Code = 'preprocess_v1'

#--------------------------------------------------------------------------------#

# Global Variables

filters=['u', 'g', 'r', 'i', 'z']

# GIT New CCD specifications

read_noise=12.0
gain=1.04
data_max=55000

# HCT CCD specifications

read_noise_new = 5.75
gain_new = 0.28
data_max_new = 700000

#--------------------------------------------------------------------------------#

#image header keyword

Date='DATE-OBS'
Filter='FILTER'
OBJECT='OBJECT'
CCDSECTION = 'CCDSEC'
RA = 'RA'
DEC = 'DEC'

# Working Directory 


DIR_PHOT = 'For_Phot'

print('The current working directory is:\n %s' % cwd)

os.chdir(cwd)

# Functions for file handling

# function for removing files: argument: file_name

def remove_file(file_name):
	
	try:
		os.remove(file_name)
	except OSError:
	    pass

# Function for removing files having similar names: argument: common_text

def remove_similar_files(common_text):
	
	for residual_file in glob.glob(common_text):
		remove_file(residual_file)
		
		
for text in ['*b_flat*', 'list_*']:
    remove_similar_files(common_text=text)

# Function to group similar files: arguments: text_list, common_text, exceptions

def group_similar_files(text_list, common_text, exceptions=''):
	list_files=glob.glob(common_text)
	if exceptions !='':
		list_files = filter(lambda x: not re.search(exceptions, x), list_files)
	list_files.sort()
	if len(text_list) !=0:
		with open(text_list, 'w') as f:
			for file_name	in list_files:
				f.write(file_name + '\n')

	return list_files


def text_list_to_python_list(text_list):
	if os.path.exists(text_list):
		with open(text_list, 'r+') as f:
			python_list=f.read().split()

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
    print('#'+('-'*5) + str(text) +('-'*5)+'#')
    print('#'+'-'*(10+len(text))+'#')
            
            

# Load IRAF Packages:

iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.crutil(_doprint=0)
iraf.images(_doprint=0)
iraf.ccdred.instrument='ccddb$kpno/direct.dat'

def edit_header(textlist_files, object_name):
    
    list_files = text_list_to_python_list(textlist_files)
    
    for file_name in list_files:
        hdulist = fits.open(file_name, mode = 'update')
        file_header = hdulist[0].header
        OBJECT = file_header['OBJECT']
        #RA = file_header['RA']
        #DEC = file_header['DEC']
        
        
        
        list_keywords = ['OBJECT', 'CCDSEC', 'RA', 'DEC']
        dict_header = {'OBJECT': OBJECT, 'CCDSEC': CCDSECTION, 'RA': RA, 
        			   'DEC': DEC}
       
        for keyword in list_keywords:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all = True)

        list_update_keywords = ['OBJECT', 'RA', 'DEC', 'COMMENT', 'RN', 'GAIN']
        dict_update_header = {'OBJECT': object_name, 'RA': RA_object, 
                              'DEC': DEC_object, 'COMMENT': Software_Code, 'RN': read_noise,
                              'GAIN': gain}

        for keyword in list_update_keywords:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all = True)
            file_header.append(card=(keyword, dict_update_header[keyword]))
        hdulist.flush()
        hdulist.close() 

        
            
def zero_combine(textlist_bias, master_bias, rd_noise, ccd_gain):

    '''
    About: combining bias frames


    textlist_bias: list of bias as a text file.
    master_bias: name of output master bias file (.fits)
    rd_noise: provide readnoise. check ccd_property script for finding readnoise and gain.
    ccd_gain: provide ccd_gain. check ccd_property script for finding readnoise and gain.
    '''
    
    
    task=iraf.noao.imred.ccdred.zerocombine
    task.unlearn()
    
    task.combine='median'
    task.reject='ccdclip'
    task.ccdtype=''
    task.process='no'
    task.delete='no'
    task.clobber='no'
    task.scale='none'
    task.statsec=''
    task.nlow=0
    task.nhigh=1
    task.nkeep=1
    task.mclip='yes'
    task.lsigma=3.
    task.hsigma=3.
    task.snoise=0
    task.pclip=-0.5
    task.blank=0
    task.mode='ql'
    
    #remove_file(master_bias)
    task(input='@'+ textlist_bias, output = master_bias, rdnoise=rd_noise,
        gain=ccd_gain)

def bias_subtract(textlist_tbs, master_bias='mbias.fits', prefix_str='b_'):

    task=iraf.noao.imred.ccdred.ccdproc
    task.unlearn()
    
    task.ccdtype=''
    task.max_cac=0
    task.noproc='no'
    task.fixpix='no'
    task.oversca='no'
    task.trim='no'
    task.zerocor='yes'
    task.darkcor='no'
    task.flatcor='no'
    task.illumco='no'
    task.fringec='no'
    task.readcor='no'
    task.scancor='no'
    task.readaxi='line'
    task.fixfile=''
    task.biassec=''
    task.trimsec=''
    task.zero=master_bias
    task.dark=''
    task.flat=''
    task.fringe=''
    task.minrepl=1
    task.scantyp='shortscan'
    task.nscan=1
    task.interac='no'
    task.functio='legendre'
    task.order=1
    task.sample='*'
    task.naverag=1
    task.niterat=1
    task.low_rej=3
    task.high_re=3
    task.grow=0
    task.mode='ql'
    
    output_filename=prefix_str+'//'+'@'+textlist_tbs
    remove_file(output_filename)
    
    task(images='@'+textlist_tbs, output=output_filename)
    
def image_statistics(ctext, text_list='list_bsflat'):

    textlist_bsflat='list_bsflat'
    group_similar_files(textlist_bsflat, common_text=ctext)
    
    task=iraf.images.imutil.imstatistics
    task.unlearn()
    
    task.lower='INDEF'
    task.upper='INDEF'
    task.nclip=0
    task.lsigma=3
    task.usigma=3
    task.binwidt=0.1
    task.cache='no'
    task.mode='ql'
    
    task(format='no', images='@'+textlist_bsflat, field='mean', 
        Stdout=text_list+'_mean')
    task(format='yes', images='@'+textlist_bsflat, 
         field='image,npix,mean,midpt,stddev,min,max', Stdout=text_list+'_stat')
    
    list_mean=text_list_to_python_list(text_list+'_mean')
    #print list_mean
    
    return list_mean
    
def flat_normalize(ctext, textlist_mean='list_bsflat', prefix_str='n'):

    list_bsflat=group_similar_files('', common_text=ctext)
    #print list_bsflat
    #print len(list_bsflat)
    list_bsflat_mean=image_statistics(ctext, textlist_mean)
    
    task=iraf.images.imutil.imarith
    task.unlearn()
    print('FILENAME'+'     '+'MEAN_VALUE')
    for index in range(0, len(list_bsflat)):
        print list_bsflat[index],'    ',list_bsflat_mean[index]
        output_filename=prefix_str+list_bsflat[index]
        remove_file(output_filename)
        task(operand1=list_bsflat[index], op='/', 
             operand2=float(list_bsflat_mean[index]), result=output_filename)
        

def flat_combine(ctext):
    
    band=ctext[4:9]
    #band=ctext[4]
    #print ("Band is:", band)
    textlist_nbflat='list_nbflat'
    group_similar_files(textlist_nbflat, common_text=ctext)
    
    task=iraf.noao.imred.ccdred.flatcombine
    task.unlearn()
    
    task.combine='median'
    task.reject='ccdclip'
    task.ccdtype=''
    task.process='no'
    task.subsets='no'
    task.delete='no'
    task.clobber='no'
    task.scale='mode'
    task.statsec=''
    task.nlow=1
    task.nhigh=1
    task.nkeep=1
    task.mclip='yes'
    task.lsigma=3
    task.hsigma=3
    task.rdnoise=float(read_noise_new)
    task.gain=float(gain_new)
    task.snoise=0
    task.pclip=-0.5
    task.blank=1
    task.mode='ql'
    
    output_filename='mflat'+band+'.fits'
    remove_file(output_filename)
    task(input='@'+textlist_nbflat, output=output_filename)

def flat_correction(ctext, master_flat, exception='flat', prefix_str='f'):
    
    list_bsobject=group_similar_files('',common_text=ctext, exceptions=exception)
    print list_bsobject
    
    task=iraf.images.imutil.imarith
    task.unlearn()
    
    for image in list_bsobject:
        output_filename=prefix_str+image
        remove_file(output_filename)
        task(operand1=image, op='/', operand2=master_flat, result=output_filename)
        
        
def crmedian(ctext, prefix_str='c'):
    
    list_cosmic=group_similar_files('', common_text=ctext)
    
    task=iraf.noao.imred.crutil.crmedian
    task.unlearn()
    
    task.crmask=''
    task.median=''
    task.sigma=''
    task.residua=''
    task.var0=0
    task.var1=0
    task.var2=0
    task.lsigma=10
    task.hsigma=3
    task.ncmed=5
    task.nlmed=5
    task.ncsig=25
    task.nlsig=25
    task.mode='ql'
    
    for image in list_cosmic:
        output_filename=prefix_str+image
        remove_file(output_filename)
        task(input=image, output=output_filename)
        
def astroscrappy(file_list, gain, rd_noise, data_max, sig_clip=5.5):
    cosmics_no = []
    mask_list = []
    for i in range(len(file_list)):
        hdu = fits.open(file_list[i])[0]
        data = hdu.data
        header = hdu.header
        object_ = header['OBJECT']
        crmask, cleanarray = cr.detect_cosmics(data, sigclip = sig_clip, gain = gain,
                             readnoise = rd_noise, satlevel = data_max, sepmed=False,
                            cleantype='medmask', fsmode='median')
        print('Number of affected pixels is %d for file %s'%(np.sum(crmask), 
                                                             file_list[i])) 
        cosmics_no.append(np.sum(crmask))
        
        proc_cr = cleanarray/gain
        clean_hdu = fits.PrimaryHDU(proc_cr)
        clean_hdu.header = header
        clean_hdu.writeto('c'+file_list[i], overwrite = True)
        
        mask_list.append(crmask)
        
    return [cosmics_no, mask_list]        
    

#--------------------------------------------------------------------------------#
# Do you want to remove files from previous run    

remove_prev_files=True

if remove_prev_files:
    for text in ['list_', 'mbias*']:
        remove_similar_files(common_text=text)

    
#--------------------------------------------------------------------------------#        

# list of bias frames
list_bias = group_similar_files('list_bias', '*bias*.fits', exceptions = 'bias_snspec')


# Edit header 
edit_header('list_bias', 'BIAS')


display_text('ZEROCOMBINE TASK IS BEING PERFORMED TO MAKE A MASTER BIAS')


# Run zerocombine task of IRAF
zero_combine(textlist_bias = 'list_bias', master_bias = 'mbias.fits', rd_noise=read_noise, ccd_gain=gain)

display_text('MASTER BIAS FRAME is CREATED. CONGRATS')

#--------------------------------------------------------------------------------#

#------------------------------------THE END-------------------------------------#
