#!/usr/bin/env python
# coding: utf-8

# In[15]:


# Import the required modules and functions

import os, re, shutil
import glob
from astropy.io import fits
from pyraf import iraf


# In[16]:


#  Functions for Handling files 

cwd = '/home/anirban.dutta/SN2022erq_Reduction/SA113/'
os.chdir(cwd)
OBJECT_STUDY = 'SA113_466'
Spec_standard = 'Feige110'
Phot_standard = 'PG'


def remove_file(filename):
    os.remove(filename)
    
def remove_similar_files(common_text):
    for filename in glob.glob(common_text):
        remove_file(filename)
        
def group_similar_files(text_list, common_text, exceptions = ''):
    list_files = glob.glob(common_text)
    if exceptions != '':
        list_files = filter(lambda x: not re.search(exceptions, x), list_files)
    
    list_files.sort()
    if len(text_list) != 0:
        with open(text_list, 'w') as f:
            for filename in list_files:
                f.write(filename+'\n')
                
    return list_files

for text in ['*.list', 'c_*', 'flat_*', 'bias_*', 'Fe*.fits', '*SN*.fits', '*Gr*', 'OBJECT_*', '*2019-*', '2020-*']:
    remove_similar_files(common_text = text)


# In[17]:


def imcopy(common_text, prefix_str):
    
    task = iraf.images.imutil.imcopy
    task.unlearn()
    
    task.verbose = 'yes'
    task.mode = 'ql'
    
    list_files = 'imcopy.list'
    input_filename = '@'+list_files+'[1]'
    output_filename = prefix_str + '//' + '@' + list_files
    group_similar_files(text_list = list_files, common_text = common_text)
    task(input = input_filename, output = output_filename)

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
    print('#'+('-'*5) +str(text)+ ('-'*5)+'#')
    print('#'+'-'*(10+len(text))+'#')

def hselect(common_text):
    
    task = iraf.images.imutil.hselect
    task.unlearn()
    
    task.fields = '$I,NAXIS1,NAXIS2,OBJECT,FILTER,GRISM,EXPTIME'
    task.expr = 'yes'
    task.missing = 'INDEF'
    task.mode = 'ql'
    
    list_files = 'hselect.list'
    output_filename = 'imageslog.list'
    input_filename = '@' + list_files
    print output_filename
    print input_filename
    group_similar_files(text_list = list_files, common_text = common_text)
    
    #task(images = input_filename, Stdout = output_filename)
    

def file_rename(common_text, logfile = 'log.list'):
    
    list_str = map(lambda x:x.lower(),['DARK','FLAT','Bias_img', 'Bias img', 'BIAS', \
            'bias snspec', 'bias_img', 'bias_snspec', 'eve_flat', 'morning_flat', 'Biasim', 'Bias_imgexit', 
            'Bias_Speco', 'M-fLAT', 'Bias_imaging', 'J0137+149', 'Evening_Flat'])
    list_files = 'images.list'
    list_images = group_similar_files(list_files, common_text)
    with open(logfile, 'w') as f:
        f.write("{0:>4s}{1:>20s}{2:>15s}{3:>15s}{4:>17s}{5:>15s}\n\n".format('FILENAME', 'NAXIS1', 'NAXIS2', 'OBJECT', 'FILTER', 'EXPTIME'))
        for i, filename in enumerate(list_images, 1):
            hdul = fits.open(filename)
            NAXIS1 = hdul[0].header['NAXIS1']
            NAXIS2 = hdul[0].header['NAXIS2']
            OBJECT = hdul[0].header['OBJECT'].lower()
            FILTER = hdul[0].header['FILTER'].replace(' ', '')
            #RISM = hdul[0].header['GRISM'].lower().replace(' ', '')
            EXPTIME = hdul[0].header['EXPTIME']
            DATE_OBS = hdul[0].header['DATE-OBS'].split('T')
            #ATE = hdul[0].header['DATE-OBS']
            DATE = DATE_OBS[0]
            TIME = DATE_OBS[1]
            f.write("{0:>4s}{1:>10d}{2:>11d}{3:>15s}{4:>17s}{5:>15f}{6:>25s}\n".format(filename, NAXIS1,NAXIS2,OBJECT,FILTER,EXPTIME,DATE_OBS))
        
	    
            if EXPTIME < 0.005 and NAXIS1 == 250 and NAXIS2 == 3500:
                newname = DATE+'-'+'bias_snspec_'+str(i)+'.fits'
                os.rename(filename, newname)
                print filename,'-->', newname
            elif EXPTIME < 0.005 and NAXIS1 == 4096 and NAXIS2 == 4108:
                newname = DATE+'-'+'bias_'+str(i)+'.fits'
                os.rename(filename, newname)
                print filename,'-->',newname
            elif  EXPTIME > 0.01 and bool([x for x in list_str if(x in OBJECT)]):
                newname = DATE+'-'+ 'flat_'+FILTER+'_'+str(i)+'.fits'
                os.rename(filename, newname)
                print filename,'-->',newname
                
            elif NAXIS1 == 4096 and NAXIS2 == 4108:
				newname = DATE+'-'+OBJECT+'_'+FILTER+'_'+str(i)+'.fits'
				os.rename(filename, newname)
				print filename,'-->',newname 

            elif Phot_standard.lower() in OBJECT:
				newname = DATE+'-'+OBJECT+'_'+FILTER+'_'+str(i)+'.fits'
				os.rename(filename, newname)
				print filename,'-->',newname


            elif Spec_standard.lower() in OBJECT and NAXIS1 == 250 and NAXIS2 == 3500:
                if GRISM == '4grism7':
                    newname = DATE+'-'+OBJECT+'_Gr7_'+str(i)+'.fits'
                    os.rename(filename, newname)
                    print filename,'-->',newname
                elif GRISM == '3grism8':
                    newname = DATE+'-'+ OBJECT+'_Gr8_'+str(i)+'.fits'
                    os.rename(filename, newname)
                    print filename,'-->',newname

                        
            elif  OBJECT == 'fear':
                newname = DATE+'-'+ 'FeAr_'+str(i)+'.fits'
                os.rename(filename, newname)
                print filename,'-->',newname
            
            elif OBJECT == 'fene':
                newname = DATE+'-'+'FeNe_'+str(i)+'.fits'
                os.rename(filename, newname)
                print filename,'-->',newname
            elif OBJECT  in OBJECT_STUDY.lower():
                if GRISM == '4grism7':
                    newname = DATE+'-'+ OBJECT_STUDY+'_Gr7_'+str(i)+'.fits'
                    os.rename(filename, newname)
                    print filename,'-->',newname
                elif GRISM == '3grism8':
                   newname = DATE+'-'+OBJECT_STUDY+'_Gr8_'+str(i)+'.fits'
                   os.rename(filename, newname)
                   print filename,'-->',newname
            
display_text('Using The Imcopy Task in IRAF')       
#imcopy(common_text = '*.fits', prefix_str = 'c_')
display_text('Renaming files')
file_rename(common_text = '*.fits')




