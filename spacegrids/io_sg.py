#encoding:utf-8

""" io related
"""

import numpy as np

from config import *

# use_scientificio is set in config
use_scientificio = False

if use_scientificio is True:  
  try:
    import Scientific.IO.NetCDF
    print 'Using Scientific.IO.NetCDF'
  except:
    print 'no Scientific io. Reverting to scipy'
    from scipy.io import netcdf  
    use_scientificio = False
else:
  from scipy.io import netcdf         

import os

from fieldcls import *

warnings.formatwarning = warning_on_one_line

def netcdf_file(filepath,mode = 'r'):
  """
  Wrapper for open Netcdf functions from ScientificIO or Scipy
  """

  if use_scientificio is True:
    try:
      file = Scientific.IO.NetCDF.NetCDFFile(filename = filepath, mode = mode)
    except IOError:
     raise IOError('Cannot open %s using Scientific.IO'%filepath)
      
    else:
      return file  


  else:
    # Scipy way:

    try:
      file = netcdf.netcdf_file(filename = filepath, mode = mode)
    except IOError:
      raise IOError('Cannot open %s using Scipy'%filepath)
    else:
      return file




def msk_read(filepath='masks/msk', crop = 1):

  """
  Reads a text file containing a mask pointed to by filepath and returns a corresponding array.
  Due to the way these masks are stored for the UVic model, cropping is needed, as indicated by the crop flag in the arguments.
This is the lowest level mask read function in sg.

  """

  str_data = []
  
  with open(filepath,'r') as fobj:
    str_data = fobj.readlines()
  
  data = [] 
  for eachline in str_data:
    data_line = []
    for elem in eachline:
      try:
        data_line.append(int(elem))  
      except:
        pass
    data.append(data_line)
   

  if crop:
    return np.flipud(np.array(data))[1:-1,1:-1]
  else:  
    return np.flipud(np.array(data))

def read_masks(dir_path, msk_shape=0,grids = False, msk_val =2, parent = 'orphan'):
  
      """
      Reads mask and returns a list of field objects containing masks.
      """

      if not(grids):
        print 'read_masks: Provide grid --> no masks loaded.'
        return

      if isinstance(grids,gr):
        grids = [grids]
   
      try:
        L = os.listdir(dir_path)
      except:
        print "No mask dir."
        L=[]      

      masks = {}
      for l in L:
        try:
          fpath = os.path.join(dir_path,l)
	  msk = msk_read(fpath)

          if msk_shape:
# only test shape if msk_shape not default value of 0
            if (msk.shape != tuple(msk_shape)):
              print "Warning: mask shape does not agree: " + l,
              print msk.shape,
              print msk_shape
          msk = np.array(msk,dtype = np.float32)
          
          if msk_val:
            msk[msk != msk_val] = np.nan
            msk[msk == msk_val] = 1 

          for g in grids:
            try:
              print 'Trying to join mask and grid to create field, ignore possibe error.'
              mskob = field(name = l, value = msk, grid =g)
              break
            except:
              pass

          masks[l] = mskob
 
        except:	  
          print "No mask."

      return masks


def locate(top = '/home/',fname = projnickfile):
  """
  Locates all files with filename fname. Helper function to info function.
  
  Inputs: 
  top		(default '/home/') the start dir
  fname		(default projname) the filename to look for. 
  
  Outputs:
  paths		all paths to dirs containing fname
  
  
  """
  paths = []
  for root, dirs, files in os.walk(top=top):
    if fname in files:
      paths.append(root)
  
  return paths







