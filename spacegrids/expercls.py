#encoding:utf-8

""" exper class
"""

import types
import datetime
import glob
import warnings

from config import *

# use_scientificio is set in config
#if use_scientificio is True:  
#  try:
#    import Scientific.IO.NetCDF
#    print 'Using Scientific.IO.NetCDF'
#  except:
#    print 'no Scientific io. Reverting to scipy'
#    from scipy.io import netcdf  
#    use_scientificio = False
#else:
#  from scipy.io import netcdf         

import os
import copy

from fieldcls import *
from io_sg import *

warnings.formatwarning = warning_on_one_line

# ---------------- Class definition for experiments -----------

class exper:
  
  def __repr__(self):
    return self.name

  def show(self):

    loaded_fields = self.vars.keys()
    l = len(loaded_fields)


    REP = report()
    REP.echoln(self.name)
    REP.line()
    if loaded_fields:
      REP.echoln(loaded_fields, width = 30, cols = 3)

    if len(loaded_fields) != 1:
      plural = 's'
    else:
      plural = ''

    REP.echo('exper using %1.2f Mb. %i field%s loaded.  '%(self.nbytes/1024./1024.  ,len(loaded_fields)  , plural ) )

    print REP.value
    

  def __init__(self,path=home_path, name = 'test',cstack = [], params = {}, descr = 0, parent = 'orphan'):
# --> belongs to class exper

    self.path = os.path.join(path,name) 
    self.name = name
    self.cstack = cstack
    self.parent = parent
      
    if not(descr):
      self.descr = name
    else:
      self.descr = descr  

    # this is needed because otherwise all exper objects have the same params object!
    self.params = copy.deepcopy(params)
    # check for other places where this effect might occur.

          
# The variables associated with this experiment. It is a dictionary of name vs struct    
    self.vars = {}

    self.nbytes = reduce( lambda x,y: x + y, [e.nbytes for e in cstack]  )
    

  def __call__(self, varnames ='*'):

# --> belongs to class expers
    return self.get(varnames)

  def __setitem__(self,varname, arr):

    self.vars[varname] = arr


  def __getitem__(self, varnames ='*'):

# --> belongs to class expers    
    return self.get(varnames)



  def get(self, varnames):
  
    fields = []
    VN = sublist(self.vars.keys(), varnames)
    if VN:
      for vn in VN:
	fields.append(self.vars[vn])
	
    if not(fields):
      try:
        fields = [field(varnames, parent = self)]
      except:
        print "No fields found for %s. Try P.load(\'%s\')"%(self.name, varnames)
        return 
    
    if len(fields) == 1:
# if only single field found, return that field and not a list of fields.
      fields =fields[0]

    return fields    	

  def list_vars(self):
    RP = report()
    if len(self.vars) == 0:
      RP.echo('No variables.')
    else:

      RP.echoln('name  \t descr') 
     
      RP.echoln('-'*20)
      for k in self.vars.keys():
        RP.echo(k,' \t ')
        if hasattr(self.vars[k],'descr'):
          RP.echoln(self.vars[k].descr)
        else:
          RP.echoln()

    return RP      


  def __delitem__(self,i):
  
    del self.vars[i]

  def delvar(self, varnames, msg = ''):
# this delvar is a method of class exper
  
    if not(isinstance(varnames, types.ListType)):
      varnames = [ varnames ]
  
    for varname in varnames:
      if varname in self.vars:
        del self.vars[varname]
        print 'Deleting existing '+varname + ' ... ',  
      else:
        if msg:
          print 'var ' + varname +' not in list. ' + msg,


    self.update_nbytes()

#  def cdf(self):
# --> method of class exper        
#    print 'Netcdf variables available (loaded marked **):'
#    print_box(mark_sublist(self.var_names,self.vars.keys()))
      


  def write(self, path = None, name = None , history = 'Created from Spacegrids ' , insert_dual = True ):

    """
    Write method of exper class.

    Creates Netcdf file and writes all loaded field to it, along with their coord objects.

    """

    if name is None:
      name = self.name

    if not name.split('.')[-1] in ['nc','cdf']:
      name = name +'.nc'
    if not path is None:
      name = os.path.join( path , name ) 
   
    
    if len(self.vars) > 0:
      print 'Writing experiment %s to file %s'%(self.name, name)

      try:
        # using io_sg version of netcdf_file. May use ScientificIO or Scipy
        file_handle = netcdf_file(name , 'w')

      except IOError:
        warnings.warn('Cannot open%s'%name)
      else:

        for fld in self.vars.values():

          file_handle = fld.cdf_insert(file_handle, insert_dual = insert_dual)

        file_handle.history = history + '%s'%str(datetime.datetime.now())
 
#    var_cdf.units = self.units

        file_handle.close()

    else:
      print 'No fields to write for experiment %s'%self.name
  


      
  def load(self,varnames, squeeze_field = True, ax=None, name_suffix='_cat', new_coord_name = 'gamma', new_coord= None ):
    """
    Field load method of exper class.

    Load a variable or list of variables contained in varnames. Takes either a single string or a list of strings. If multiple files inside a directory contain the same variable, this method will attempt to concatenate them (e.g. in the case where there are different time slices).

    Inputs:
    varnames		list of the variable names to load
    squeeze_field	switch to squeeze field on loading (default True)	

    The following arguments are passed on to the concatenate function:
    ax			
    name_suffix
    new_coord_name
    new_coord
 
    if self.path is to a file (likely to be an experiment file), the variable will be loaded from that file.
    if self.path is to a directory (likely to be an experiment dir), the variable will be loaded from Netcdf files inside that directory.

    """  

# --> this load is a method of class exper

    if not(isinstance(varnames, list)):
      varnames = [varnames ]
    
    for varname in varnames: 
# Big loop.

     
      # Prepare the paths to all the netcdf files into a list     
      if os.path.isfile(self.path):
        # this exper object was created from a project containing an experiment file
        paths = [self.path]

      else:            
          # experiment corresponds to directory containing Netcdf files. Construct list of full paths to those files
          paths = []
          for root, dirs, files in os.walk(top = self.path):
            for fname in files:
              if fname.split('.')[-1] in ['nc','cdf']:
                paths.append(os.path.join(root,fname))


      # test if var already in list.
      self.delvar(varname, msg = "")  	
            
#	Try to find netcdf var in any of the found files, and then read into field object.

      fnames = []
      F = []
      for filepath in paths:

        try:       
          file = netcdf_file(filepath,'r')
          fnm = os.path.split(filepath)[1]
        except IOError:
          raise IOError('Cannot open %s'%filepath)
        else:

          if varname in file.variables:
           
            F.append(cdfread(filepath,varname,self.cstack,self.axes))
            fnames.append(fnm)
#          break
          file.close()

     
      if F == []:
     
        print '%s for %s could not be read.'%(varname,self.name) 

#        self.vars[varname] = None
      else:
        num_files = len(paths)
        # determine whether to use plural of word 'file' in stdout message.
        if num_files > 1:
          plural = 's'
        else:
          plural = ''
        print 'OK. Fetched field %s for %s. %s file%s found.'%(varname, self.name,num_files,plural)

        new_field = concatenate(F,ax=ax, name_suffix=name_suffix, new_coord_name = new_coord_name, strings = fnames, new_coord= new_coord   )

        # insert field into experiment
        if squeeze_field:
          self.insert(what = (varname, squeeze( new_field ) ) ) 
        else:
          self.insert(what = (varname, new_field ) ) 

    self.update_nbytes() 


  def insert(self,what):
    """
    Insert field in field list of this exper object under key varname.  If argument what is a number, it will be inserted as a parameter.

    Input: what a 2 tuple (pair) of name and value: (name, value). Value can be a field or a single value. Name must be a string, but can be None in the case of a field, where the field name will then be used. For example what = ('temp',TEMP), where TEMP is a field. If value is a single value (e.g. int or float), a name must be provided.

    Argument what can also be a list of (name,value) pairs, in which case the entire collection of pairs will be inserted.   

    """

    if what is None:
      return

    if isinstance(what, list):
      for pair in what:
        self.insert(pair)

    else:

      name = what[0]
      value = what[1]

      if isinstance(value,field):  

        if name is None: 
          name = value.name
        else:
          value = value.copy(name = name)
        # insert field 
        self.vars[name] = value

      else:
      # it is assumed field is a parameter
        if name is None:
          raise Exception('Provide name if trying to insert what as parameter.')

        self.params[name] = value
        


    


  def available(self):
    """
    Method of exper class that returns a list of all available Netcdf variable names (strings).
    """
  
    if os.path.isfile(self.path):
        # this exper object corresponds to a file (not directory).
      paths = [self.path]

    else:            
        # this exper object corresponds to a directory (not file).
        paths = []
        for root, dirs, files in os.walk(top = self.path):
          for fname in files:
            if fname.split('.')[-1] in ['nc','cdf']:
              paths.append(os.path.join(root,fname))

           
# Try to find netcdf var in any of the found files, and then read into field object.

    variables = []
    for filepath in paths:
      try:
        f = netcdf_file(filepath,'r')
      except IOError:
        print 'Cannot open',filepath
      else:
        for var_in_cdf in f.variables.keys():
          if not var_in_cdf in variables:
            variables.append(var_in_cdf)

        f.close()
    return variables


  def ls(self, width = 20):
    """
    Variable list method of exper class.
    Examine which fields (Netcdf variables) are available of experiment object.

    """  
 
    RP = report()
    av_list = self.available()
    av_list.sort()
    RP.echo(av_list,delim = '\t ',maxlen=300, width = width)

    print RP.value


  def update_nbytes(self):
     coord_bytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.cstack]  )
     if len(self.vars) > 0:
       field_bytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.vars.values()]  )    
     else:
       field_bytes = 0
     self.nbytes  = coord_bytes + field_bytes

     
    

# ----- end exper class ------------






def isexpdir(path, file_extensions = cdf_file_extensions):

  """
 
  Tests whether the subdirectories in the path contain data files recognised as known data files and returns the list of those subdirectories (relative path to path argument) that contain these known files. To be used by adexp functionality and such.
file_extensions is the list of known filenames in the form of glob expressions, e.g. ['*.nc','*.cdf'] (the default). 

  
  """

    # examine all subdirectories of path. Create copy for manipulation.

 
  if os.path.isdir(path):
    L = os.listdir(path)
    Lc = copy.deepcopy(L)
  elif os.path.isfile(path):
    return [path,]
  else:
    raise Exception('No such file or directory %s.'%path)
    
  # go through list L of subdirectories of path (e.g. /home/me/PROJECTS/test_project/) to see which ones are experiment directories (i.e. contain .nc and .cdf files).
  for l in L:
    try:
      # look for files ending in .nc or .cdf
      
      exp_path = os.path.join(path , l)
      # Try globbing <path>.nc and <path>.cdf. E.g. /home/me/PROJECTS/test_project/DPC/*.nc for a globfpath
      globfpaths = [os.path.join(exp_path , e) for e in file_extensions]
 
      # Test whether any files of the extensions (e.g .nc) in file_extensions occur in this directory:
      if not(reduce(lambda x,y: x + y, [glob.glob(e) for e in globfpaths])):	
         # if no files with these extensions are found, delete the directory from the list:
        Lc.remove(l)
    except:
     # bad directory anyway
      Lc.remove(l)

  globfpaths = [os.path.join(path , e) for e in file_extensions]
  raw_files = reduce(lambda x,y: x + y, [glob.glob(e) for e in globfpaths])
  files = [e for e in raw_files if os.path.isfile(e)]

  # Returning list of directory and file names (not full paths) believed to correspond to experiments.
  return Lc + [ os.path.split(e)[-1] for e in   files ]

