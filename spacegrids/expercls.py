#encoding:utf-8

"""The Exper class and associated functions. Exper represents experiment data sets. 

An Exper object corresponds to (collections of) Netcdf File(s).

The general workflow starts with the creation of a Project object (providing a path to a project directory). This will lead sg to look through the project directory for subdirectories and Netcdf files (based on the file suffix). An Exper object is created for each subdir and Netcdf file found, and added to the project. So an Exper object may represent either a subdirectory of a project directory or a specific Netcdf file inside that project directory (recorded in the 'path' attribute). Each Exper with name attribute 'foo' can then be accessed via P['foo']. This procedure provides groundwork by interpreting axis and coordinate date found in the Netcdf files (inside the subdirectories in the case where the Exper object represents a subdirectory) and adding this information to the Exper objects (see Attributes). Once this structure is established, specific data sets can be loaded across experiments using P.load.

  Examples:

  >>> import spacegrids as sg
  >>> D = sg.info_dict()
  >>> P = sg.Project(D['my_project'])
  >>> E = P['DPO']
  >>> E.show()
  DPO
  ----------
  Exper using 0.01 Mb. 0 fields loaded.  

"""

import types
import datetime
import glob
import warnings

from _config import *

import os
import copy

from fieldcls import *
from _iosg import *

warnings.formatwarning = warning_on_one_line

# ---------------- Class definition for experiments -----------

class Exper(object):
  """
  Represents experiment data sets. Corresponds to (collections of) Netcdf File(s).

  Attributes:
    axes: (list of Ax objects) link to axes list belonging to project (same list for all Exper objects). Generally constructed from Netcdf data
    coords: (dictionary of name vs Coord objects) for easy reference, see cstack
    cstack: (list of Coord objectd) the coordinate stack of all Coord objects constructed for this experiment
    descr: (str) an optional description of this experiment.
    name: (str) the name of this experiment. Generally derived from data file name
    nbytes: (int) approximate memory usage of experiment in bytes
    params: (dictionary of name value pairs) collection of single value named parameters (e.g. co2 vs 280)
    path: (str) full path on filesystem to experiment directory or file 
    vars: (dictionary of name vs Field objects) contains the loaded Fields (e.g. a 3D dataset of temperature)  


  The general workflow starts with the creation of a Project object (providing a path to a project directory). This will lead sg to look through the project directory for subdirectories and Netcdf files (based on the file suffix). An Exper object is created for each subdir and Netcdf file found, and added to the project. So an Exper object may represent either a subdirectory of a project directory or a specific Netcdf file inside that project directory (recorded in the 'path' attribute). Each Exper with name attribute 'foo' can then be accessed via P['foo']. This procedure provides groundwork by interpreting axis and coordinate date found in the Netcdf files (inside the subdirectories in the case where the Exper object represents a subdirectory) and adding this information to the Exper objects (see Attributes). Once this structure is established, specific data sets can be loaded across experiments using P.load.

  Examples:

  >>> import spacegrids as sg
  >>> D = sg.info_dict()
  >>> P = sg.Project(D['my_project'])
  >>> E = P['DPO']
  >>> E.show()
  DPO
  ----------
  Exper using 0.01 Mb. 0 fields loaded.  
  """
   
  def __repr__(self):
    return self.name

  def show(self):
    """Shows summary of key Exper specifics.

    Shows the experiment name, its memory usage and the number of Field objects loaded.
    """
    loaded_fields = self.vars.keys()
    length = len(loaded_fields)


    REP = Report()
    REP.echoln(self.name)
    REP.line()
    if loaded_fields:
      REP.echoln(loaded_fields, width = 30, cols = 3)

    REP.echo('Exper using %1.2f Mb. %i field%s loaded.  '%(self.nbytes/1024./1024.  ,length  , plural(length) ) )

    print REP.value
    
# --> belongs to  Exper
  def __init__(self,path=home_path, name = 'test',cstack = [], params = {}, vars = {}, descr = None):
    """
    Initialize Exper instance.

    Args:
      path: (str) path to dir containing experiment. Used to construct Exper dir.
      name: (str) the name of this experiment. Generally derived from data file name
      cstack: (list of Coord objectd) the coordinate stack of all Coord objects constructed for this experiment
      params: (dictionary of name value pairs) collection of single value named parameters (e.g. co2 vs 280)
      vars: (dictionary of name vs Field objects) contains the loaded Fields (e.g. a 3D dataset of temperature). Generally empty on init.
      descr: (str) an optional description of this experiment.    
    """

    self.path = os.path.join(path,name) 
    self.name = name
    self.cstack = cstack
      
    if descr is None:
      self.descr = name
    else:
      self.descr = descr  

    # this is needed because otherwise all Exper objects have the same params object!
    self.params = copy.deepcopy(params)
    # check for other places where this effect might occur.

          
# The variables associated with this experiment. It is a dictionary of name vs struct    
    self.vars = copy.deepcopy(vars)

    self.nbytes = reduce( lambda x,y: x + y, [e.nbytes for e in cstack]  )
    

  def __setitem__(self,varname, arr):
    """Simple insertion of name vs Field into vars dictionary. Generally use insert method.
    """
    self.vars[varname] = arr

# --> belongs to Exper 
  def __getitem__(self, varnames ='*'):
    """Obtain Field from vars dictionary. Shortcut to get method. See get.
    """

    return self.get(varnames)



  def get(self, varnames):
    """
    Args:
      varnames: (str or list) simple filter (or list of -) used by fnmatch.fnmatch 

    Returns:
      None if no match, or one Field if only one Field name matches varnames pattern, otherwise list of Field objects. 

    Examples:
    >>> E = P['DPO']
    >>> P.load(['O_temp','O_sal'] )  
    >>> E['O_*'] 
    [O_temp, O_sal]
    """


    # fields is the list of Fields to be returned, whose names match the pattern
    fields = []
    # The names of the matching fields:
    VN = sublist(self.vars.keys(), varnames)
    if VN:
      for vn in VN:
	fields.append(self.vars[vn])
	
    if not(fields):
      print "No fields found for %s. Try P.load(\'%s\')"%(self.name, varnames)
      return 
    
    if len(fields) == 1:
# if only single Field found, return that Field and not a list of fields.
      fields =fields[0]

    return fields    	

  def list_vars(self):
    """Show report of loaded variables in stdout.
    """

    RP = Report()
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
    """
    Delete an item from the vars Field dict attribute of Exper using del.

    Args:
      i: (int) index of item to be deleted.

    Returns:
      None
    """
    del self.vars[i]

  def delvar(self, varnames, msg = ''):
    """
    Delete items from the vars Field dict attribute of Exper using del.

    Args:
      varnames: (str or list) items to delete

    Returns:
      None
    """
  
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


  def write(self, path = None, name = None , history = 'Created from Spacegrids ' , insert_dual = True , *args, **kwargs):
    """Write Exper to Netcdf file.

    *args, **kwargs are passed on to netcdf io function (Dataset() ).

    Args:
      path: (str) path to the directory where the file will go (e.g. 'data/' or '/home/me/', default pwd).
      name: (str) file name (e.g. "foo.nc")
      history: (str) Brief history or general description of the data.
      insert_dual: (Boolean) Flag determining whether to include the duals of the Coord objects in the file.


    Returns:
      None

    Creates Netcdf file and writes all loaded Field to it, along with their Coord objects (and their duals if requested).

    Examples:
    >>> E = P['DPO']
    >>> E.write()  # yields DPO.nc in pwd

    >>> E = P['DPO']
    >>> E.write(path='TMP/',name='foo.nc')  # yields TMP/foo.nc with respect to pwd
    """

    if name is None:
      name = self.name

    if not name.split('.')[-1] in ['nc','cdf']:
      name = name +'.nc'
    if not path is None:
      name = os.path.join( path , name ) 
   
    
    if len(self.vars) > 0:
      print 'Attempting to write experiment %s to file %s'%(self.name, name),

      try:
        # using io_sg version of netcdf_file. May use ScientificIO or Scipy
        file_handle = netcdf_file(name , 'w', *args, **kwargs)

      except (IOError, RuntimeError):
        print '... FAIL'
        warnings.warn('Cannot open %s. File not written.'%name)
      else:

        for fld in self.vars.values():

          file_handle = fld._cdf_insert(file_handle, insert_dual = insert_dual)

        file_handle.history = history + '%s'%str(datetime.datetime.now())
 
#    var_cdf.units = self.units

        file_handle.close()
        print '... OK'

    else:
      print 'No fields to write for experiment %s'%self.name
  


      
  def load(self,varnames, squeeze_field = True, ax=None, name_suffix='_cat', new_coord_name = 'gamma', new_coord= None,slices = None, slice_suffix = '_sliced' , *args, **kwargs):
    """
    Load a variable or list of variables contained in varnames into Exper. 

    Takes either a single string or a list of strings. If multiple files inside a directory contain the same variable, this method will attempt to concatenate them (e.g. in the case where there are different time slices).

    if self.path is to a file (an Experiment file), the variable will be loaded from that file.
    if self.path is to a directory (an experiment dir), the variable will be loaded from Netcdf files inside that directory.

    *args, **kwargs are passed on to the 'netcdf_file' function that handles opening of the file. 

    Args:
      varnames: (str or list) var name or list of the var names to load
      squeeze_field (boolean):	Flag to squeeze Field on loading (default True)	
      ax (Ax): passed on to the concatenate function			
      name_suffix (str): passed on to the concatenate function
      new_coord_name (str): passed on to the concatenate function
      new_coord (Coord): passed on to the concatenate function
      slices: (tuple of slice, Coord and Ax objects) slices to take. No slicing if None.     
      slice_suffix: (str) suffix to add to variable name in case of slicing


    Returns:
      None

    Raises:
      IOError: if path to Netcdf file not valid (rare under automatic sg usage).

    Examples:

    >>> P.load('O_temp',slices=(Z,0,X,50))
    >>> E.show()
    DPO
    ----------
    O_temp_sliced 
    >>> TEMP = E['O_temp_sliced']
    >>> TEMP.grid
    (latitude)
    """  

# --> this load is a method of Exper

    if not(isinstance(varnames, list)):
      varnames = [varnames ]
    
    for varname in varnames: 
# Big loop.

     
      # Prepare the paths to all the netcdf files into a list     
      if os.path.isfile(self.path):
        # this Exper object was created from a project containing an experiment file
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
            
#	Try to find netcdf var in any of the found files, and then read into Field object.

      fnames = []
      F = []
      for filepath in paths:

        try:       
          file = netcdf_file(filepath,'r', *args, **kwargs)

          # extract the filename fnm from the path:
          fnm = os.path.split(filepath)[1]
        except IOError:
          raise IOError('Cannot open %s'%filepath)
        else:

          if varname in file.variables:
           
            F.append(cdfread(filepath,varname,self.cstack,self.axes,slices =slices, slice_suffix = slice_suffix))
            fnames.append(fnm)

          file.close()
     
      if F == []:
     
        print '%s for %s could not be read.'%(varname,self.name) 

      else:
        num_files = len(paths)

        print 'OK. Fetched Field %s for %s. %s file%s found.'%(varname, self.name,num_files,plural(num_files) )

        new_field = concatenate(F,ax=ax, name_suffix=name_suffix, new_coord_name = new_coord_name, strings = fnames, new_coord= new_coord   )

   
        # insert Field into experiment
        if squeeze_field:
          self.insert(what = (new_field.name, squeeze( new_field ) ) ) 
        else:
          self.insert(what = (new_field.name, new_field ) ) 

    self.update_nbytes() 


  def insert(self,what):
    """
    Insert Field in Field list (vars attribute), or into params att if argument is a number.

    The "what" argument is a 2 tuple (pair) of name and value: (name, value). Value can be a Field or a single value. Name must be a string, but can be None in the case of a Field, where the Field name will then be used. For example what = ('temp',TEMP), where TEMP is a Field. If value is a single value (e.g. int or float), a name must be provided.

    Argument what can also be a list of (name,value) pairs, in which case the entire collection of pairs will be inserted.   

    At the moment, 'insert' forces key van object name to be consistent. Might change in future. Also copies field object to do this.

    Args:
      what: (length 2 tuple or list thereof) name and value: (name, value).

    Returns:
      None
    """

    if what is None:
      return

    if isinstance(what, list):
      for pair in what:
        self.insert(pair)

    else:

      # value here refers to the name, value pairs, not Field value
      name = what[0]
      value = what[1]

      if isinstance(value,Field):  

        if name is None: 
          name = value.name
        else:
          value = value.copy(name = name)
        # insert Field 
        self.vars[name] = value

      else:
      # it is assumed Field is a parameter
        if name is None:
          raise Exception('Provide name if trying to insert what as parameter.')

        self.params[name] = value
        


    


  def available(self):
    """
    Obtain list of all available Netcdf variable names (strings) for this Exper.
   
    Args:
      No args.

    Returns:
      List of strings of variable names in experiment Netcdf file(s).

    Raises:
      IOError: when Netcdf file cannot be opened.

    Called by ls method.

    See also:
      ls method.
    """
  
    if os.path.isfile(self.path):
        # this Exper object corresponds to a file (not directory).
      paths = [self.path]

    else:            
        # this Exper object corresponds to a directory (not file).
        paths = []
        for root, dirs, files in os.walk(top = self.path):
          for fname in files:
            if fname.split('.')[-1] in ['nc','cdf']:
              paths.append(os.path.join(root,fname))

           
# Try to find netcdf var in any of the found files, and then read into Field object.

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
    Show list of all available Netcdf variable names (strings) for this Exper to screen.
   
    Args:
      width: (int) column width formatting for screen output.

    Returns:
      None

    Raises:
      IOError: when Netcdf file cannot be opened.

    Calling available method.

    See also:
      available method.


    Variable list method of Exper.
    Examine which fields (Netcdf variables) are available of experiment object.
    """  
 
    RP = Report()
    av_list = self.available()
    av_list.sort()
    RP.echo(av_list,delim = '\t ',maxlen=300, width = width)

    print RP.value


  def update_nbytes(self):
     """Recalculate and update memory usage of this Exper.
     """

     coord_bytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.cstack]  )
     if len(self.vars) > 0:
       field_bytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.vars.values()]  )    
     else:
       field_bytes = 0
     self.nbytes  = coord_bytes + field_bytes

     
    

# ----- end Exper  ------------






def isexpdir(path, file_extensions = cdf_file_extensions):
  """
  Tests whether the subdirectories in the path contain data files. 

  Returns the list of those subdirectories (relative path to path argument) that 
  contain these known files. To be used by adexp functionality and such.
  file_extensions is the list of known filenames in the form of glob expressions, e.g.   
  ['*.nc','*.cdf'] (the default).   

   
  Args:
    path: (str) path to directory containing experiment directories or files (usually a project dir). Generally computed by sg
    file_extensions: (str) file extensions to look for in path (default .cdf and .nc)

  Returns:
    List of directory and file names (not full paths) believed to correspond to experiments.

  Raises:
    RuntimeError: when path not valid.

  path can be relative to pwd or full path.

  Examples:
  >>> sg.isexpdir('/home/me/PROJECTS/test_project/')
  ['DPO', 'DPC', 'Lev.cdf']
  """

    # examine all subdirectories of path. Create copy for manipulation.

 
  if os.path.isdir(path):
    L = os.listdir(path)
    Lc = copy.deepcopy(L)
  elif os.path.isfile(path):
    return [path,]
  else:
    raise RuntimeError('No such file or directory %s.'%path)
    
  # go through list L of subdirectories of path (e.g. /home/me/PROJECTS/test_project/) to see which ones are experiment directories (i.e. contain .nc and .cdf files).
  for it in L:
    try:
      # look for files ending in .nc or .cdf
      
      exp_path = os.path.join(path , it)
      # Try globbing <path>.nc and <path>.cdf. E.g. /home/me/PROJECTS/test_project/DPC/*.nc for a globfpath
      globfpaths = [os.path.join(exp_path , e) for e in file_extensions]
 
      # Test whether any files of the extensions (e.g .nc) in file_extensions occur in this directory:
      if not(reduce(lambda x,y: x + y, [glob.glob(e) for e in globfpaths])):	
         # if no files with these extensions are found, delete the directory from the list:
        Lc.remove(it)
    except:
     # bad directory anyway
      Lc.remove(it)

  globfpaths = [os.path.join(path , e) for e in file_extensions]
  raw_files = reduce(lambda x,y: x + y, [glob.glob(e) for e in globfpaths])
  files = [e for e in raw_files if os.path.isfile(e)]

  # Returning list of directory and file names (not full paths) believed to correspond to experiments.
  return Lc + [ os.path.split(e)[-1] for e in   files ]

