#encoding:utf-8

""" Projects represent experiment data organized in project directories.

Netcdf data must be organised in a projects directory tree. Create a directory named PROJECTS inside your home directory (name default), to contain all the data for all the projects. Then create a subdirectory for each of your projects, say one is called test_project. When using sg.info() and sg.Project with the nonick = True argument (see documentation), this is all that is needed to allow the sg module to recognise these subdirectories as representing projects. 
  This will allow the sg Python code to regard all Netcdf files placed directly inside test_project as belonging to individual experiments. Alternatively, instead of a single file, a subdirectory may represent an experiment. So directories inside the test_project directory are associated with output from a single experiment each, where all Netcdf files contained within are interpreted as belonging to that same experiment.   

  To get started in ipython (bring up with sg.overview() ):

  >>> import spacegrids as sg		
  >>> D = sg.info(nonick = True)  
  >>> P = sgPproject(D['my_project'] , nonick = True)  
  >>> P.load(['temperature','u'])  
  >>> # obtain the axes under their names T, X, Y, Z in namespace:
  >>> for c in P['some_experiment'].axes:
  >>>   exec c.name + ' = c'	
 
"""

import types
import os
import re

from _config import *

from fieldcls import *
from expercls import *


# ---- temporary classes for old code warnings -----

class project(object):
  """Legacy class to generate error when initialized by depricated code.
  """

  def __init__(self,path=home_path, expnames = ['*'], varnames = [],msk_grid = 'UVic2.8_t', name = None, nonick = False, descr = None, verbose = False):
    """
    In case user tries to invoke class by old name
    """

    raise Exception('sg.project now depreciated. Use sg.Project. All class names now use capital convention.')


# ---------------- Class definition for projects -----------

class Project(object):
  """Represents experiment data organized in project directories.

  Data directory tree:

  Before using projects, Netcdf data must be organised in a projects directory tree. Create a directory named PROJECTS inside your home directory (name default), to contain all the data for all the projects. Then create a subdirectory for each of your projects, say one is called test_project. When using sg.info() and sg.Project with the nonick = True argument (see documentation), this is all that is needed to allow the sg module to recognise these subdirectories as representing projects. 
  This will allow the sg Python code to regard all Netcdf files placed directly inside test_project as belonging to individual experiments. Alternatively, instead of a single file, a subdirectory may represent an experiment. So directories inside the test_project directory are associated with output from a single experiment each, where all Netcdf files contained within are interpreted as belonging to that same experiment.   

  To get started in ipython (bring up with sg.overview() ):

  >>> import spacegrids as sg		
  >>> D = sg.info(nonick = True)  
  >>> P = sgPproject(D['my_project'] , nonick = True)  
  >>> P.load(['temperature','u'])  
  >>> # obtain the axes under their names T, X, Y, Z in namespace:
  >>> for c in P['some_experiment'].axes:
  >>>   exec c.name + ' = c'	
  """

  def __repr__(self):
    return self.name
 
  def __unicode__(self):
    return self.name
    

  def __init__(self,path=home_path, expnames = ['*'], varnames = [],msk_grid = 'UVic2.8_t', name = None, nonick = False, descr = None, verbose = False):

    """
    Initialize Project object to query Netcdf experiment files belonging to project.

    Args:
      path: (str) path to the Project directory (e.g. /home/me/PROJECTS/test_project/). 
      expnames: (str or list of str) filters on experiment names to load
      varnames: (list of str) usually empty, variable list (better to do after init)
      msk_grid: (str) grid to use for masks
      name: (str) project name. Filled using projname content (or dir name) if None.
      nonick: (Boolean) Switch whether to look for a projname text file inside the Project directory to obtain Project nickname. True/ False
      descr: (str) longer description of project.

    Returns:
      Project object to be used to query and load Netcdf files.

    Examples:

    >>> D = sg.info()	# obtain dictionary of nicknames vs full paths
    >>> P = sg.Project(D['my_project'])
    
    
    Or with hard coded path:

    >>> P = sg.Project('/home/me/PROJECTS/test_project/')

    Or when ignoring the projname file inside the project directory: 

    >>> P = sg.Project('/home/me/PROJECTS/test_project/', nonick=True)
    """

    global ax_disp_names

    self.path = path 
# The variables associated with this experiment. It is a dictionary of name vs struct   
    self.expers = {}

    if os.path.isfile(path):
      self.name = 'joep'
    else:

      if name is None:

        # if nonick is selected, a Project nickname file (usually projname) is not required.
        if nonick:
          try:
            with open(os.path.join(self.path, projnickfile  )  ) as f:
              name = read_projnick(f)      
          except:
            # following only works on unix/ linux systems:
            name = end_of_filepath(path)

            if path[-1] != '/':
              path = path + '/'

        else:

          with open(os.path.join(self.path, projnickfile  )  ) as f:
            name = read_projnick(f)      

        self.name = name
 
    # sniff out the netcdf situation:
    
       # find a list of the directories in the Project path that qualify as experiment data dirs, i.e. contain .nc or .cdf files: 
    DL = isexpdir(self.path) 
    # Filter experiment names list against expnames filter argument:

    expnames = sublist(DL,expnames)     

    if verbose:
      print 'Adding: ',
      print expnames

    if descr is None:
      self.descr = name
    else:
      self.descr = descr  
 
    if expnames:
      self.adexp(expnames)
      # after the experiments have been loaded, their Coord elements are all different objects. So check if some of them should be the same object:
      for le in self.expers:
        for re in self.expers:
        # if axes between experiments are different objects but have identical attribute values, make the objects identical. If they only have the same axis attribute value, make them equivalent. 

          find_equal_axes( self.expers[le].cstack , self.expers[re].cstack )

      for e in self.expers:
        # make a handy dictionary for easy reference
        self.expers[e].coords = {c.name:c for c in self.expers[e].cstack}

      # at this stage of the loading process, the Coord elements in the experiment cstack attributes (lists) have an axis attribute, but they are only the strings (e.g. 'X','Y','Z') obtained from the netcdf file. Next, we will replace these string instances with axis objects (defined as a class in sg) with their name attributes being the original strings (e.g. 'X'). Here, as with coords above, we need to be carefully to make those objects the same objects (in memory) if they are equal according to the weaksame method criteria (equal name, values etc, but possibly different objects in memory). So: if a.weaksame(b) True, then b is replaced with a in the list). 
      # So, replace axis string attribute with object attribute and obtain corresponding accumulated list of axis objects. We do this using the make_axis function:
       
      L = make_axes([ c for e in self.expers for c in self.expers[e].cstack ])

      # NOTE THAT ALL EXP OBJECTS HAVE AXES POINTING TO THE SAME LIST IN MEM!!
      for e in self.expers:
        self.expers[e].axes = L      


      for e in self.expers:
        self.expers[e].cstack =  find_set_dual(self.expers[e].cstack)

      # axis objects have now been created based on information in the netcdf files.
      # the exper objects now have a list of axis objects as an attribute, these should normally be identical for all experiments.
      # each Coord object in each cstack attached to each exper object now also has an axis object associated with it.

      # now that the axis objects have been created, we can attempt to assign display names to them, to be used in plotting etc.

     
      for ax_elem in L:
        if ax_elem.name in ax_disp_names:
          ax_elem.long_name = ax_disp_names[ax_elem.name]
        

    if msk_grid:
      # this is a temporary fix
#      y_name = 'yt'
#      x_name = 'xt'
      try:
        y_name, x_name = tgrid_names[msk_grid]
        print 'Trying to read mask using ' + y_name +', ' + x_name + ' for mask.'
        self.masks = read_masks(self.path +mask_dir, grids = [self[ename].coords[y_name]*self[ename].coords[x_name] for ename in self.expers] ) 
      except:    
        print '...no masks read in Project init.'

    if varnames:
      self.load(varnames)

    if len(self.expers.values() ) > 0:
      self.nbytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.expers.values()]  )
    else:
      self.nbytes = 0.


  def ls(self):
    """
    Display experiments belonging to Project in columns.
    """
    RP = Report()

    RP.echo(self.expers.keys(),delim='\t')   
 
    return RP

  def get(self, expnames, varnames):
    """
    Fetch fields from a Project.

    expnames will be matched against the available project experiment names using fnmatch, allowing patterns such as wild cards. The same with varnames.

    Args:
      expnames: (str or list of str) containing experiment name patterns to load from. 
      varnames: (list of str) of variable name patterns to match against 

    Returns:
      List of Field objects from Project with matching attributes.
    """
    
    fields = []
  
    # Search in masks first

    VN = sublist(self.masks.keys(), varnames)
	
    if VN:
      for vn in VN:
	fields.append( self.masks[vn]  )

    # and then in fields	    
    if not(fields): 
      EN = sublist(self.expers.keys() , expnames)    
      if EN:
        for en in EN:
  	  VN = sublist(self.expers[en].vars.keys(), varnames)
	
	  if VN:
	    for vn in VN:
	      fields.append( self.expers[en].get(vn)  )

	    
    if not(fields):
      print "No fields found. Try P.load(%s)"%varnames
    else:
      if len(fields) == 1:
        # if only single Field found, return that Field and not a list of fields.
        fields =fields[0]

    return fields    

  # --> method belongs to Project	    
  def ad_mask(self,sub_dir = mask_dir,grid = False, msk_val = 2):
    """Read text file masks from file and add to project.
    """
#    print read_masks(self.path + sub_dir, grid = grid, msk_val = msk_val)
    self.masks = read_masks(self.path + sub_dir, grid = grid, msk_val = msk_val)

    
  def __getitem__(self,val):
    """
    Returns an exper object by the name of argument 'val'.
    """
# --> method belongs to Project
 
    if val in self.expers.keys():
      return self.expers[val]
    elif (val == slice(None, None, None)) | (val == slice(0, 9223372036854775807, None)):
      return [self.expers[k] for k in self.expers]

    else:
      raise Exception("io: exper %s not loaded." % val)


# --> belongs to Project      
  def show(self): 
    """
    Display summary of experiments and fields loaded in Project.
    """

    RP = Report('---- %s ----\n'%self.name)
#    RP.echo('Experiments: ')
    explist = self.expers.keys()
    explist.sort()
   
    RP.echoln(explist,delim='\t')

    if not(self.expers.keys()):

      RP.echo('No experiments in Project.')
      print RP.value
      return

    RP.echoln()
    RP.echo('Project using %1.2f Mb.'%(self.nbytes/1024./1024.) )

    print RP.value
    
  def incdf(self, varname, expname = ''):
    """
    Check whether variable name is in the Netcdf files.

    Args:
      varname: (str) variable name to check
    
    Returns:
      True/ False depending on whether the variable name is in the netcdf database.
    """

 # --> belongs to Project

    exps = self.expers 
    if expname == '':
      if len(exps.keys()) > 0:
        exp = exps[exps.keys()[0]]
      else:
        print "Load exps first."
        return False
    else:
      exp = exps[expname]
    
    if varname in exp.var_names:
      return True
    else:
      return False

 # --> belongs to Project
 
  def delvar(self,varname):
    if self.expers: 
      for expname in self.expers.keys():
        self.expers[expname].delvar(varname) 

        self._update_nbytes()

    else:
      print 'No experiments.'

  def delexp(self, expnames, msg = " "):
  
    if not(isinstance(expnames,types.ListType)):   
      expnames = [ expnames ]
    
    for expname in expnames:

      if expname in self.expers:
        del self.expers[expname]
        print 'Deleting existing '+ expname   	      
      else:
        if msg:
          print 'No experiment by that name. ' + msg

    self._update_nbytes()


# --> belongs to Project

  def adexp(self, expnames=['*'], descr = 'An experiment.'):
    """Adds an experiment to this Project. 

    Also looks for a sub-directory "masks" for masks. If it can't find that, it will use the mask inherited from the Project.

    Args:
      expnames: (str or list of str) of experiment name patterns to match
      descr: (str) optional description of experiment    
    """

    # Convert argument expnames to list if it is not yet a list, e.g. 'flx_BL' to ['flx_BL']
    if not(isinstance(expnames,types.ListType)):
      expnames = [ expnames ]
    
    # find a list of the directories in the Project path that qualify as experiment data dirs. 
    DL = isexpdir(self.path)
    
    # match the argument expnames against DL via wildcard expansion etc.
    expnames = sublist(DL,expnames)     
    
    for expname in expnames:
       
      # removes this exp if it is already present.
      self.delexp(expname,'')  	

# Create new experiment object projectd on name
 
      print 'Creating exp object ' + expname +' from ' + self.path

      cstack = cdfsniff(os.path.join(self.path, expname) )

      self.expers[expname] = Exper(path = self.path,name = expname, cstack = cstack, descr = descr)

    self._update_nbytes()
    return 

  def load(self, varnames, descr = 0,chk_loaded = False, ax=None, name_suffix='_cat', new_coord_name = 'gamma', new_coord= None, slices = None ):
# -->This load belongs to Project. It calls the load method of the Exper. 
# slices argument value is just an example
    """
    Load Field into Project from Netcdf files.
    
    Loads a variable (or list of variables) for all experiments belonging to the Project instance. E.g. P.load('salinity').

    Attempts to concatenate if the same variable name is found in multiple files.
 
    Args:
      varnames: (str or list of str) of variable names to load
      descr: (str) optional description     
      chk_loaded: (Boolean)  if True (default False), makes load method check whether the Field is already loaded. If so, it is not loaded.
      ax (Ax): passed on to the concatenate function from Exper.load
      name_suffix: (str) suffix to apply to name
      new_coord_name (str): passed on to the concatenate function
      new_coord (Coord): passed on to the concatenate function
	
    Examples:
   
    >>> P.load('salinity')
    >>> E = P['CNTRL_2200p']	# pick this particular experiment
    >>> S = E['salinity']		# create Field object and assign it to S
    """

# errno is returned by this method and gives information about errors encountered in this method.    
    errno = 0

# list functionality only exists at the level of the Project .
    if not(isinstance(varnames, list)):
      varnames = [ varnames]
   
    for varname in varnames:
 
# Test if experiments have been registered with the Project:    
      if not(self.expers):
        print "No experiments." 
# errno signifies that there are no experiments to load for.
        errno = 1
      else:

# See if chk_loaded is set. If so, check whether the Field is already loaded for the fist experiment.
# If so, do not attempt to load again (but keep what's been loaded).
        if chk_loaded:
# chk_loaded flags whether it should be checked whether the Field is already loaded. Usually False.
          if self.loaded(varname):
            # This error number signifies a warning, and means that the variable is already loaded.
            errno = -1
	    continue
# Now the actual netcdf calling part.

	    # in case netcdf information has been pre-added:
   	for k in self.expers.keys():
            
	  self.expers[k].load(varname,ax=ax, name_suffix=name_suffix, new_coord_name = new_coord_name, new_coord= new_coord, slices = slices)  

    self._update_nbytes()
    return    

  def write(self, path = None, name = None, force = False , history = 'Created from Spacegrids ', insert_dual = True  ): 
    """
    Write entire loaded Project to disk. 
    """     


    if path is None:
      path = os.path.join(home_path,"PROJECTS")

      if name is None:
        name = self.name
     
      path = os.path.join(path,name)

    if not os.path.exists(path):
      os.mkdir(path)
      print 'Creating dir %s'%path  

    if not(force):
      if os.path.samefile(path , self.path ):
        print 'Refusing to write Project %s to own path at risk of overwriting data. Use different path or set force = True'%name
        return

    # write Project name text file:
    with open(os.path.join( path, projnickfile ),'w') as f:    
      f.write(name)

    for e in self.expers.values():
      e.write(path = path, insert_dual = insert_dual)
      


  def loaded(self,varname):
    """Checks whether variable by name of (str) varname is loaded.
    """
    flag = False
    
    if (self.expers):
      ex=self.expers
      if varname in ex[ex.keys()[0]].vars.keys():
	flag = True
	
    return flag	
 
  def _update_nbytes(self):
    """
    Update estimate of number of Mb used by Project.
    """
    if len(self.expers.values() ) > 0:
      self.nbytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.expers.values()]  )

    else:
      self.nbytes = 0



  def cat(self, fld_name, new_ax_name = 'exper', new_coord_name = None, new_coord = None):

    """
    Concatenate the Field fld_name for all experiments where that Field is loaded.
  
    If new_coord coord object argument is given, other arguments are overridden and that argument is used as the new coord along which the different experiments are lined up.

    if no new_coord argument is given, a new Coord will be constructed.
 
    Args:
      fld_name: (str) name of Field objects to concatenate
      new_ax_name: (str) name of new Ax object to create to concat along if desired
 
      new_coord_name: (str) name to be used for creation of new Coord in case ax not in Field grid
      new_coord: (Coord) overrides default behaviour if set by concatenating along that new Coord.
    """

    e_loaded = [e for e in self.expers.keys() if fld_name in self.expers[e].vars]
    if e_loaded ==[]:
      warnings.warn('Field %s not loaded. Returning None.'%fld_name) 
      return     

    # If new_coord has been provided, use that:
    if new_coord is not None:
      fields = [self.expers[e].vars[fld_name] for e in e_loaded]
    
      if len(fields) != len(new_coord):
        raise Exception('Provide equal amount of loaded Field instances to length new_coord')

      # EXIT POINT
      return concatenate(fields,  new_coord = new_coord)
      

    # if new_coord has not been provided, construct one.

    if new_coord_name is None:
      new_coord_name = new_ax_name + '_crd'

    e_loaded.sort()

    fields = [self.expers[e].vars[fld_name] for e in e_loaded]

    new_axis = Ax(new_ax_name)   
    new_coord_value = np.arange(len(fields))
    new_coord = Coord(name = new_coord_name , value = new_coord_value, axis = new_axis, direction = new_axis.name, strings = e_loaded )

    # EXIT POINT  
    return concatenate(fields,  new_coord = new_coord)


  def param2gr(self, param_name, func, name_filter = None, sort = True, new_Ax_name = None, add2cstack = True):
    """
    Use a param from the experiments to construct a Coord, and so construct a new concatenated Field defined on a new grid with this new Coord.

    func is a function that must yield a Field and take an Exper object as argument.
    it should return None for elements that are not desired. This method is somewhat different from its near-opposite method insert, where the function argument is allowed to yield fields or single values (float/ int).

    The idea is: a parameter belongs to an Exper, so any function taking an Exper argument and returning a Field of certain dimensions allows us to concatenate these Field objects along a new Coord constructed from the parameters. 

    Args:
      param_name: (str) parameter name
      func: (function) that must yield a Field and take Exper object as argument
      name_filter: (str) pattern to filter only certain Exper names.
      sort: (Boolean) if True, sorts the Coord
      new_Ax_name: (str) name for new Ax object corresponding to new param-based Coord, can be None (it will then be constructed)
      add2cstack: (Boolean) determines whether new Coord will be added to the Exper coordinate stacks.

    Returns:
      The concatenated Field.
    """

    if name_filter is None:
      expers = self.expers.values()
    else:
      exp_keys = simple_glob([k for k in self.expers.keys() ], name_filter )
      expers = [self.expers[k] for k in exp_keys]


    if new_ax_name is None:
      new_ax_name = param_name

    new_coord_name = new_ax_name + '_crd'

    pairs = [(E.params[param_name], func(E) ) for E in expers if param_name in E.params]

    pairs = [e for e in pairs if e[1] is not None   ] 

    if sort:
      pairs.sort()


    new_ax = Ax(new_ax_name)

    if add2cstack:
      for E in self.expers.values():
        i = new_ax.sameindex(E.axes)
        if i is not None:
          new_ax = E.axes[i]
          break


      # ax has to be added to ONLY 1 EXP. This is because axes attributes of all experiments point to the same list of axis objects!!!
      if i is None:
        if not new_ax in self.expers.values()[0].axes:
          self.expers.values()[0].axes.append(new_ax)  

    def func(x): return np.nan if x is None else x 

    new_coord = Coord(name = new_coord_name, value = np.array( [func(e[0]) for e in pairs] ), axis = new_ax, direction = new_ax.name)

    if add2cstack:
      for E in self.expers.values():
        i = new_coord.sameindex(E.cstack)
        if i is not None:
          new_coord = E.cstack[i]
          break
            
      if i is None:
        for E in self.expers.values():
          if not new_coord in E.cstack:
            E.cstack.append(new_coord)  


    return concatenate(fields = [e[1] for e in pairs], new_coord = new_coord)

  def insert(self, func):
    """
    Apply function func to each Exper object in self and insert the result into Project. 

    The result of func must be a tuple (name, value) or a list of such tuples (as taken by the Exper insert method). Conform the Exper insert method, name must be a string, and value a Field or ordinary float/ int kind of datatype.

    Args:
      func: (function) taking Exper argument and returning list of (name,value) pairs.  
    """
    for E in self.expers.values():
      if E is not None:  
        E.insert( func(E)  )


  def pattern2gr(self,fld_name,pattern, parname = None, name_filter = None,nomatch_fill = None):
    """Use rexexp patterns on experiment file names to extract parameters and apply param2gr.
 
    The first grouping in the regexp will be used to yield the parameter.

    Args: (as in param2gr, except one)
      pattern: (str) regular expression to apply to file names.
      
    Returns:
      The concatenated Field.

    Examples:

    >>> regexp = '[\w_]+_(\d+)ppm[\w]+' # the grouping between brackets yields the value
    >>> bigfield_DP=P.pattern2gr(fld_name = 'temp',pattern = regexp,parname='co2', name_filter = 'DP*', nomatch_fill = 280.)
    >>> # this is for filenames with co2 values in their file name.
    """
# expect [(parname, value),]  or None

    self.insert(parse_fname_func(pattern,parname=parname, nomatch_fill = nomatch_fill))

    return self.param2gr(parname , lambda x:x[fld_name], name_filter = name_filter) 

# ---------------- End  definition for projects -----------


# -------- Exper helper functions to be used with insert method ------------

def read_control_func(filename):
  """
  function that yields a useful function that can be used in the Project insert method: it reads a control file in the MOM 2 control file format (as used by UVic) and outputs it in the list of (name, value) tuples format required by the Exper insert method. Think of it as inserting a parameter to each Exper, although it can be a Field too.

  The argument filename is the file name of the control file to be read, usually control.in.

  For example: P.insert(sg.read_control_func('control.in')) goes through all experiment directories and looks for, and parses, the file control.in, and then inserts the result as Exper parameters.

  The user could construct similar functions for their own use: the function must take an experiment as argument, and return None or a list of (name, value) pairs (2 tuples).
  """

  def get_controls(E):

    path = os.path.join(E.path, filename)
    try:
      La, L = parse_control_file(path)
    except:
      La = None
    return La
  

  return get_controls

def parse_fname_func(pattern, parname = None,value_type = float, nomatch_fill = None):
  prog = re.compile(pattern)  
  """Takes a regular expression, and the intended parameter name, to use it to yield another function that extracts the parameter from experiment (file) names.

  Args:
    pattern: (str) regular expression to match names against
    parname: (str) parameter name
    value_type: (float or int etc)
    nomatch_fill: value to fill points that don't match with.

  Returns:
    List of (name, value) pairs to be used by insert method.

  Examples:

  >>> P.insert(sg.parse_fname_func('[\w_]+_(\d+)ppm[\w]+',parname='co2'))
  """

  if parname is None:
    def get_fname(E):
      M =prog.match(E.name);


      if M is not None:
        parname = M.groups()[0]
        value = value_type(M.groups()[1])
      else:
        return 

      return [(parname, value),]
  

    return get_fname

  else:
    def get_fname(E):
      M =prog.match(E.name);


      if M is not None:
        value = value_type(M.groups()[0])
      else:
        return [(parname, nomatch_fill),]

      return [(parname, value),]
  

    return get_fname



# -----------------------------
# functions to be used for creating new coords

def avg_temp(P, varname = 'O_temp'):
  """Example function to create avg ocean temperature from expers 
  """
  for c in P.expers.values()[0].axes:
    exec c.name + ' = c'

  for k in P.expers.keys():
    if varname in P.expers[k].vars:
      mTEMP = P.expers[k].vars[varname]/(X*Y*Z)

      P.expers[k].insert(mTEMP, name = 'avg_temp')

#  return P    

# ------- general functions -----------

def overview():
  """Prints quick overview of sg to screen.
  """
  RP = Report()
  RP.echoln("Quick start instructions.")
  RP.line()
  RP.echoln()
  RP.echoln("# To get started: ")
  RP.echoln("D = sg.info(nonick = True)	# look for projects in ~/PROJECTS/")
  RP.echoln("P = sg.Project(D['my_project'] , nonick = True)	# init Project")
  RP.echoln("P.load(['temperature','u'])     # load some variables into mem")
  RP.echoln("for c in P['some_experiment'].axes:   ")
  RP.echoln("  exec c.name + ' = c'	# bring axes names into namespace")
  RP.echoln("for c in P['some_experiment'].cstack:   ")
  RP.echoln("  exec c.name + ' = c'	# bring coord names into namespace")
  RP.echoln()
  RP.echoln("# *** Some examples ***   ")
  RP.echoln("TEMP = P['some_experiment']['temperature']   ")
  RP.echoln("U = P['some_experiment']['u']   ")
  RP.echoln("TEMP_sliced = TEMP[Y,:50]     ")
  RP.echoln("m_TEMP = TEMP_sliced/(X*Y)   ")
  RP.echoln("TEMP_regridded = TEMP(U.grid)	   ")
  RP.echoln("   ")

  print RP.value

def ls(rootdir = os.environ['HOME'], projdirname = 'PROJECTS',fname = projnickfile, nonick = False):
  """
  Simple function to display all Project directories vs their projname nicknames.

  See info method for details.

  Calls info_dict.  
  """

  D = info_dict(rootdir = rootdir, projdirname = projdirname,fname = fname, nonick = nonick)
 
  RP = Report()
  RP.echoln("Projects rooted in %s"%projdirname)
  RP.line()
  project_names = D.keys()
  project_names.sort()
  RP.echo(project_names,width = 17,maxlen=100)

  print RP
  


def info(rootdir = os.environ['HOME'], projdirname = 'PROJECTS',fname = projnickfile, nonick = False, verbose = True):
  """
  Simple function to take inventory of all Project directories based on contents of projname text files inside project directories.  

  This way, no specific paths need to be used, and projects can be referred to by their nicknames defined in a file called projname in each directory containing the experiment directories.
 
  Args: 
    rootdir: (str) dir in which to look for main PROJECTS dir
    projdirname: (str) dir in which to look for specific Project directories. Default is 'PROJECTS'. So by default, projects are expected in '~/PROJECTS/').
    fname: (str) the filename to look for and read the content of to determine nickname of specific Project. Disabled with nonick = True
  
  Returns:
    D, dictionary of Project nicknames vs full paths		
    
  If nonick is False:  finds all dir paths with file called fname in them, reads that file for each dir to find the nickname of that Project and builds a dictionary of nicknames vs paths. Otherwise, just takes inventory of all sub directories of projects (~/PROJECTS)
  
  Example of use when projects are in default location ~/PROJECTS/ (spacegrids loaded as sg):
  
  >>> D=sg.info()
  
  Say this gives keys 'test' and 'TKE'. To open the Project with nickname test:

  >>> P = sg.Project(D['test'])
  
  and start working with P. No specific path had to be identified, only a sufficiently specific top.
  """

  D = info_dict(rootdir = rootdir, projdirname = projdirname,fname = fname, nonick = nonick)
  
  if verbose:

    ls(rootdir = rootdir, projdirname = projdirname,fname = fname, nonick = nonick)
  
  return D

def info_dict(rootdir = os.environ['HOME'], projdirname = 'PROJECTS',fname = projnickfile, nonick = False):
  """
  See info method.
  """


  top = os.path.join(rootdir, projdirname)

  
#  first find all the file paths using locate function defined above
  if nonick:
    paths = locate(top = top,fname = None)

  else:
    paths = locate(top = top,fname = fname)

# initialise dictionary  
  D = {}


  if nonick:

    for path in paths:
      # fname contains the Project nickname file (usually projname)       
      fp = os.path.join(path,fname)
      try:
        with open(fp,'r') as f:
# read the first line and make sure the return char \n is deleted.    
          projnick = read_projnick(f)

      except:
          projnick = end_of_filepath(path)
# ad to dictionary    
      D[projnick] = os.path.join(path,'')

  else:

    for path in paths:
      # fname contains the Project nickname file (usually projname)
      fp = os.path.join(path,fname)
      with open(fp,'r') as f:
# read the first line and make sure the return char \n is deleted.    
        projnick = read_projnick(f)
  
# ad to dictionary    
        D[projnick] = os.path.join(path,'')

  
  return D


    
def read_projnick(f):
  """Reads project nickname file (usually "projname") in project directory.
  """  
  return f.readline().rstrip()







