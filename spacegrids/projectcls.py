#encoding:utf-8

""" project class
"""

import types
import os

from config import *

from fieldcls import *
from expercls import *



# ---------------- Class definition for projects -----------

class project:

  def __repr__(self):
    return self.show()
    

  def __init__(self,path=home_path, expnames = ['*'], varnames = [],msk_grid = 'UVic2.8_t', name = None, descr = None, verbose = False):

    global ax_disp_names

    self.path = path 
# The variables associated with this experiment. It is a dictionary of name vs struct   
    self.expers = {}

    if os.path.isfile(path):
      self.name = 'joep'
    else:

      if name is None:

        with open(os.path.join(self.path, projnickfile  )  ) as f:
          name = read_projnick(f)      
          self.name = name

          
    

    # sniff out the netcdf situation:
    
       # find a list of the directories in the project path that qualify as experiment data dirs, i.e. contain .nc or .cdf files: 
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
      # after the experiments have been loaded, their coord elements are all different objects. So check if some of them should be the same object:
      for le in self.expers:
        for re in self.expers:
        # if axes between experiments are different objects but have identical attribute values, make the objects identical. If they only have the same axis attribute value, make them equivalent. 

          find_equal_axes( self.expers[le].cstack , self.expers[re].cstack )

      for e in self.expers:
        # make a handy dictionary for easy reference
        self.expers[e].coords = {c.name:c for c in self.expers[e].cstack}

      # at this stage of the loading process, the coord elements in the experiment cstack attributes (lists) have an axis attribute, but they are only the strings (e.g. 'X','Y','Z') obtained from the netcdf file. Next, we will replace these string instances with axis objects (defined as a class in sg) with their name attributes being the original strings (e.g. 'X'). Here, as with coords above, we need to be carefully to make those objects the same objects (in memory) if they are equal according to the & method criteria (equal name, values etc, but possibly different objects in memory). So: if a&b True, then b is replaced with a in the list). 
      # So, replace axis string attribute with object attribute and obtain corresponding accumulated list of axis objects. We do this using the make_axis function:
       
      L = make_axes([ c for e in self.expers for c in self.expers[e].cstack ])
      for e in self.expers:
        self.expers[e].axes = L      


      for e in self.expers:
        self.expers[e].cstack =  find_set_dual(self.expers[e].cstack)

      # axis objects have now been created based on information in the netcdf files.
      # the exper objects now have a list of axis objects as an attribute, these should normally be identical for all experiments.
      # each coord object in each cstack attached to each exper object now also has an axis object associated with it.

      # now that the axis objects have been created, we can attempt to assign display names to them, to be used in plotting etc.

     
      for l in L:
        if l.name in ax_disp_names:
          l.display_name = ax_disp_names[l.name]
        

    if msk_grid:
      # this is a temporary fix
#      y_name = 'yt'
#      x_name = 'xt'
      try:
        y_name, x_name = tgrid_names[msk_grid]
        print 'Trying to read mask using ' + y_name +', ' + x_name + ' for mask.'
        self.masks = read_masks(self.path +mask_dir, grids = [self[ename].coords[y_name]*self[ename].coords[x_name] for ename in self.expers] ) 
      except:    
        print '...no masks read in project init.'

    if varnames:
      self.load(varnames)

    if len(self.expers.values() ) > 0:
      self.nbytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.expers.values()]  )
    else:
      self.nbytes = 0.


  def ls(self):
    """
    Method of project class to display experiments in columns.
    """
    RP = report()

    RP.echo(self.expers.keys(),delim='\t')   
 
    return RP

  def get(self, expnames, varnames):
    """
    Project method to fetch fields from a project.
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
# if only single field found, return that field and not a list of fields.
        fields =fields[0]

    return fields    

  # --> method belongs to project	    
  def ad_mask(self,sub_dir = mask_dir,grid = False, msk_val = 2):

#    print read_masks(self.path + sub_dir, grid = grid, msk_val = msk_val, parent = self)
    self.masks = read_masks(self.path + sub_dir, grid = grid, msk_val = msk_val, parent = self)

    
  def __getitem__(self,val):
    """
    Returns an exper object by the name of key k
    """
# --> method belongs to project
 
    if val in self.expers.keys():
      return self.expers[val]
    elif (val == slice(None, None, None)) | (val == slice(0, 9223372036854775807, None)):
      return [self.expers[k] for k in self.expers]

    else:
      raise Exception("io: exper %s not loaded." % val)
    
  def __call__(self, expnames='first', varnames='show'):
 
  # --> method belongs to project
  
    if varnames == 'show':
# summary functionality invoked with choosing argument varnames = 'show'
      print "Project path " + self.path + '\n'

      print "Experiments: "
      print 20*"-"
      print_box(self.expers.keys())
      print '\n'
      self.show()
      return
    elif varnames == 'cdf':
      self.cdf()
      return
   
    if expnames == 'first':
      expnames = [self.expers.keys()[0]]

           
    if not isinstance(varnames, types.ListType):
      varnames = [varnames]
    
    if not isinstance(expnames,types.ListType):
      expnames = [expnames]

    if varnames == [None]:
      return_val = []
      for ee in expnames:
        if ee in self.expers.keys():
          return_val.append(self.expers[ee]) 
      if len(return_val) == 1:
	return_val = return_val[0]
      return return_val      
    else:
  
      fields = self.get(expnames = expnames, varnames = varnames)
	
      return fields  

# --> belongs to project      
  def show(self): 
    """
    Display summary of experiments and fields loaded in project.
    """

    RP = report('---- %s ----\n'%self.name)
#    RP.echo('Experiments: ')
    explist = self.expers.keys()
    explist.sort()
   
    RP.echoln(explist,delim='\t')

    if not(self.expers.keys()):

      RP.echo('No experiments in project.')
      return RP.value

#    exp_name1 = self.expers.keys()[0]
#    vars = self.expers[exp_name1].vars.keys()
    
#    dims = []
#    for thisvar in vars:
      
#      shapes_list = []
#      for expname in self.expers.keys():
#        if thisvar not in self.expers[expname].vars or self.expers[expname].vars[thisvar] is None:
#          shape_str = ''
#        else:
#          shape_str =  str(self.expers[expname].vars[thisvar].shape)
        
#        if not(shape_str in shapes_list):
#          shapes_list.append(shape_str)

#      dims.append(reduce(lambda x,y: x+','+y, shapes_list))
    
#    vars = dict(zip(vars, dims  ))
    
#    lvars = [e[0] + ' ' + e[1] for e in zip(vars, dims  )]

#    RP.echoln('Loaded fields (and sizes):')
#    RP.line()
#    if lvars:
#      RP.echo(lvars,width =30, cols=3)
#    else:
#      print 'NO VARS'

    RP.echoln()
    RP.echo('Project using %1.2f Mb.'%(self.nbytes/1024./1024.) )

    return RP.value
    
  def incdf(self, varname, expname = ''):
    """
    Check whether variable name is in the netcdf database.

    Input:
    varname	variable name to check
    
    Output:
    True/ False depending on whether the variable name is in the netcdf database.

    """

 # --> belongs to project

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

 # --> belongs to project
 
  def delvar(self,varname):
    if self.expers: 
      for expname in self.expers.keys():
        self.expers[expname].delvar(varname) 

        self.update_nbytes()

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

    self.update_nbytes()


# --> belongs to class "project"

  def adexp(self, expnames=['*'], descr = 'An experiment.'):
    """
    Adds an experiment class to this project. Also looks for a sub-directory "masks" for masks. If it can't find that, it will use the mask inherited from the project.
    """

    # Convert argument expnames to list if it is not yet a list, e.g. 'flx_BL' to ['flx_BL']
    if not(isinstance(expnames,types.ListType)):
      expnames = [ expnames ]
    
    # find a list of the directories in the project path that qualify as experiment data dirs. 
    DL = isexpdir(self.path)
    
    # match the argument expnames against DL via wildcard expansion etc.
    expnames = sublist(DL,expnames)     
    
    for expname in expnames:
       
      # removes this exp if it is already present.
      self.delexp(expname,'')  	

# Create new experiment object projectd on name
 
      print 'Creating exp object ' + expname +' from ' + self.path

      cstack = cdfsniff(os.path.join(self.path, expname) )

      self.expers[expname] = exper(path = self.path,name = expname, cstack = cstack, descr = descr, parent = self)

    self.update_nbytes()
    return 

  def load(self, varnames,filename=None, descr = 0,chk_loaded = False):
# -->This load belongs to the project class. It calls the load method of the exper class. 
# slices argument value is just an example
    """
    Method of project class. E.g. P.load('salinity') for project P.
    Loads a variable (or list of variables) for all experiments belonging to the project instance.
    E.g. 
    
    P.load('salinity')
    E = P['CNTRL_2200p']	# pick this particular experiment
    S = E['salinity']		# create field object and assign it to S
    
   
    chk_loaded = True (default False) makes load method check whether the field is already loaded. If so, it is not loaded.
    
    """

# errno is returned by this method and gives information about errors encountered in this method.    
    errno = 0

# list functionality only exists at the level of the "project" class.
    if not(isinstance(varnames, list)):
      varnames = [ varnames]
   
    for varname in varnames:
 
# Test if experiments have been registered with the project:    
      if not(self.expers):
        print "No experiments." 
# errno signifies that there are no experiments to load for.
        errno = 1
      else:

# See if chk_loaded is set. If so, check whether the field is already loaded for the fist experiment.
# If so, do not attempt to load again (but keep what's been loaded).
        if chk_loaded:
# chk_loaded flags whether it should be checked whether the field is already loaded. Usually False.
          if self.loaded(varname):
            # This error number signifies a warning, and means that the variable is already loaded.
            errno = -1
	    continue
# Now the actual netcdf calling part.

	    # in case netcdf information has been pre-added:
   	for k in self.expers.keys():
            
	  self.expers[k].load(varname,filename = filename)  

    self.update_nbytes()
    return    

  def write(self, path = None, name = None, force = False , history = 'Created from Spacegrids ', insert_dual = True  ): 
    """
    Write entire loaded project to disk. 
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
        print 'Refusing to write project %s to own path at risk of overwriting data. Use different path or set force = True'%name
        return

    # write project name text file:
    with open(os.path.join( path, projnickfile ),'w') as f:    
      f.write(name)

    for e in self.expers.values():
      e.write(path = path, insert_dual = insert_dual)
      

  def getbody(self, varname, expnames = '*', loadflag = 1, verbose =0):

# --> This getbody belongs to project class.
# gets a list of numpy arrays of body belonging to varname for all experiments.
# loadflag = 1 means: load the variable into the class if not present. 

    expers = subdict(self.expers, expnames)

    bodies = []
    
    if not(self.expers):
      print "No experiments." 
    else:    
      for expname in expers.keys():
        bodies.append(expers[expname].getbody(varname, loadflag))
        if verbose:
	  print 'Retrieved ' + varname + ' for ' + expname
	
    return bodies

  def loaded(self,varname):
  
    flag = False
    
    if (self.expers):
      ex=self.expers
      if varname in ex[ex.keys()[0]].vars.keys():
	flag = True
	
    return flag	
    
  def report(self, varname):
    
    if not(self.expers):
      print "No experiments." 
    else:    

      for expname in self.expers.keys():
        if varname in self.expers[expname].vars:   

          print expname + ' ' +  varname + ' of len ' + str(len(self.expers[expname].getbody(varname)))
    
        else:
	  print varname +' not in experiment ' + expname

  def update_nbytes(self):
    """
    Update estimate of number of Mb used by project.
    """
    if len(self.expers.values() ) > 0:
      self.nbytes = reduce( lambda x,y: x + y, [e.nbytes for e in self.expers.values()]  )

    else:
      self.nbytes = 0



  def cat(self, fld_name, new_ax_name = 'exper', new_coord_name = None, new_coord = None):

    """
    Concatenate the field fld_name for all experiments where that field is loaded.
  
    If new_coord coord object argument is given, other arguments are overridden and that argument is used as the new coord along which the different experiments are lined up.

    if no new_coord argument is given, a new coord will be constructed.
 
    """

    e_loaded = [e for e in self.expers.keys() if fld_name in self.expers[e].vars]
    if e_loaded ==[]:
      warnings.warn('Field %s not loaded. Returning None.'%fld_name) 
      return     

    # If new_coord has been provided, use that:
    if new_coord is not None:
      fields = [self.expers[e].vars[fld_name] for e in e_loaded]
    
      if len(fields) != len(new_coord):
        raise Exception('Provide equal amount of loaded field instances to length new_coord')

      # EXIT POINT
      return concatenate(fields,  new_coord = new_coord)
      

    # if new_coord has not been provided, construct one.

    if new_coord_name is None:
      new_coord_name = new_ax_name + '_crd'

    e_loaded.sort()

    fields = [self.expers[e].vars[fld_name] for e in e_loaded]

    new_axis = ax(new_ax_name)   
    new_coord_value = np.arange(len(fields))
    new_coord = coord(name = new_coord_name , value = new_coord_value, axis = new_axis, direction = new_axis, strings = e_loaded )

    # EXIT POINT  
    return concatenate(fields,  new_coord = new_coord)


  def param2gr(self, param_name, func, name_filter = None, sort = True, new_ax_name = None):
    """
    Use a param from the experiments to construct a coord, and so construct a new concatenated field defined on a new grid with this new coord.

    func is a function that must yield a field and take an exper object as argument.
    it should return None for elements that are not desired. This method is somewhat different from its near-opposite method insert, where the function argument is allowed to yield fields or single values (float/ int).

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

    new_ax = ax(new_ax_name)
    new_coord = coord(name = new_coord_name, value = np.array( [e[0] for e in pairs] ), axis = new_ax, direction = new_ax)

    return concatenate(fields = [e[1] for e in pairs], new_coord = new_coord)

  def insert(self, func, param_name= None):
    """
    Apply function func to each exper object in self and insert the result. The result must be a tuple (name, value) or a list of such tuples (as taken by the exper insert method). Conform the exper insert method, name must be a string, and value a field or ordinary float/ int kind of datatype.

   
    """
    for E in self.expers.values():
      
      E.insert( func(E)  )

  

# ---------------- End class definition for projects -----------


def read_control_func(filename):
  """
  function that yields a useful function that can be used in the project insert method: it reads a control file in the MOM 2 control file format (as used by UVic) and outputs it in the list of (name, value) tuples format required by the exper insert method.

  The argument filename is the file name of the control file to be read, usually control.in.

  For example: P.insert(sg.read_control_func('control.in')) goes through all experiment directories and looks for, and parses, the file control.in, and then inserts the result as exper parameters.

  The user could construct similar functions for their own use: the function must take an experiment as argument, and return None or a list of (name, value) pairs.


  """

  def get_controls(E):

    path = os.path.join(E.path, filename)
    try:
      La, L = parse_control_file(path)
    except:
      La = None
    return La
  

  return get_controls

def info(rootdir = os.environ['HOME'], projdirname = 'PROJECTS',fname = projnickfile, verbose = True):
  """
  Simple function to take inventory of all project directories so that no specific paths need to be used, and projects can be referred to by their nicknames defined in a file called projname in each directory containing the experiment directories.
 
  Inputs: 
  top		(default '~/PROJECTS/') the start dir
  fname		(default projname) the filename to look for and read the content of
  
  Outputs:
  D, dictionary of project nicknames vs paths		
    
  finds all dir paths with file called fname in them, reads that file for each dir to find the nickname of that project and builds a dictionary of nicknames vs paths
  
  Example of use when projects are in default location ~/PROJECTS/ (spacegrids loaded as sg):
  
  D=sg.info()
  
  Say this gives keys 'test' and 'TKE'. To open the project with nickname test:

  P = sg.project(D['test'])
  
  and start working with P. No specific path had to be identified, only a sufficiently specific top.

  
  """


  top = os.path.join(rootdir, projdirname)

  
#  first find all the file paths using locate function defined above
  paths = locate(top = top,fname = fname)

# initialise dictionary  
  D = {}
  for path in paths:
    fp = os.path.join(path,fname)
    with open(fp,'r') as f:
# read the first line and make sure the return char \n is deleted.    
      projnick = read_projnick(f)
  
# ad to dictionary    
      D[projnick] = os.path.join(path,'')
  
  if verbose:

    RP = report()
    RP.echoln("Projects rooted in " + top)
    RP.line()
    project_names = D.keys()
    project_names.sort()
    RP.echo(project_names,width = 17,maxlen=100)

    print RP
  
  return D
    
def read_projnick(f):
  
  return f.readline().rstrip()


# -----------------------------
# functions to be used for creating new coords

def avg_temp(P, varname = 'O_temp'):

  for c in P.expers.values()[0].axes:
    exec c.name + ' = c'

  for k in P.expers.keys():
    if varname in P.expers[k].vars:
      mTEMP = P.expers[k].vars[varname]/(X*Y*Z)

      P.expers[k].insert(mTEMP, name = 'avg_temp')

#  return P    




