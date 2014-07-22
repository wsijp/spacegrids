#encoding:utf-8

""" Coord, Gr, Ax, AxGr and Field class
"""

import numpy as np
import numpy.ma as ma
from scipy.interpolate import griddata

import inspect
import datetime
import glob
import copy
import warnings

from config import *

import os
import copy

#from spacegrids.decorators import field2cumsum

from utilsg import *
from iosg import *

warnings.formatwarning = warning_on_one_line

# ---- decorators -----

def field2cumsum(func):
  """
  Decorator to dispatch argument to other methods depending on argument.

  Args:
    context: (dict) data type vs method name dictionary used in decorator construction.

  Returns:
    Function: the decorator.
  """
  
  def dispatcher(*args,**kwargs):
    """
    Decorator function. Expect args rather than kwargs. 

    args will (calling object, other, some_more_args_maybe)

    """

    caller = args[0]
    other = args[1]
      

    if isinstance(other,Field) or isinstance(other,VField):
      return other.vcumsum(coord = caller)

    # otherwise go about normal business:
    return func(*args, **kwargs )

  return dispatcher









# ----- most general classes ------------

class Named(object):
  """
  Base class for most other sg classes, representing objects with copy and same methods.

  The "same" method indicates when Named objects are "the same", namely when their "name" attribute is the same. This method coincides with the "weaksame" method. "Weaksame" is generally a weaker condition in the derived classes. The "same" method allows the implementation of the "samein" and "sameindex" methods at this abstract level, with generally the "same" method overriden in derived classes.

  This class provides a copy method that is used by the derived classes.

  Attributes:
      name: (str) name of Object
  """



  def __init__(self,name='scalar' ,long_name= ''):  
    """
    Initialisation of Name object. 

    Args:
      name: (str) name of Object
      long_name: (str) longer description (e.g. for display)


    Returns:
      Named object
    """
  
    self.name = name
    self.long_name = long_name


  def copy(self, equiv = True, *args, **kwargs):
    """
    Copy method for Named. See __init__ for arguments.

    Most child classes should inherit this method.

    Args (in addition to __init__):
      equivs: (Boolean) If True, the copies will be equivalent. 

    Returns: 
      a copy of the Directional object.

    Copy methods in sg work as follows: when no value is selected for an argument, a copy of the self attribute will be used. Otherwise, the **kwargs argument value will be used.
    """

    # keys to exclude:    
    forbid = ['self','frame']

    # keys value dict to examine:
    allow = {key:self.__dict__[key] for key in inspect.getargspec(self.__init__)[0] if key not in forbid  }

    # new kwargs to use in new object .__init__ construction
    new_kwargs = copy.deepcopy({})

    for key in allow:
      if key in kwargs:
        # if given in kwargs, override
        new_kwargs[key] = kwargs[key]
      else:
        # otherwise use value in self attribute
        new_kwargs[key] = self.__dict__[key]  

    new_kwargs = self._copy_cleanup(**new_kwargs)

    # initialize new object:
    result = self.__class__(**new_kwargs)

    if equiv:
      # decide whether to copy over the equivs list
      if hasattr(self,'equivs'):
        result.equivs = self.equivs

    return result

  def _copy_cleanup(self, **new_kwargs):
    """
    Override this method for specific work to be done before new object init in copy method (e.g. duals in Coord).
    """
    return new_kwargs



  def __and__(self,other):
    """Shorthand for weaksame method. See weaksame.
    """
    return self.weaksame(other = other)



  def weaksame(self,other):
    """
    Tests if two Directional objects have the same name.

    Weak test to see if two Directional objects are similar.

    Args:
         other: other Directional object to compare self with
       
    Returns: 
         True/ False

    **See also**
    same method
    samein method
    sameindex method
    """

    if (self.name == other.name):
      return True
    else:
      return False


  def same(self,other):
    """Method to check whether this Named object has identical main attributes to argument other. 

    Placeholder identical to weaksame: to be overriden in child classes.

    Args:
      other: (Named) to check against

    Returns:
      True/ False

    Attributes checked:
      name: via str ==
 
    See also:
      samein method
      same_index method
    """

    return self.weaksame(other)




  def samein(self,L):
    """Tests whether this Directional is the same as any element in list L, under 'same' method.

    Uses: same method.

    Args:
      L: (list of Directional objects) to test against

    returns:
      True/ False

    See also:
      same method
      same_index method
    """
    return reduce(lambda x,y:x or y, [self.same(it) for it in L] )

  def sameindex(self,L):
    """Find index of this Directional in list L of Directional objects, under 'same' method.

    Uses: same method.

    Args:
      L: (list of directional objects) to search

    returns:
      None or Integer, the index of the first item it in the list that satisfies self.same(it)

    See also:
      same method
      samein method
    """

    for i,it in enumerate(L):
      if self.same(it):
        return i
     
    return None    








class Directional(Named):
  """
  Base class for derived Coord and Ax classes, representing "direction" (e.g. "latitude" or "depth").

  An abstract equivalence relationship is defined among directional objects by adding them to each other's equiv attribute (a list of equivalent Directionals). This relationship is generally used to indicate whether two Directional objects have the same direction (e.g. X,Y), but could represent other relationships depending on the user.

  The same method is differentiated from the weaksame method (unlike the parent class), with the more strict additional condition that in addition to "name", the "direction" attribute also needs to be the same. therefore, two Directional objects are considered "same" when both "name" and "direction" match. They are "weaksame" when only the "name" matches.

  This base class is closely related to the Ax class. The Coord class is also derived from it.


  Attributes:
      name: (str) name of Object
      direction: (str) name of direction in which object points
      long_name: (str) longer description (e.g. for display or in Netcdf)
  """

  def __repr__(self):
    """Display alias if attribute present, name otherwise.
    """
    if hasattr(self,'alias'):
      return self.alias
    else:
      return self.name

  def __init__(self,name='scalar',direction ='scalar',long_name= '' ):  

    """
    Initialisation of Directional object. 

    Args:
      name: (str) name of Object
      direction: (str) name of direction in which object points
      long_name: (str) longer description (e.g. for display or in Netcdf)

    Returns:
      Directional object
    """
  
# choosing the name ID creates an identity object. ID*b = b for all Coord elements b.
# could implement the identity Field in __call__

# Metric could be a class. Objects of this class could be constructed by a method of the Coord class (Coord objects then spawn metric objects).
    self.equivs = [self]
    self.name = name
    self.direction = direction
    self.long_name = long_name

  def same(self,other):
    """Method to check whether this Directional object has identical main attributes to argument other. 

    Overrides Name class same method.    

    Args:
      other: (Directional) to check against

    Returns:
      True/ False

    Attributes checked:
      name: via str ==
      direction: via str == 

    See also:
      samein method
      same_index method
    """

    return (self.name == other.name) and (self.direction == other.direction)




# ----- equivalence related ---------  

  def make_equiv(self,other):
    """
    Register equivalence of two Coord objects. 

    Args:
      other: (Coord)

    Returns:
      None

    See also:
      is_equiv
      eq_index
      eq_in

    Examples:
    >>> depth.is_equiv(longitude) # generally different directions.
    False
    >>> depth.make_equiv(longitude) # don't do this in real work
    >>> depth.is_equiv(longitude) # uphysically:
    True 
    """  

#      self.equivs.append(other)
#      other.equivs.append(self)    

    for e in set(self.equivs):
      e.equivs.append(other)

    for e in set(other.equivs):
      e.equivs.append(self)        

    self.equivs = list(set(self.equivs))
    other.equivs = list(set(other.equivs))

    return


  def is_equiv(self,other):
    """
    Test for equivalence (under make_equiv). e.g. xt is equivalent to xu

    Args:
      other: (Coord or Ax)

    Returns:
      True when equivalent, False otherwise.

    Examples:
    >>> depth.is_equiv(longitude) # generally different directions.
    False

    See also:
      make_equiv
      eq_index
      eq_in
    """
    
    if (other in self.equivs) | (self in other.equivs):

      return True     
    else:
      # Warnings helpful for debugging and spotting potential problems:

      if ( self.same(other) ):
        warnings.warn('Warning (severe): %s.is_equiv(%s) is False, but %s.same(%s) is True! ' % (self,other,self,other) )
      elif ( self.weaksame(other) ):
        warnings.warn('Warning (severe): %s.is_equiv(%s) is False, but %s.weaksame(%s) is True! ' % (self,other,self,other) )

      # then go about normal business:
      return False



  def eq_index(self, collection ):
    """ Find index of first element in collection equivalent to self under is_equiv.

    Args:
      collection: (e.g. List) order collection of Directional objects

    Returns:
      Integer if equivalent found (index to equiv element), None otherwise.


    See also:
      make_equiv
      is_equiv
      eq_in
    """

    for i, e in enumerate(collection):
      if self.is_equiv(e):
        return i
    return 


  def eq_in(self, collection):

    """ Determines whether Coord (self) is equivalent to any of the constituent Coord objects of the argument Gr or GrAx, and returns equivalent object.

    Uses: eq_index

    Args:
      collection: (Gr or AxGr) object to be checked.

    Returns:
      The equivalent object  when crd is equivalent to one of the Coord objects in argument, None otherwise.

    See also:
      eq_in method of Ax, GrAx
      is_equiv
      make_equiv
      eq_index
      eq_in   
    """


    i = self.eq_index(collection)
    if i is not None:
      return collection[i ]
    else:   
      return




# Belongs to Directional

  @field2cumsum
  def __or__(self,other):
    """
    Shorthand calling make_equiv (to register equivalence with other Directional object). 

    Args:
      other: (Directional)

    Returns:
      None
    """

    self.make_equiv(other)


  def __xor__(self,other):
    """
    Shorthand calling is_equiv (to test equivalence with other Directional object). 

    Args:
      other: (Directional)

    Returns:
      None
    """

    return self.is_equiv(other)



# ----- multiplication related ---------      

  def __pow__(self,n):
    """
    Repeated multiplication of object with itself.
    """
    return reduce(lambda x,y: x*y, n*[self])


class Membered(Named):
  """

  """


  def weaksame(self,other):
    """
    Member-wise weaksame comparison.
    """

    if len(self) == len(other):
      for i,c in enumerate(self):
        if not(c.weaksame(other[i]) ):    
          return False

      return True

    else:

      return False




class Valued(Named):
  """
  Base class for classes that contain a ndarray value attribute. 
  """

  def __init__(self,name='scalar',value = np.array([0]),long_name =''):  

    """
    Initialisation of Valued object. 

    Args:
      name: (str) name of object.
      value: (Numpy ndarray) 
      long_name: (str) a longer description. 
    """
  
    self.name = name
    self.value = value
    self.long_name = long_name
    self.shape = value.shape


  def __repr__(self):
    """Display alias if attribute present, name otherwise.
    """
    if hasattr(self,'alias'):
      return self.alias
    else:
      return self.name

  def __getitem__(self,i):
    """Obtain item from value atttribute.
    """
    return self.value[i]

  def array_equal(self,other):
    """ test whether Valued objects contain identically valued ndarrays in value attributes.

    This is a common method that should be inherited by child classes. 

    Args:
      other: (Valued object) the Valued to compare with

    Returns:
      True/ False (using np.array_equal)

    Raises:
      ValueError: when argument is not a Coord object
    """

    if not isinstance(other,Valued):
      raise TypeError('Error: provide Valued argument (%s provided).'%other)

    return np.array_equal(self.value,other.value) 


  def same(self,other):
    """
    Tests whether this Valued object contains identical name and value to argument object. 

    Overrides Named same method and is a stronger condition. Generally to be overriden in child classes.

    Args:
      other: (Valued) object to compare against.

    Returns: 
      True/ False
    """

    # the following test and warning is to do with little things that we don't want to trip over with errors:
    # self is a Coord, so it has a value attribute (this should be put into the specs!), but other might not:
    if not hasattr(other, 'value'):
      warnings.warn('Valued method %s.same(%s) on argument without Ax attribute: returning False.'%(self.name, other.name))      
      return False

    return (self.name == other.name) and self.array_equal(other) 



  def __neg__(self):
    """
    Obtain version of Valued with negative values (e.g. -xt). 

    Returns:
      Valued copy with value is -self.value
    """

    return self.copy(value =-self.value)


  def __add__(self, other):
    """ Addition of value attributes
    """
    return self.copy(value = self.value + other.value)

  def __sub__(self, other):
    """ Substraction of value attributes
    """
    return self.copy(value = self.value - other.value)


  def __mul__(self, other):
    """ Multiplication of value attributes
    """
    return self.copy(value = self.value*other.value)

  def __div__(self, other):
    """ Division of value attributes
    """
    return self.copy(value = self.value/other.value)



# -------- Coord  ----------------

class Coord(Directional, Valued):
  """
  Representing distrete coordinate collection. Corresponds to dimension variable in Netcdf.

  Coord objects are defined by a name, value (generally a numpy array) and units. The value of a Coord is a 1D ndarray containing the locations of data points. Coord is the basic building block of Gr (grid) objects.
  Examples of Coord objects are xt or yt, corresponding to the tracer grid cells in the x and y directions. Coord objects have a corresponding dual. For xt it is xt_edges and vice versa. The dual generally contains the edges of the grid cells. If no dual argument is given, the Coord object is its own dual.

  Coord objects c1 and c2 are considered equal, c1&c2 yields True, when the name, value (numpy array) and units attributes are equal.

  Being a container of Coord objects, the Gr object (grid) is closely related to Coord. A shorthand for Gr construction is via multiplication of Coord objects. If two Coord objects coord1 and coord2 are not equivalent (generally when they point in different directions, e.g. X and Y), their product is a shorthand for the creation of a 2D grid coord1*coord2 = Gr((coord1, coord2)). By induction, products containing n elements yield Gr objects of dimension <=n. See class documentation.

  The Coord class contains methods related to distances, with the following dependencies. dist is fundamental. delta_dist depends on dist. der depends on dist. d depends on dual Coord and delta_dist. vol depends on d.


  Attributes:
    axis: (Ax) axis along which Coord is defined.
    dual: (Coord) link to other Coord object, usually representing edges of cells.
    direction: (str) name of direction in which Coord is defined.
    equivs: (list of Coords) list of links to equivalent Coord objects. Usually coincides with all Coord's defined in same direction.
    len: (int) length of the value attribute (ndarray)
    long_name: (str) a longer description. Usually taken from identically named variable in Netcdf file.
    metadata: (dict) name vs value dict of metadata to save into Netcdf. Usually obtained from Netcdf read.
    name: (str) name of object (e.g. 'xt' or 'depth').
    nbytes: (int) approximate memory usage of Coord object.
    strings: (list of strings or None) str labels for coord points. If not None, needs to be of same length as value
    units: (str) units (e.g. "meter")
    value: (Numpy ndarray) 1D array corresponding to spatial points (e.g. 10 degrees S for a point).

  Examples:    
    >>> coord1 = sg.fieldcls.Coord(name = 'test1',value =np.array([1.,2.,3.]) ,axis =sg.fieldcls.Ax('X'),direction ='X', metadata = {'hi':5} )
    >>> coord2 = sg.fieldcls.Coord(name = 'test2',value =np.array([1.,2.,3.,4.]),axis =sg.fieldcls.Ax('Y'),direction ='Y', metadata = {'hi':7})
    >>> coord3 = sg.fieldcls.Coord(name = 'test3',value =np.array([1.5,2.5,3.5,4.5]),axis =sg.fieldcls.Ax('X'),direction ='X', metadata = {'hi':7})
  """


  def same(self,other):
    """
    Coord method to check whether this Coord has identical main attributes (except units) to argument other Coord.

    The axis attributes may sometimes be a str. In this case, == is applied and a warning is issued.

    Args:
      other: (Coord) Coord to check against

    Returns:
      True/ False

    Attributes checked:
      value: via Coord array_equal method
      name: via str ==
      axis: same method if Ax object, str == otherwise (axis attribute can sometimes str by choice, although this is not great)
      direction: via str == method

    See also:
      samein method
      same_index method
    """

    # not checking for units

    # the following test and warning is to do with little things that we don't want to trip over with errors:
    # self is a Coord, so it has an axis attribute (this should be put into the specs!), but other might not:
    if not hasattr(other, 'axis'):
      warnings.warn('Coords method %s.same(%s) on argument without Ax attribute: returning False.'%(self.name, other.name))      
      return False

    if (isinstance(self.axis,str) or isinstance(self.axis,unicode)) and (isinstance(other.axis,str) or isinstance(other.axis,unicode)):
      # Both of the axis attributes are a str. Comparing apples with apples.

      warnings.warn('Coords %s, %s have str axis attribute %s, %s'%(self.name, other.name,self.axis, other.axis))

      return self.array_equal(other) and (self.name == other.name) and (self.axis == other.axis ) and (self.direction == other.direction  )

    elif (isinstance(self.axis,str) or isinstance(self.axis,unicode)) or (isinstance(other.axis,str) or isinstance(other.axis,unicode)):
      # only one of the axis attributes is a str. The other would be an Ax. Apples and oranges: tricky, but should be uncommon.

      warnings.warn('!! One of Coords %s, %s has str axis, not the other! Proceeding using str compare. (%s, %s)'%(self.name, other.name,self.axis, other.axis))


      return self.array_equal(other) and (self.name == other.name) and (str(self.axis) == str(other.axis) ) and (self.direction == other.direction  )


    else:
      # Both of the axis attributes are an Ax. Comparing apples with apples. Best case of 3.
      return self.array_equal(other) and (self.name == other.name) and (self.axis.same(other.axis) ) and (self.direction == other.direction  )


  def weaksame(self,other):
    """
    Tests whether Coord objects contain identical Coord values, name and direction. 

    Overrides Directional weaksame method and is a stronger condition.

    Args:
      other: (Coord) object to compare against.

    Returns: 
      True/ False
    """

    # the following test and warning is to do with little things that we don't want to trip over with errors:
    # self is a Coord, so it has a value attribute (this should be put into the specs!), but other might not:
    if not hasattr(other, 'value'):
      warnings.warn('Coords method %s.weaksame(%s) on argument without Ax attribute: returning False.'%(self.name, other.name))      
      return False

    return (self.name == other.name) and (self.direction == other.direction) and self.array_equal(other) 




    
  def sort(self,*args,**kwargs):
    """
    Sorts Coord value. Passes arguments on to Numpy sort method.
    """

    self.value.sort(*args,**kwargs)


  def __init__(self,name='scalar',value = np.array([0]), dual = None,axis = '?',direction ='scalar', units = None,long_name ='?', metadata = {} , strings = None):  

    """
    Initialisation of Coord object. 

    If argument dual is provided, the dual attribute of that Coord object is also set to self (the Coord objects become each other's dual).

    Args:
      name: (str) name of object (e.g. 'xt' or 'depth').
      value: (Numpy ndarray) 1D array corresponding to spatial points (e.g. 10 degrees S for a point).
      dual: (Coord) link to other Coord object, usually representing edges of cells.
      axis: (Ax) axis along which Coord is defined.
      direction: (str) name of direction in which Coord is defined.
      units: (str) units (e.g. "meter")
      long_name: (str) a longer description. Usually taken from identically named variable in Netcdf file.
      metadata: (dict) name vs value dict of metadata to save into Netcdf. Usually obtained from Netcdf read.
      strings: (list of strings or None) str labels for coord points. If not None, needs to be of same length as value

    Raises:
      ValueError: if argument strings is not None and if len(strings) != len(value)
    """
  
# choosing the name ID creates an identity object. ID*b = b for all Coord elements b.
# could implement the identity Field in __call__

# Metric could be a class. Objects of this class could be constructed by a method of the Coord class (Coord objects then spawn metric objects).
    self.equivs = [self]
    self.name = name
    self.value = value
    self.axis = axis
    self.units = units
 
    self.direction = direction
    self.long_name = long_name

    self.fld = None
    self.len = len(value)
    self.metadata = metadata
    self.nbytes = self.value.nbytes

    if (not axis is None) and isinstance(axis, Ax)  :
      (self.equivs).append(axis)

    if strings is not None:
      if len(value) != len(strings):
        raise ValueError('Provide strings argument of equal length to value argument for Coord %s if providing strings argument. %s vs %s'%(name, len(value) , len(strings)))
    self.strings = strings

    self.give_dual(dual)
 

  def _copy_cleanup(self,**new_kwargs):
    """
    When a Coord is copied and is self-dual, the new copy must be dual to itself, not the original calling Coord that is being copied. 

    This method is called by the copy method of the parent class. If this method were not called, the copy of a self-dual Coord would be dual to its original calling Coord rather than itself.

    Note that this method works only on self-dual coords: copies of other-dual Coord objects remain dual to that existing dual Coord.
    """

    if 'dual' in new_kwargs:

      if new_kwargs['dual'] is self:
        new_kwargs['dual'] = None

    return new_kwargs

  def give_axis(self, axis):
    """
    Give axis attribute value and set equivalence to axis.

    Args:
      axis: (Ax) the axis to provide
    """
    self.axis = axis
    self.make_equiv(axis)


  def give_dual(self,dual = None):
    """
    Provide Coord with dual Coord (e.g. latitude vs latitude_edges). 

    If Coord is provided, dual attribute will be set to that Coord and vice versa. Called upon initialization, and therefore by copy method.

    Args:
      dual: (Coord or None) Coord to make dual 
    """

    if dual is None:
      self.dual = self
    else:
      self.dual =dual
      dual.dual = self


  def __len__(self):
    """Obtain length of value attribute (1D ndarray)
    """
    return self.len

# belongs to Coord 
  def __call__(self, other = None, index = 0):
    """Shorthand for cast. See cast method.
    """
    return self.cast(other=other, index = index)

  def cast(self, other = None, index = 0):
    """
    Broadcasts Coord onto grid.

    Calling a Coord object with a grid Gr object as argument yields an array with the coord values broadcast onto that grid.

    In other words, the resulting array A is defined on that grid where A[:,i,...] = self.value[i] for all i. This leads to an expansion of the Coord value useful for grid operations such as interpolation.

    Args:
      other: (Gr, Field, None) Gr to broadcast on or Field to slice.
      index: (int) index at which to take a slice in case of Field argument.

    Returns:
       Field with Coord broadcast onto that grid if Coord is in that Gr. If not, None is returned.
       
    Examples:
      >>> R = xt(xt*yt*zt)   # obtain Field R
      >>> R.shape == (len(xt),len(yt),len(zt))
      True
      >>> R.grid
      (xt, yt, zt)
      >>> # Here, the value of R is constant in yt and zt, but equal to xt along the xt axis.
    """

#    if not(other) or (other is self):
#      if not(self.fld):
#        self.fld = Field(name = self.name, value = self.value, grid = self**2, units = self.units)
    
#      return self.fld
#    else:

    if isinstance(other,Gr):
      if self in other:
        return Field(name = self.name, value = (self**2)(other)(self.value), grid = other)
      else:
        return 
#      elif isinstance(other,Field):
#          return self._bigslice(other,index)


# ----- addition related --------- 

  def __add__(self,other):

    """
    Refine Coord by combining grid points from both. Only implemented for self-dual Coord objects.
    """

    if ((self.dual == self) and (other.dual == other)):

      result = self.copy()

      result.value = merge(self.value,other.value)
      
      return result

    else:

      print '+ only implemented for self-dual Coord objects (e.g. time), returning None.'
      return None  


      
  def __neg__(self):
    """
    Obtain version of Coord with negative values (e.g. -xt). Includes a negative dual (edges). 

    Returns:
      Coord copy with value is -self.value

    The dual is also made negative.
    """

    neg_crd = self.copy(value =-self.value)

    if self.dual is self:
      neg_crd.give_dual()
   
    else:
     
      neg_dual = self.dual.copy(value =-self.dual.value)
      neg_crd.give_dual(neg_dual)

    return neg_crd

  def __xor__(self,other):
    """
    Shorthand for several methods depending on argument type. 

    Becomes is_equiv if argument Coord, Gr.eq_in on self if argument is Gr and self.der on argument for Field argument.

    Args:
      other: (Coord, Gr or Field)

    Returns
      None unless other is Field

    Raises: TypeError
    """


    if isinstance(other,Coord) | isinstance(other,Ax):
      return self.is_equiv(other)
      
    elif isinstance(other,Gr):
      if other.eq_in(self):
        return other[other.eq_index(self)  ]
      else:   
        return

    elif isinstance(other,Field):
      return self.der(other)
    else:
      raise TypeError('Coord error in %s^%s with Coord %s. Provide Coord,Gr or Field object for right multiplicant (now %s).' % (self,other,self,other) )
      return

  def eq_index(self,grid):
    if grid.eq_in(self):
      return grid.eq_index(self)  
    else:   
      return




# ----- multiplication related ---------      


# --> belongs to Coord 
  def __mul__(self,other):
    """
    Multiplication of Coord with Coord, Ax, Gr, (V)Field.
    
    A shorthand for grid (Gr object) construction is via multiplication of Coord objects. If two Coord objects coord1 and coord2 are not equivalent (generally when they point in different directions, e.g. X and Y), their product is a shorthand for the creation of a 2D grid coord1*coord2 = Gr((coord1, coord2)). If coord1 and coord2 are equivalent (point in the same direction), coord1*coord2 yields Gr((coord1,)). By induction, products containing n elements yield Gr objects of dimension <=n. See class documentation.

    In case the right multiplicant (argument) is a (V)Field, the product yields the zonal integral of that (V)Field.

    Args:
      other: (Coord, Ax, Gr or (V)Field)

    Returns
      Gr if argument Coord/ Gr 
      None, Coord if argument Ax (so that Ax multiplication commutes)
      (V)Field if argument field (zonal integral)

    Raises: TypeError    
    """

    if isinstance(other,Coord):   
      if other.is_equiv(self):
        return Gr((self,))
      else:
        return Gr((self,other))

    elif isinstance(other,Ax):
      # Ax and Coord objects commute
      return other*self

    elif isinstance(other,Gr):   
      if other.eq_in(self):
        new_other = list(other)
       
        new_other[other.eq_index(self)] = self
        return Gr(new_other)
      else:
        return Gr( [self] + list(other))
    elif isinstance(other,Field):
#        print 'Warning (benign): converting left multiplicant to Gr object from Coord object.'
      return (self**2)*other

    elif isinstance(other,VField):
       # case whether vector Field is multiplied by Coord (yielding vsum).
       # this commutes:
       return other*self

    else:
      raise TypeError('Coord error in %s*%s with Coord %s: provide Coord, Gr or Field object as right multiplicant (now %s). If multiplicant appears to be a Coord of other multiplicant Field, check whether its definition is stale (reloaded sg since its creation). '% (self,other,self,other) )




  def start_zero(self):
    """
    Returns a copy of this Coord where the coordinate values start at 0.

    This is achieved by subtraction of self.value[0] from self.value
    """
    return self.copy(name = self.name + '_zero'  , value = self.value - self.value[0])
 
  def cdf_insert(self,file_handle, miss_default = 9.96921e+36):
    """
    Netcdf insert method of Coord 
    Inserts Coord as variable into Netcdf file.

    Input: file_handle file handle of opened Netcdf file.

    Inserts Coord as variable into Netcdf file.

    Args:
      file_handle: (file object). Points to opened file.
      miss_default: (float) value to use as default for missing values/ NaNs

    Returns:
      file_handle
    """

    # make a copy of self content and deal with missing values.
    value = copy.deepcopy(self[:])

    miss_val = miss_default
    if 'FillValue' in self.metadata:
      miss_val = self.metadata['FillValue']
    elif 'missing_value' in self.metadata:
      miss_val = self.metadata['missing_value']

    try:
    
      value[  np.isnan( value[:]  ) == True   ] = miss_val
    except:
      try:
        value[  np.isnan( value[:]  ) == True   ] = miss_default

      except:
        warnings.warn('Could not set missing value for Coord %s.'%self.name)

    file_handle.createDimension(self.name,len(self))


    var_cdf = file_handle.createVariable(self.name, value.dtype.char, (self.name,)   )
    
    for k in self.metadata:
      setattr(var_cdf,k, self.metadata[k]) 

    var_cdf[:] = value

    return file_handle


  def write(self, path = None, name = None , history = 'Created from Spacegrids '  ):
    """
    Writes Coord data to Netcdf file.


    Args:
      path: (str) path to the directory where the file will go (e.g. 'data/' or '/home/me/', default pwd).
      name: (str) file name (e.g. "foo.nc")
      history: (str) Brief history or general description of the data.
      insert_dual: (Boolean) Flag determining whether to include the duals of the Coord objects in the file.


    Returns:
      None
   
    No error is raised in case of problems (only a message is displayed).
    """

    if name is None:
      name = self.name

    if not name.split('.')[-1] in ['nc','cdf']:
      name = name +'.nc'
    if not path is None:
      name = os.path.join( path , name ) 

    print 'Writing Field to file %s'%name

    try:
      # using io_sg version of netcdf_file. May use ScientificIO or Scipy
      file_handle = netcdf_file(name , 'w')
    except IOError:
      print 'Cannot write ', name
    else:

      file_handle = self.cdf_insert(file_handle)

      file_handle.history = history + '%s'%str(datetime.datetime.now())
 
#    var_cdf.units = self.units

      file_handle.close()



  def finer(self,factor = 5.):
    """
    Method of Coord. Refine the coordinate point interval by a given factor.

    Args:
      factor: (float) factor by which to refined Coord value.

    Returns:
      Coord with value refined according to factor.
    """

    result = []

    for i in range(0,len(self)-1):
      result += list(np.arange(self[i],self[i+1],(self[i+1] - self[i])/factor))


    finer_coord = self.copy(name = self.name + '_fine',value = np.array(result))  
    finer_coord|self
    return finer_coord


  def _bigslice(self, F = None, index = 0):
    """
   Coord method that takes slice of field along Coord (self) at argument index.    

    Args:
      F: (Field) field to slice
      index: (int) location of slice


    Returns:
      Field representing sliced object, with self no longer appearing in Gr.
    """

    if not(F):
      warnings.warn('Warning Coord part: Provide Field.')

    if self in F.grid:
      sl = slice(None,None,None)
      sl2 = slice(index,index+1,1)
      L = []
      for e in F.grid:
        if self is e:
          L.append(sl2)
        else:
          L.append(sl)

      return F.copy(name = F.name , value = np.squeeze(F.value[L]), grid = F.grid/self, units = F.units)
   
    else:
      warnings.warn('Warning Coord slice: Coord not in Field grid. Returning Field.' )   
      return F

# belongs to  Coord
  def coord_shift(self,F,shift, keepgrid = False, nan_val = np.nan):
    """
    Coord method that shifts the coordinates and value of a field by a number of indices. 

    The shifted Coord in the grid of the Field argument is replaced with a (different) shifted Coord: disable this behaviour with argument keepgrid = True. Calls roll function. 

    Args:
      F: (Field) field to shift
      shift: (int) magnitude (corresponding to array index) of shift
      keepgrid: (Boolean, default False) grid not updated if True

    Returns:
      The shifted Field.
    """

    return roll(F,shift = shift,coord = self,mask=True, keepgrid = keepgrid, nan_val = nan_val)

# belongs to  Coord 
  def trans(self,F):
    """
    Gives the change in Field F upon a Coord shift of 1 index in the direction of the self Coord.

    F - self.coord_shift(F,shift=1,keepgrid = True)

    Args:
      F: (Field) field to act on

    Returns:
      Field: the transformed Field.
    """

    # select keepgrid = True to avoid substraction errors relating to different grids
    return F - self.coord_shift(F,shift=1,keepgrid = True)

  def sum(self,F, land_nan = True):
    """
    Method of Coord  that sums Field F along self Coord direction. Not weighted with grid cell width. Uses masked arrays to handle nan values. nan values can be used to eliminate areas from summing area.


    Args:
      F: (Field) field of certain dimension n to sum
      land_nan: (Boolean)

    Returns:
      Field of dimension n-1 or float if n=1. 

    Raises:
      ValueError: when Coord (self) is not in F.grid (grid of Field).
    """

    if not self in F.grid:
      raise  ValueError('Coord sum method of %s: Coord must be in grid of argument Field. Make sure Coord object is identical to one of Coord objects in Field grid. (Also watch for stale objects.)'%self.name)   

    value = np.array(ma.sum(ma.masked_array(F[:],np.isnan(F[:])), axis = (F.grid).index(self) ))
    find_land = np.array(ma.sum(ma.masked_array(np.ones(F[:].shape),np.isnan(F[:])), axis = (F.grid).index(self) ))

    if land_nan:
      if not value.ndim == 0:
        value[find_land == 0.] = np.nan


    if len(F.grid) == 1:
      return float(value)
    else: 
      return F.copy(name = F.name,value = value, grid = F.grid/self )
 

# ---> belongs to Coord class

  def roll(self,shift = 0):
    """Yields copy of Coord (self) with value shifted by (int) argument shift using numpy.roll().

    Args:
      shift: (int) number of array positions to shift by
 
    Returns:
      Coord: shifted Coord.
    """

    return self.copy(name = self.name + '_rolled',value = np.roll(self.value,shift = shift))

  def flip(self,F):
    """
    Reverse order of elements along axis of this Coord (a mirror, or flip). 

    Grid remains unchanged: strictly, this will lead to an inconsistency between the Field data and the grid, assuming this is what the user wants.


    Args:
      F: (Field) Field to mirror
 
    Returns:
      Field of equal dimension, mirrored along Coord.
    """
  
    SI = self._slice_index(F.grid, slice_obj = slice(None,None,-1))

    return F.copy(name = F.name,value = F[SI],grid = F.grid, units = F.units)



  def _slice_index(self,grid , slice_obj = slice(1,None,None)):
    """
    Yields a list of slice objects that can be used to slice along the axis of this (self) Coord.


    Args:
      grid: (Gr) grid context of slicing
      slice_obj: Slice object 

    Returns:
      List of slice objects of equal length to (argument) grid. These objects can be used to slice Field objects defined on that grid.


    Examples:

    >>> coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) )
    >>> coord2 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]) )
    >>> K = coord1(coord1*coord2)
    >>> coord1.slice_index(K)
    [slice(1, None, None), slice(None, None, None)]
    """
    sl = slice(*(None,))
    L = []
    for e in grid:
      if self is e:
        L.append(slice_obj)
      else:
        L.append(sl)

    return L

# --> belongs to Coord 

  def cumsum(self,F, upward = False,land_nan = True):
    
    """
    Compute cumulative sum (integral) of input Field F along axis of F corresponding to this Coord object.

     If argument upward is set to true, summing takes place with increasing array index. If it is set to False, summing takes place with decreasing array index starting at index -1. Values of nan are set to 0, and therefore not counted.

    Args:
      F: (Field) Field to sum
      upward: (Boolean) flag to set direction of cumsum
      land_nan: (Boolean) flag to set land to nan in the resulting array

    Returns:
      Field on same grid containing the cumsum.

    Raises:
      ValueError: when Coord (self) is not in F.grid (grid of Field).

    Examples:

    >>> coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) )
    >>> coord2 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]) )
    >>> R = coord1.cumsum(sg.ones(coord1*coord2)  );R.value
    array([[ 3.,  3.,  3.,  3.],
           [ 2.,  2.,  2.,  2.],
           [ 1.,  1.,  1.,  1.]])

    >>> R = coord1.cumsum(sg.ones(coord1*coord2) , upward = True );R.value
    array([[ 1.,  1.,  1.,  1.],
           [ 2.,  2.,  2.,  2.],
           [ 3.,  3.,  3.,  3.]])
    """

# nan values are set to 0. They are not counted.
    if not self in F.grid:

      raise ValueError('Coord %s cum_sum method: Coord must be in grid of argument Field %s. Make sure Coord object is identical to Coord objects in Field grid.'%(self.name, F.name)  )


    if upward:
      # use the copy method of the Field to obtain a similar Field, but with some attributes different (namely, those set in the argument).
      Fc = F.copy()
        
    else:
      Fc = self.flip(F.copy() )
    
    # do not count land in sums:  
    land_i = np.isnan(Fc[:]) 
    Fc[land_i] = 0.

    result_array = np.array(np.cumsum(Fc[:], axis = (Fc.grid).index(self) ))

    if land_nan == True:
      # if desired, set land to nan in resulting array:
      result_array[land_i] = np.nan

    if upward:
      return F.copy(value = result_array )
    else:
      # need to flip back:
      return self.flip(F.copy(value = result_array))



  def vsum(self,F):
    """
    Method of Coord .
    Sums Field along self Coord, weighted with grid cell width (using self.d(), called by self.vol(F.grid)). Note: due to possible dependence of one Coord on the other, only use mean method of grid. There is no mean method for Coord objects.

    Method of Coord  that sums Field F along self Coord direction, weighted with the grid cell width. Calls sum method. See sum method.

    Calculation is self.sum(F*(self.vol(F.grid))). For grids with grid cell width depending on coordinates, use corresponding Gr methods.     

    Args:
      F: (Field) field of certain dimension n to sum
      
    Returns:
      Field of dimension n-1 or float if n=1. 

    Raises:
      ValueError: when Coord (self) is not in F.grid (grid of Field).
    """

    return self.sum(F*(self.vol(F.grid)))     



  def vcumsum(self,F,upward =True):
    """
    Compute cumulative sum, weighted with grid cell width, of input Field F along axis of F corresponding to this Coord object.

     If argument upward is set to true, summing takes place with increasing array index. If it is set to False, summing takes place with decreasing array index starting at index -1. Values of nan are set to 0, and therefore not counted. Calls cumsum method. See cumsum.  For grids with grid cell width depending on coordinates, use corresponding Gr methods.     

    Calculation is self.cumsum(F*(self.vol(F.grid)) )   

    Args:
      F: (Field) Field to sum
      upward: (Boolean) flag to set direction of cumsum
      land_nan: (Boolean) flag to set land to nan in the resulting array

    Returns:
      Field on same grid containing the cumsum.

    Raises:
      ValueError: when Coord (self) is not in F.grid (grid of Field).
    """

    return self.cumsum(F*(self.vol(F.grid)) , upward = upward)   

# belongs to  Coord
  def dist(self, fact = 1.):
    """
    Distance (signed) along Coord from a certain fixed point (e.g. from ocean surface or from equator along y-direction).


    Args:
      fact: (float) magnification factor if required (e.g. radius of Earth).

    Returns:
      Field on rid self**2 containing factor*self.value as value. Represents distances of points from a certain point along Coord.

    See also:
      d method 
      delta_dist method
      der method
      vol method      
    """

    return Field(name='distance_'+self.name,value = fact*self[:],grid = (self**2), units = self.units) 


  def der(self, F):

    """
    Coord derivative method on Field F. 

    If Coord non-cyclical, the first derivative element is nan and the second is the derivative at the first element of the original Coord. 


    Args:
      F: (Field) Field to take derivative of

    Returns:
      Field on same grid containing the derivative.

    Raises:
      ValueError: when Coord (self) is not in F.grid (grid of Field).

    See also:
      d method 
      delta_dist method
      dist method 
      vol method      
    """

    if self in F.grid:
      dF = self.trans(F)
      ds = self.trans(self.dist())

      return dF/ds

    else:

      raise  ValueError("Field argument grid does not contain Coord %s."%self.name)


# belongs to  Coord  
  def delta_dist(self, fact = 1.):
    """
    Method to calculate the distance between adjacent elements of Coord.
    Appropriate to vertical direction.
    To be over-ridden for hor coords x,y => derive classes XCoord, YCoord

    Calls dist method and applies trans method.

    Returns an array as len(result) == len(grid)-1

    Calculates self.trans(self.dist())*fact

    Args:
      fact: (float) magnification factor if required (e.g. radius of Earth).

    Returns:
      Field: containing the distances between the adjacent coord points (i.e. i and i+1).

    See also:
      d method 
      der method
      dist method 
      vol method      

    Examples:

    >>> coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) )
    >>> R = coord1.delta_dist();R.value
    array([ nan,   1.,   1.])
    """
       
    return self.trans(self.dist())*fact

# --> belongs to  Coord
  def d(self):
    """
    Calculates width of grid cell in direction of Coord (self) using the dual of self (e.g. zt_edges). 

    Yields grid cell widths. Can be used to compute volumes.

    To be overriden in x and y direction to accomodate for sphere.

    Calculates self.dual.delta_dist() where it is defined

    Returns:
      Field: containing the distances between the adjacent coord cell edges.

    See also: 
      delta_dist method
      der method
      dist method 
      vol method      

    Examples:

    >>> coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]))
    >>> coord1_edges = sg.fieldcls.Coord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = coord1 ) # specifying the dual here registers coord1_edges also as the dual attribute of coord1
    >>> coord1.d().value    # the distance between the cell edges
    array([ 1.,  1.,  1.])
    """

    # calculate distances between adjacent points in dual:
    ret_Field = self.dual.delta_dist()
    if self != self.dual:
      # for non-self dual Coord objects: truncate to achieve equal length to self:
      ret_Field.value = ret_Field[1:]

    # update these attributes, as the original copy was based on self.dual:
    ret_Field.grid = self**2
    ret_Field.shape = ret_Field.grid.shape()

    return ret_Field


  def vol(self, gr):

    """
    Generalized volume method related to .d() method: self.d() if self in gr, None otherwise.

    Determines widths of cells along self Coord. grid argument acts as filter: aborts if self not in grid. The grid argument becomes much more critical in some derived classes (e.g. XCoord), where auxhiliary coordinates are picked from Gr and need to be present.

    See .d() method.
 
    Args:
      gr: (Gr) grid to test against.

    Returns:
       Field or None: self.d() if self in gr, None otherwise.

    See also: 
      d method 
      delta_dist method
      der method
      dist method 
    """

    if self not in gr:
      print 'Coord must be in grid argument, returning None.'
      return
    else:
      return self.d()




# -------- End Coord  ----------------


# The following contains two Coord subclasses XCoord and YCoord.
# x is longitude and re-entrant and y is latitude.

class XCoord(Coord):
  """
  Specialized Coord class for representing longitudinal direction in spherical coordinates. A re-entrant geometry is assumed. See Coord.

  The XCoord class contains methods related to distances, with the following dependencies. dist and delta_dist are fundamental (no link, in contrast to Coord). The rest is the same as Coord: der depends on dist. d depends on dual Coord and delta_dist. vol depends on d.
  """


  def roll(self,shift = 0):
    """Yields copy of XCoord (self) with value shifted by integer using numpy.roll().

    Values are modulo 720 degrees from -360. 

    Overrides .copy method of Coord

    Args:
      shift: (int) number of array positions to shift by
 
    Returns:
      Coord: shifted Coord.

    Examples:

    >>> xcoord1 = sg.fieldcls.XCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) )
    >>> xcoord1.roll(1).value
    array([-357.,    1.,    2.])
    """

 
    value = np.roll(self.value,shift = shift)
    if shift > 0:
      value[:shift] -= 360.
    elif shift < 0:
      value[shift:] += 360.
      value -= 360.
   
    return self.copy(name = self.name + '_rolled',value = value)

# belongs to XCoord 
  def coord_shift(self,F,shift, keepgrid = False):
    """
    XCoord method that shifts the coordinates and value of a field by a number of indices. 

    Overides Coord coord_shift method. Here, mask is False, so that the array elements are (1D) rotated: this is the simple way in which this method differs from its Coord counterpart.

    Note that because this method overrides, the trans method also behaves differently because it calls coord_shift, even though the trans method itself is not overriden! 

    The shifted Coord in the grid of the Field argument is replaced with a (different) shifted Coord: disable this behaviour with argument keepgrid = True. Calls roll function. 

    Args:
      F: (Field) field to shift
      shift: (int) magnitude (corresponding to array index) of shift
      keepgrid: (Boolean, default False) grid not updated if True

    Returns:
      The shifted Field, where elements at the very beginning or end of the value array re-appear at the other side of the array (1D-rotation).
    """

    return roll(F,shift = shift,coord = self, mask = False, keepgrid = keepgrid)


  def delta_dist(self,y_coord,fact = R):
    """
    Computes distances between adjacent spherical longitudinal coordinate points of XCoord (self) taking into account latitudinal positions.

    Overrides delta_dist method of Coord class. Does not call dist method (unlike Coord method), but assumes self.value to be an array of longitudinal polar coord values in degrees.

    Args:
      y_coord: (YCoord) latitudinal coordinate positions (usually component in grid context).
      fact: (float) factor by which to multiply result: this should be radius of sphere.

    Returns:
      Field: defined on grid y_coord*self (2D), containing the longitudinal distances.

    See also: 
      d method 
      der method
      dist method 
      vol method      
    """
 
    # crdvals is in degrees longitude
    crdvals = np.roll(self[:],1)
    crdvals[0] -= 360.
    crdvals -= self[:]

    val = np.array([ -fact*np.cos(np.radians(y))*(np.pi/180.)*crdvals for y in y_coord  ])
   
    return Field(name='delta_'+self.name,value = val,grid = y_coord*self, units = self.units) 


  def der(self, F,y_coord):

    """
    XCoord derivative method on Field F. 

    XCoord is cyclical, the first derivative element is not nan. 

    Overrides Coord der method. Differs in that this method takes an extra y_coord argument, which is required for the calculation of cell distances.

    Args:
      F: (Field) Field to take derivative of
      y_coord: (YCoord) latitudinal Coord (usually component in grid context).

    Returns:
      Field on same grid containing the derivative.

    Raises:
      ValueError: when Coord (self) is not in F.grid (grid of Field).

    See also: 
      d method 
      delta_dist method
      dist method 
      vol method      
    """


# Cyclical coords uses different method for ds: it takes the y_coord arg.

    if self in F.grid:
      dF = self.trans(F)
      ds = self.delta_dist(y_coord)

      return dF/ds

    else:

      raise  ValueError("Field argument grid does not contain XCoord %s."%self.name)



  def dist(self, y_coord, fact = R):
    """
    Distance (signed) along XCoord (self) in only one direction (increasing index) from a fixed point.

    The fixed point is usually 0 longitude, but depends on the Coord.

    Overrides dist method of Coord class. Assumes self.value to be an array of longitudinal polar coord values in degrees.


    Args:
      y_coord: (YCoord) latitudinal Coord (usually component in grid context).
      fact: (float) factor by which to multiply result: this should be radius of sphere.


    Returns:
      Field: of shape (len(yt_coord),len(self)), defined on grid y_coord*self, containing the distances.

    See also: 
      d method 
      delta_dist method
      der method
      vol method      
    """

    crdvals = self[:]
    
    return Field(name='distance_'+self.name,value = np.array([ fact*np.cos(np.radians(y))*crdvals for y in y_coord  ])*np.pi/180.,grid = y_coord*self, units = self.units) 




  def d(self,y_coord):
    """
    Calculates width of grid cell in direction of XCoord (self) using the dual of self (e.g. xt_edges). 

    Yields grid cell widths. Can be used to compute volumes.

    Overides the d method of the Coord class.

    Calculates self.dual.delta_dist() where it is defined

   Args:
      y_coord: (YCoord) latitudinal coordinate positions (usually component in grid context).
 
    Returns:
      Field: of shape (len(yt_coord),len(self)), defined on grid y_coord*self, containing the distances between the adjacent coord cell edges.

    See also:
      delta_dist method
      der method
      dist method 
      vol method      

    Examples:

    >>> y_step=30;
    >>> xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    >>> ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )
    >>> xcoord1_edges = sg.fieldcls.XCoord(name = 'testx_edges',direction ='X',value = np.arange(0.,360+45.,90.) -45. , dual = xcoord1  )
    >>> K = xcoord1.d(ycoord1)
    >>> K.shape
    (7, 4)
    """

    # This d method overides the standard d method of the Coord class.
    # It takes the y coordinate object in argument y_coord.
    # When constructing classes derived from the Coord class, use
    # this naming convention for Coord arguments: name argument
    # {x,y,z}_coord to take {x,y,z}coord object. This can then be
    # used to determine volume elements in grid objects (using the inspect module).


    ret_Field = self.dual.delta_dist(y_coord)
    ret_Field.value = ret_Field[:,1:] # truncate field value
    ret_Field.grid = y_coord*self
    ret_Field.shape = ret_Field.grid.shape() # we have truncated the field value, so recalc

    return ret_Field


# --> belongs to XCoord

  def vol(self,gr):
    """
    Generalized volume method related to d method of XCoord: yields self.d(y_coord) if self and a y-coord y_coord in gr, None otherwise.

    Determines widths (1D "volumes") of cells along self Coord. grid argument acts as filter: aborts if self not in grid. The grid argument is more critical in derived classes (e.g. x_coord), where auxhiliary coordinates are picked from Gr and need to be present.

    Overrides vol method of Coord.

    Picks auxiliary coordinate (as when x-widths depend on y) from grid argument gr (y grid is chosen on the same grid as x-coord) depending on what inspect module finds in d method interface
.

    Args:
      gr: (Gr) grid to test against.

    Returns:
      Field or None: self.d() if self in gr, None otherwise.

    Raises:
      RuntimeError if no matching (e.g. y_coord) Coord can be found for d method argument.

    See also:
      d method 
      delta_dist method
      der method
      dist method 
    """

    # Depends on the use of {x,y,z}_coord convention in arguments to d() method of classes derived from Coord class (e.g. XCoord takes y_coord argument).

    if self not in gr:
      print 'Coord must be in grid argument, returning None.'
      return


    coord_types = {'x_coord':XCoord,'y_coord':YCoord,'z_coord':Coord}

    coord_store = {}
# Determine the type of each coord in self
    for r in gr:
      for i in coord_types:
        if isinstance(r,coord_types[i]):
          coord_store[i] = r

   
      # get the Coord-derived objects that need to be passed to each d method of Coord (e.g. xt.d(yt))

    try:
      coords = [coord_store[c] for c in inspect.getargspec(self.d)[0] if c in ['x_coord','y_coord','z_coord']]
    except:
      raise RuntimeError('Error in Coord argument matching for %s. Required Coord likely absent in grid argument %s.'%(self.name,gr))
     

    # Use splat operator * to pass coords list on as argument. Create Field dV.
    # usually this is just the YCoord element of grid.  
    dV = self.d(*coords)
        
    return dV     
  


class YCoord(Coord):
  """
  Specialized Coord class for representing latitudinal direction in spherical coordinates. See Coord.

  The YCoord class contains methods related to distances, with dependencies identical to Coord (as opposed to XCoord). dist is overriden, working through into der, delta_dist, d and vol (even though these affected methods are inherited from Coord).  """

  def dist(self, fact = R):
    """
    Distance (signed) along YCoord (self) in only one direction (increasing index) from a fixed point.

    That fixed point is usually 0 latitude (yielding positive and negative distances), but depends on the Coord.

    Overrides dist method of Coord class. Assumes self.value to be an array of latitudinal polar coord values in degrees. This affects the Coord (parent) methods that depend on it: der, delta_dist, d and vol.


    Args:
      fact: (float) factor by which to multiply result: this should be radius of sphere.


    Returns:
      Field: of shape (len(self),), defined on 1D grid self**2, containing the distances.

    See also: 
      d method 
      delta_dist method
      der method
      dist method 
      vol method      
    """

    return Field(name='distance_'+self.name,value = fact*self[:]*np.pi/180.,grid = (self**2), units = self.units) 



# -------- Ax  definition ----------



class Ax(Directional):
  """
  Axis. Represents direction: e.g. the longitudinal direction, X, or the vertical, Z.

  Coord objects have an attribute that points to an Ax object, representing its direction.

  Attributes: (identical to Directional parent class)
      name: (str) name of Object
      direction: (str) name of direction in which object points
      long_name: (str) longer description (e.g. for display or in Netcdf)
  """
 


  def vcumsum(self,other, upward=True):
    """
    Calls the Coord cumsum method by picking the right Coord from other.grid.

    Fails if Ax not in other.grid Ax objects.

    See Coord.cumsum
    """

    # The product picks the Coord. e.g. X*(xt, yt) is xt 
    return (self*other.grid).vcumsum(other, upward=upward)

# --> belongs to Ax    

  def __xor__(self,other):
    """
    If argument is a Coord, this method tests for equivalence. e.g. xt ~ xu
    If argument is a Field, this method takes the derivative along self Coord axis.
    """
    
  
    # This method works recursively to ensure associativity of the ^ relationship. E.g. a0^a1 is True and a1^a2 is True => a0^a3 is True

    if isinstance(other,Coord) | isinstance(other,Ax):
      if (other in self.equivs) | (self in other.equivs):

        return True
      
      else:
        return False

    elif isinstance(other,Field):
#      return (self*(other.grid)).der(other)
      return self.der(other)

    elif isinstance(other, VField):
      return VField( [self^e for e in other ] )
    else:
      raise Exception('Ax error in %s^%s with Ax %s. Provide Coord, Ax or Field for right member (now %s).' % (self,other,self,other) )
      return

  def der(self,F):
    """
    Derivative method of Ax . Uses entire grid, in case some coords depend on other coords. e.g. x-differentiation requires knowledge of y-position due to nature of polar coords.
    """

    return (F.grid).der(crd = self*F.grid, F = F)      


# --> belongs to Ax 
  def __mul__(self,other):
    """
    Ax multiplication method.

    Yields for arg.:
      Ax/ AxGr: AxGr according to same is_equiv-based rules as Coord 
      Gr: Coord if it is in Gr object and has attribute equal to self, None otherwise.
      Field: calls the Gr.vsum method.
 
    Examples:

    >>> X*(X*Y) 
    X*Y
    >>> X*(latitude*longitude)
    longitude
    """


    if isinstance(other,Ax):   
      # --> multiplication with Ax object: behaves as Coord multiplication.
      if other.is_equiv(self):
        return AxGr((self,))
      else:
        return AxGr((self,other))

    elif isinstance(other,Gr):   
      # --> multiplication with Gr object: yields equivalent Coord in Gr or None.
      if other.eq_in(self):

        return other[other.eq_index(self)]

      else:
#        raise Exception('Axis not in Coord grid.')
        return None 
        
    elif isinstance(other,AxGr):
      if self in other:
        return other
      else:
        return AxGr([self] + list(other) )


    elif isinstance(other,Field) or isinstance(other,VField):
      # --> multiplication with Field object: yields grid method on Field, which is vsum
      # reduce to Coord via multiplication and then to grid method via power of 2:
      return ((self*other.grid)**2).vsum(other)

    else:
      raise Exception('Ax error in %s*%s with Ax %s: provide Coord, Gr or Field object as right multiplicant (now %s). If multiplicant appears to be a Coord of other multiplicant Field, check whether it is stale --> update from Exper Coord stack to be synchronous with the grid of that Field. ' % (self,other,self,other))




# -------- End Ax  definition ----------


class AxGr(tuple, Membered):
  """
  Same as Gr, but containing Ax objects instead of Coord.
  """

  def __repr__(self):    
    rp = '('
    for i in self:
      rp += i.name +','
    
    return rp+')'


  def __and__(self,other):
    """
    Shorthand to element-wise weaksame comparison.
    """

    return self.weaksame(other)


  def weaksame(self,other):
    """
    Member-wise weaksame comparison.
    """

    if len(self) == len(other):
      for i,c in enumerate(self):
        if not(c.weaksame(other[i]) ):    
          return False

      return True

    else:

      return False



  def copy(self):
    """
    Create new AxGr object via element-wise copy calls.
    """

    return AxGr( [ e.copy() for e in self  ] )

  def __div__(self,other):
    """
    Division of AxGr. Same rules as Gr: E.g. xt*yt*zt/yt yields xt*zt

    Examples:

    >>> (X*Y)/Z
    (X,Y,)
    >>> (X*Y)/X
    (Y,)
    """

    if isinstance(other,Ax):
      other = AxGr((other,))
    elif isinstance(other,Coord):
      other = Gr((other,))

    result = list(self)
    for it in self:
      if other.eq_in(it):
        result.remove(it)
    return AxGr(result)




  def __mul__(self,other):
    """
    Multiplication of Ax grids.


    Multiplication can take other arguments than just Ax grids. If a Field is provided as right multiplicant, the Field is summed over the left multiplicant grid, weighted with grid cell widths (the equivalence of integration over the grid space). If the right multiplicant is a Coord object, it is converted to a single-element grid (gr) object before multiplication. 

    if right multiplicant is Gr object, operation picks elements from right multiplicant that are equivalent with Ax objects in AxGr object left multiplicant and yields a product in the order of the left multiplicant.

    Raises:
      TypeError: when inapparopriate type is used.
  

    Examples:

    >>> (X*Z)*(zt*yt*xu)   # note the order of the output elements
    (xu, zt)
    """


    if isinstance(other,Coord):

      if self.eq_in(other):
        return other
      else:
        return 

    elif isinstance(other,Ax):
      # multiplication between Gr Ax and Ax objects
      if other in self:
        return self
      else:
        return AxGr(list(self) + [other])


      return other*self

    elif isinstance(other,AxGr):

      return reduce(lambda x,y: x*y , list(self) + list(other))

    elif isinstance(other,Gr):
      # if right multiplicant is Gr object, operation picks elements from right multiplicant that are equivalent with Ax objects in AxGr object left multiplicant and yields a product in the order of the left multiplicant.
      L = []

      for it in self:
        if other.eq_in(it):
          L.append(other[other.eq_index(it)])

      return Gr(L)    

    elif isinstance(other,Field) or isinstance(other,VField):
     
      return (self*other.grid).vsum(other)

    else:
      raise TypeError('gr type error %s*%s with Gr %s (grid): provide Field, Gr or Coord object or np array as right multiplicant (now %s).' % (self,other,self,other) )


  def eq_in(self, crd):
    """ Determines whether Coord crd is equivalent to any of the constituent Ax objects of this AxGr 

    Gr and AxGr have an eq_in method that acts on a Coord or Ax.

    Args:
      crd: (Coord) object to be checked.

    Returns:
      True when crd is equivalent to one of the Ax objects, False otherwise.

    See also:
    eq_in method of Gr     
    """

    for i in self:
      if crd.is_equiv(i): return True
    return False




#------------------------- end Ax and AxGr  -------------------------------


# -------------- grid  --------------------------


class Gr(tuple):

  """
  Represents Coord grids. Consists of a tuple of Coord objects, with additional methods. Gr objects g1 and g2 are considered weaksame, g1.weaksame(g2) yields True, when the individual Coord elements are weaksame.
  """



  def __eq__(self,other):

    if len(self) == len(other):

      return reduce(lambda x,y: x and y, [ np.array_equal(e[:], other[i][:]) for i,e in enumerate(self)  ] )

    else:
      return False

  def same(self, other):
    """    
    Coord-wise comparison of elements using Coord same method. 
    """
    if len(self) == len(other):

      return reduce(lambda x,y: x and y, [ e.same( other[i] ) for i,e in enumerate(self)  ] )

    else:
      return False

  def __call__(self,other, method = 'linear'):

    return self.regrid(other = other, method = method)

  def regrid(self,other, method = 'linear'):
    """
    Takes another Gr object argument and yields a function F from ndarray to ndarray.

    The regrid method of a Gr object takes another Gr object, other, and yields a function F. This function F takes an array A and re-arranges the order of   the indices to match the input Gr object (other). If the length of the input object exceeds that of the calling object, F(A) also expands the array along the additional axes by creating copies of it along those axes (using the expand method). Note that the coords of the calling Gr object need to be a subset of the argument Gr object.


    Args: 
      Other Gr object (grid). 
    

    Returns: 
      A transformation on fields going from self grid to other grid.

    E.g. xt*yt(yt*xt) yields a tranpose operation on an array
    xt*yt(xu*yu) yields an interpolation acting on fields.

    yt*xt(zt*yt*xt) yields a functions transforming a 2D array corresponding to the values of a Field defined on yt*xt to a 3D array constant in the zt direction.
    
    If self is longer than other, calling will lead to a reduction. E.g.

    R=(zt*yt*xt)((yt**2))(A) where A.shape = (len(zt),len(yt),len(xt))
    Leads to a list of length len(yt) containing arrays of dimension len(zt) by len(xt).
    Then for S=array(R) we get S.shape is (len(yt), len(zt), len(xt))
    
    For R=(zt*yt*xt)((xt*yt))(A) we get a list of lists and S.shape is (len(xt), len(yt), len(zt))
    Note that yt,xt appear in different order in self and other.
    """

# might expand other to other*self in future code.
    
    if len(self) == len(other):    
# *** CASE 1 ************************

      # in this case both grids span the same space

      # check if a permutation of Coord objects exists, i.e. whether the elements of either can be rearranged to yield the other: 
      pm = self.perm(other,verbose = False)

      if pm:
        # CASE 1a ***

        # if so, a function is returned that attempts to transpose any np array according to the permutation required to go from self to other.
        # If A is a np array defined consistent with self ( A.shape is self.shape() ), then self(other)(A) is a np array consistent with other
        return lambda A: np.transpose(A,pm)
      else:
        # if no such direct permutation exists, check for a weaker conditions:
        pm = self.eq_perm(other)
        # if there is a permutation of the coords up to equivalence (case 1b below), pm is that permutation and after permuting array A, the result needs to be interpolated from the permuted self (namely self.rearrange(pm)) to other.
        # Here, "up to equivalence" means "equivalent Coord objects being considered identical in considering whether a permutation exists".
        #If A is a np array defined consistent with self, then self(other)(A) is a np array consistent with other, is interpolated onto the other.

        if pm:
          # CASE 1b ***
          return lambda A: (self.shuffle(pm))._smart_interp(np.transpose(A,pm),other, method = method)
        else:
          # CASE 1c ***
          # No luck.
          print "grid %s not equivalent to %s."%(self,other)
          return

    elif len(self) < len(other):  

# *** CASE 2 ************************

      # A grid is called on a higher dimensional grid.

      # inflate self grid by left multiplying with non-self elements
      # don't do deepcopy, it copies the individual Coord elements too!
      # instead, use the identity Coord:

      # re-arrange Coord terms in accordance with expand method (the non-self elements of the other Gr are appended on the LEFT, and in the order of other).
      # e.g. R=(yt*xt)(xt*yt*zw) yields a function yielding an array defined on the grid zw*yt*xt
      # the order of the Coord elements in self_expanded is arranged so as to perform the expansion more easily (namely, adding axes at the beginning).

      # if A is an ndarray consistent with self, then self.expand(A,other) is an ndarray consistent with the Gr self_expanded created here:
      self_expanded = (other/self)*self

      # we now have a grid of equal length to other.

      # the expanded left argument is not always in the same order as other (or even fully comprising of identical elements)
      pm = self_expanded.perm(other, verbose = False)
     
      if pm:
        # case 2a
        
        # In this case other contains only Coord elements from self_expanded.
        # return function that takes ndarray A consistent with self and returns A expanded to other (yielding array of same dimension as other) and then transposed to be consistent with other.
        # this should yield the same result as using other as argument for expand
        return lambda A: np.transpose(self.expand(A,self_expanded),pm)

      else:
        pm = self_expanded.eq_perm(other, verbose = False)
        if pm:
          # case 2b

          # line up the equivalent Coord elements in the same order for interpolation.
          return lambda A: (self_expanded.shuffle(pm))._smart_interp(np.transpose(self.expand(A,self_expanded),pm),other, method = method)
        else:
          # case 2c
          print "grids not equivalent"
          return

      return 

    else:
#**** CASE 3 ************************

      # This is the case where len(other) < len(self) => reduce method. This yields a function that slices along the Gr provided in the argument, and a permutation among those axes if they appear in different order in self and other.

# To illustrate this functionality:

#If V is an ndarray consistent with zt*yt*xt
# and we do R = np.array((zt*yt*xt)(yt*xt)(V));R = R.reshape((yt*xt*zt).shape())
# Then We get V back, but transposed onto yt*xt*zt. This is because the index grid L = yt*xt is put first by (zt*yt*xt)(yt*xt)


      # create target_grid of same dimension as self, and with other Coord elements first
      target_grid = other*(self/other)
      
      pm = self.perm(target_grid, verbose = False)      
     
      # Using reduce method. Note that reduce has the arguments the other way around. reduce is called as method of other!
      if pm:
        # case 3a
        return lambda A: other.to_slices(np.transpose(A,pm),target_grid)
      else:
        pm = self.eq_perm(target_grid, verbose = False) 
        if pm:
          # case 3b
          return lambda A: other.to_slices((self.shuffle(pm))._smart_interp(np.transpose(A,pm),target_grid, method = method),target_grid)

        else:
          # case 3c
          print 'Nope'
          return
     
    return

  def function(self,func):
    """
    Returns a Field containing the values of function argument func on the grid points defined in this grid. 

    Args:
      func: (function) function defined on domain of same dimension as grid (self).

    Returns:
      Field: Evaulated function. The Field name is the name of the function.
    """

    vfunc = np.vectorize(func)
    value = vfunc(*self.inflate())
 
    return Field(name = func.func_name, value = value, grid = self)

  def array_equal(self,other):
    """
    Grid component-wise test whether the Coord objects contain the same grid point location values. Input another grid.

    Args:
      other: (Gr) grid to compare with.

    Returns:
      List. List element corresponds to Coord in Gr and is True if Coord.array_equal True, False otherwise.

    Raises:
      ValueError if Gr objects not defined along same axes. E.g. (X,Y ) vs (Y,Z )
    """

    if not self.axis() == other.axis():
      raise ValueError('Provide grids defined along same axes.')

    return [e.array_equal(other[i]) for i,e in enumerate(self)   ]
        

  def axis(self):
    """
    Returns an AxGr object containing the axis properties of the Coord elements of this grid.
    """
    return reduce(lambda x,y:x*y,[e.axis for e in self])


  def reverse(self):
    """
    Reverse the order of the grid elements.

    Examples:

    >>> coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) )
    >>> coord2 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]) )
    >>> (coord1*coord2).reverse()
    (test2, test1)
    """

    return Gr([ self[len(self) -i -1 ] for i in range(len(self))  ])


  def nbytes(self):
    return reduce(lambda x,y:x.nbytes+y.nbytes, self)

  def __and__(self,other):
    """
    Shorthand to weaksame
    """

    return self.weaksame(other)


  def weaksame(self,other):
    """
    Application of component-wise weaksame method of Gr members.

    Args:
      other: (Gr) grid to compare with

    Returns:
      True if self and other of equal length and weaksame yields true for each self-other member pair.
    """


    if len(self) == len(other):
      L = [e.weaksame(other[i]) for i,e in enumerate(self)]
      return reduce(lambda x,y: x and y, L)
  
    else:
      return False 



  def copy(self):
    """
    Creates object with same values. A = B.copy() yields A&B = True, see __and__ method.
    """

    return Gr([e.copy() for e in self])


# ------------------------------------------
# Lower level methods:

  def to_slices(self,A,other):        
    """
    yields a list of slices along the coords defined in self. e.g.
    zt(zt*yt*xt) = [A[0,:,:],A[1,:,:],...] where A.shape is (zt*yt*xt).shape()

    Expects self coords to be subset of other, and appearing in same order in both.
    other must appear in the left side of self (i.e. self is self*(other/self)  ).
    For instance, zt*yt(zt*yt*xt) is valid,  yt*xt(zt*yt*xt) and zt(yt*xt) are not.
    The indexing in the output list (as list of lists) is of opposite order to the Coord elements in self.

    No checks are done on consistency between A or other or self and other.

    The opposite of expand. Used by call method of fields on Gr objects of lower dimension that the Field.

    Inputs: 
    A		ndarray of shape other.shape()
    other	another larger Gr (grid) object containing self

    Outputs:
    A list of nparrays being slices of input A along the self Gr.

    Note that argument is longer than self. This is opposite to __call__ method, where a longer self leads to a reduction.
    """

    # This method works with recursion. If len(self)>1, a list is built using this method on the smaller elements and indexing by the first dimension.
    # if B = self.to_slices(A,other) and A is an array, then array(B) has the same shape and values as A.
    # calling say (zt*yt).to_slices(A,zt*yt*xt) yields a list of lists. Each of those lists then contains a 1D array.



# the following code is made slightly difficult due to recursion:
#    if (force == False) and ( other != self*(other/self) ):
#      warnings.warn('Calling Gr %s on Gr %s. other must equal self*(other/self).'%(str(self),str(other) ) )

 
    if len(self) == 1:
      # single Coord Gr called on Gr. Endpoint of recursion.
    
      result = [] 
      coord = self[0]
      if isinstance(A,list):
        for i in range(len(coord)):
        # Numpy note: the colon here in A can represent several dimensions.
          result.append(A[i])   
       
      else:
        for i in range(len(coord)):
  
        # Numpy note: the colon here in A can represent several dimensions.
          result.append(A[i,:])   

      return result
    else:
      # multiple coord Gr called on Gr. Recursion until single coord Gr called on Gr.

      result = []
      subself = Gr(list(self))

      while len(subself) > 1:
        coord = subself[0]

        subself = subself/(coord**2)
        if isinstance(A,list):
          for i in range(len(coord)):
            result += (subself.to_slices(A[i],other))

        else:
          for i in range(len(coord)):
            result += (subself.to_slices(A[i,:],other))

    return result
        


  def expand(self,A,other):
    """
    Adds dimensions specified in Gr other at the beginning of array A

    input: ndarray of shape consistent with self, Gr other
    output: an ndarray of shape (other/self)*self containing identical copies of A along other/self

    Example.

    SAT = P['DPO']['A_sat']
    SAT.shape is (100,100)
    W=SAT.grid.expand(SAT[:],depth**2)
    W.shape is (19,100,100)
    W contains 19 identical copies (slices) of SAT[:] 

    Note that the other grid is appended on the left side.

    Example 2.

    (zt*yt*xt).shape() is (46, 110, 200)

    A = np.ones((xt**2).shape())

    K=(xt**2).expand(A,zt*yu*xt  )
    K.shape is (46, 110, 200)

    K=(xt**2).expand(A,zt*xt*yt  )
    K.shape is (46, 110, 200)

    Warning: method requires Gr argument, do not use coord argument. Instead, for a single coord (e.g.) xt, use xt**2.
    """

    # only use those coord elements that are in the other gr but not the self gr.
    new_coords = list(other/self)
    new_coords.reverse()

    for coord in new_coords:
      # initialize L for each coord with the A from argument, or just built below
      L = np.array([A],)
      for k in range(len(coord)-1):
        # grow the dimensions with identical copies of A
        L = np.append(L,[A,],0)
      A = L 

    return L




  def inflate(self, type = 'array', force = False):
    """
 
    Input:
    type = output type. -'array' in arguments will return a list of arrays. -'Field' in arguments will return a list of fields.

    Output: 
    A list of arrays or fields of the dimension of the grid being called.
    Each element in the list corresponds to a Coord object in the called grid, where the array equals the content of the Coord along the array index corresponding to that Coord, and is constant otherwise.


    For example, a grid defined by (yt,xt) (equal to yt*xt) yields [YT,XT] where YT = yt(yt*xt) and XT = XT(yt*xt). We refer to XT as the inflated version of xt. Here, the Coord object has been called on the grid object: this yields an array defined on the argument grid and constant in all Coord axes other than the calling Coord. The array equals the value of the calling Coord object along that axis.


    Cached for performance. Refresh with force = True.
    """

    if  not(hasattr(self,'inflated')) or (not self.inflated) or (force == True):
      # compute values and store as arrays.
        
        # This yields a list of arrays, corresponding to the inflated Coord objects.
        self.inflated = [e(self).value for e in self]
        
    if type == 'array':
      return self.inflated
    elif type =='Field':
      return [Field(name = 'inflated_'+self[i].name,value = e, grid = self ) for i, e in enumerate(self.inflated ) ]


  def _smart_interp(self,A,other, method = 'linear'):

    """
    Inputs: an array A of the shape corresponding to self.
            a destination grid.

    Outputs: an array containing A interpolated from the self grid to the destination grid.


    Smart interpolation of array A, using griddata interpolation only along Coord axes that are not equal (but must be equivalent).
    !!!Arguments must be in the right order: order(left) = order(right)!!!

  

    """
# belongs to grid object.

# it is the left element. Coord elements of self and other (grids) may be up to equivalence, but need to be in same order. The common (equal) elements will be stored in L


    L = []
    for i, it in enumerate(self):
      if it == other[i]:
        L.append(it)
      else:
        if not(it^other[i]):
          print "order/ equivalence wrong, aborting."
          return

    if L:
      # In this case, the source and destination grid have Coord elements in common. This means we need to interpolate only along the axes they do not have in common.

      # check first whether source and destination grids are equal, in which case we can simply return A.
      if len(L) == len(self):
        return A

      # In this case, source and dest grids are not equal, and contain both equal and equivalent-only Coord elements:
      L = Gr(L);
      
      result = []

    # Take slices into a list
      B=self(L)(A)      
  
#    B is a (often long) list containing the slices to be interpolated.
# array(B) will yield shape (len(coord1)*len(coord2)*...  , shape(array)) for coordi in L and array the slices to be interpolated
# This should be reshaped to list(L.shape()) + shape(array)
# where shape(array) = (self/L).shape()

      for i,b in enumerate(B):
      # perform interpolation on array b from self/L to other/L on each slice.
        srcgrid = self/L
        destgrid = other/L

        B[i] = srcgrid._interp(b,destgrid)
     
# B has now been interpolated.
      B = np.array(B)

# some commented out diagnostic prints
#      print L
#      print self/L

#      print list(L.shape())
#      print list((self/L).shape())

#      print B.shape

      B = B.reshape(list(L.shape()) + list((other/L).shape()) )

      pm = (L*(self/L)).perm(self, verbose = False)
      B = np.transpose(B,pm)

      return B
    else:
      return self._interp(A,other, method = method)

# methods belong to Gr 

  def _interp(self,A,other, method = 'linear'):


# it is assumed that self^other and that the shape of array A corresponds to the lenghts of the Coord elements of self (and therefore other).
    if not(self^other):
      print 'Arguments not equivalent. Use equivalent grids.'
      return
 
    if len(self) == 1:
    
      L = self[0][:]
      R = other[0][:]

      if L[2]<L[1]:
        # in case it's a negative scale, as in depth
        L = -L
        R = -R
       
      IA = griddata(L,A,R, method = method)
      
    
      return IA
     
    else:
  
      L = np.array([ e.ravel() for e in self.inflate() ]).transpose()
      R = np.array([ e.ravel() for e in other.inflate() ]).transpose()  

      IA = griddata(L,A.ravel(),R, method = method)
      IA=IA.reshape(other.shape())
    
      return IA


  def __div__(self,other):
    """
    Division of grids. E.g. xt*yt*zt/yt = xt*zt
    """

    if isinstance(other,Coord):
      other = Gr((other,))
    elif isinstance(other,Ax):
      other = AxGr((other,))


    result = list(self)
    for it in self:
      if other.eq_in(it):
        result.remove(it)
    return Gr(result)


  def __mul__(self,other):
    """
    Multiplication of grids.
    At the moment, xu*zt*xt*yt = (xu,zt,yt,) whereas xu*(zt*xt*yt) = (zt,xu,yt,)

    Multiplication can take other arguments than just grids. If a Field is provided as right multiplicant, the Field is summed over the left multiplicant grid, weighted with grid cell widths (the equivalence of integration over the grid space). If the right multiplicant is a Coord object, it is converted to a single-element grid (Gr) object before multiplication. 
    """

    if type(other) == np.ndarray:
      # multiplication with an array yields a Field if the sizes match.
      if self.shape() == other.shape:
        return Field(name = '', value = other, grid = self)
      else: 
        raise Exception('Gr shape error %s*%s with Gr %s: provide correct shape np array.' % (self,other,self) )
        return   
 
    elif isinstance(other,Field):

      return self.vsum(other)


    elif isinstance(other,Coord):
      if other.name == 'scalar':
        return self
      else:
        if self.eq_in(other):
          return self
        else:
          return Gr(list(self) + [other] )

    elif isinstance(other,VField):
      # multiplication of grid object with vector Field.
      # this commutes:
      return other*self

    elif isinstance(other,Ax):
      # multiplication between Gr and Ax objects is commutative
      return other*self

    elif isinstance(other,Gr):

      return reduce(lambda x,y: x*y , list(self) + list(other))

    elif isinstance(other,AxGr):
      # --> multiplication between Gr and AxGr objects DOES NOT commute: in agreement with general rules, result retains Coord element order of left multiplicant.
      
      L = []

      for it in self:
        if other.eq_in(it):
          L.append(it)

      return Gr(L)    


    else:
      raise Exception('Gr type error %s*%s with Gr %s (grid): provide Field, Gr or Coord object or np array as right multiplicant.' % (self,other,self) )


# belongs to grid 
  def strict_equiv(self, other):
    """
    Tests whether two Gr objects have equivalent Coord elements at each position.
    This is a stricter test than grid (Gr object) equivalence testing via gr1^gr2, which only tests whether both grids describe the same space (elements equivalent up to a permutation).

    """
    if len(self) == len(other):
      
      RO = True      
      for i,it in enumerate(self):
        RO *= (it^other[i])
      return bool(RO)
    else:
      return False


  def __xor__(self, other):
    """
    Checks equivalence between grids, where grids are equivalent if they define the same physical subspace, based on the equivalence definition for Coord classes. In other words, checks whether the individual Coord elements of the two grid (Gr object) arguments are equivalent up to a permutation. A stricter version of this test is strict_equiv, which allows no permutation.
    """

    if len(self) == len(other):
      
      if self.eq_perm(other):
        return True
      return False

    else:
      return False


  def eq_in(self, crd):
    """ Determines whether Coord crd is equivalent to any of the constituent Coord objects of this Gr

    Gr and AxGr have an eq_in method that acts on a Coord or Ax.

    Args:
      crd: (Coord or Ax) object to be checked.

    Returns:
      True when crd is equivalent to one of the Coord objects, False otherwise.

    See also:
    eq_in method of GrAx     
    """


    for i in self:
      if crd^i: return True
    return False

  def eq_index(self,crd):
    for i,v in enumerate(self):
      if crd^v: return i
    return -1

  def shuffle(self,permutation):
    """
    Gr method that rearranges the order of the elements of this grid Gr object via permutation arrgument. E.g.

    g1 = latitude*depth
    g1.shuffle( (1,0) ) is (depth, latitude)

    See also perm method of Gr.
    
    """
    return Gr((self[i] for i in permutation))
    
    


  def perm(self, other,verbose = False):     
    """
    yields permutation of axes going from grid self to grid other.
    E.g. for grids gr1 and gr2, g2 = g1.shuffle( g1.perm(g2) )

    Returns None if no permutation exists.

    See also shuffle.

    """  
    return find_perm(self,other,verbose = verbose)

  def eq_perm(self, other, verbose = True):      
    """
    yields permutation of axes going from grid self to grid other, where equivalent axes are treated as identical. 

    See also perm.
    """  

    if len(self) == len(other):
      perm = []
      for r in other:
        if self.eq_in(r):
          perm.append(self.eq_index(r))
        else:
          warnings.warn( 'Warning from eq_perm (often benign): inputs not permutable, returning None.')
          return

    else:
      if verbose:
        print "Message from eq_perm: inputs must be of equal length."
      return 
    return tuple(perm)

  def shape(self):
    sh = [];

    for c in self:
      if (not isinstance(c[:],str)) and (not isinstance(c[:],unicode)):
        sh.append(len(c))
      else:
        sh.append(-1)

    return tuple(sh)
 

  def dual(self):
      """
      Method of Gr object that returns a grid made up of the dual Coord elements of this Gr object.
      """
      gr_dual = Coord('scalar')
      for e in self:
        gr_dual *= e.dual
      return gr_dual

  def ones(self):

    return Field(name = 'ones', value = np.ones(self.shape() ) ,grid = self )


  def vsum(self, F):
    """
    Method of Gr object. Sum weighted with Coord grid cell widths (integration) over self grid. 

Takes Field argument and returns a Field with grid made up of remaining Coord objects or a float. E.g. if F.grid == ('zt','yt','xt'), (xt*yt).vsum(F) yields a Field defined on grid ('zt',).

    
    Note that when Coord elements with direction attribute 'X' and 'Y' both appear in the Gr object, vsum will check whether the 'X' Coord appears after the 'Y' Coord. If so, they will be interchanged when performing the calculation as otherwise no y-coord is available when the x grid cell width is required. This is a small detail.

    """


    # If X and Y directions both occur in the averaging grid, need to make sure X appears before Y because X-direction Coord grid cell width depends on Y-direction Coord.

    dirs = [e.direction for e in self]

    if ('X' in dirs) and ('Y' in dirs):
      i_X = dirs.index('X') 
      i_Y = dirs.index('Y')

      if i_X > i_Y:
        new_gr = list(self)      
        tmp = self[i_X]
        new_gr[i_X] = self[i_Y]
        new_gr[i_Y] = tmp
        del tmp
      
        return reduce(lambda x,y: y.vsum(x), [new_gr[0].vsum(F)] + new_gr[1:] )


    # Apply Coord vsum method of Coord objects in self to Field argument F, from left to right:
    return reduce(lambda x,y: y.vsum(x), [self[0].vsum(F)] + list(self[1:]) )

  def mean(self,F):
    """
    Method of Gr object.
    Determines mean of Field argument F weighted with grid cell width.
    """
    return self.vsum(F)/self.vsum(F.ones())

# --> belongs to Gr 


  def der(self,crd,F):
    """
    Method of grid object. Often the wider context of the grid needs to be known to take the derivative along a Coord, hence a Gr method.

    Input:
    crd		Coord object along which to differentiate (e.g. latitude)
    F		Field object to differentiate (e.g. temperature)


    """
    coord_types = {'x_coord':XCoord,'y_coord':YCoord,'z_coord':Coord}
 
    if crd in self:
      C = self.find_args_coord(method_name = 'der', coord_types = coord_types)    
      i = self.index(crd)  

      return crd.der(F,*C[i])

    else:

      raise Exception('Error in Gr derivative method der. %s must be in grid %s') % (crd, self)

    
  def vol(self):
    """
    Gr method that determines volumes (areas/ lengths) of grid elements, returns Field.
    """
    # Depends on the use of {x,y,z}_coord convention in arguments to d() method of classes derived from Coord  (e.g. XCoord takes y_coord argument).

    coord_types = {'x_coord':XCoord,'y_coord':YCoord,'z_coord':Coord}

    
    C = self.find_args_coord(coord_types)

   
    # Use splat operator * to pass coords list on as argument
    # cycle through coords, the list of Coord elements required as arguments for each Coord, 

    return reduce(lambda x,y : x*y, [r.d(*C[i]) for i,r in enumerate(self)]  )     
  
  def find_args_coord(self,coord_types, method_name = 'd'):

    coord_store = {}
# Determine the type of each Coord in self
    for r in self:
      for i in coord_types:
        if isinstance(r,coord_types[i]):
          coord_store[i] = r

    L = []
    for r in self:
      # get the Coord-derived objects that need to be passed to each d method of Coord (e.g. xt.d(yt))
      exec 'method = r.' + method_name
      coords = [coord_store[c] for c in inspect.getargspec(method)[0] if c in coord_types.keys()]

      L.append(coords)       

    return L

# ---------------- end Gr  ----------------





# -------------- Field  --------------------------


  
class Field(Named):
  """
  Field  to represent a scalar valued function, also dataset, defined on a grid.

  The call method allows fields to act as function defined on grid objects.
  For a 2D scalar corresponding to Field T, say defined on grid yt*xt, T(dy*dx) yields a 2D array of the scalar values.

  If Field T is naturally defined on grid yt*xt, then T(zt*yt*xt) yields a 3D array b such that b[k,:,:] = T(yt*xt) for all possible k.


   If g is a Gr (grid) or Coord object, left or right multiplication of a Field  object F with Gr Coord results in the grid-cell width weighted summing of the Field over the coords in the multiplicant g (integration, via g.vsum method), resulting in a smaller dimension Field.

  If g is a Coord object, g^F yields the derivative of F along g (via g.der method). g|F yields the grid cell width-weight cumulative sum of F over g (primitive, via g.vcumsum).

  two fields F1, F2 are considered equal, F1&F2 yields True, when their name, value (an numpy array) and Gr (grid) attribute are equal, unless they contain nan values.

  NOTE: multiplication works a bit different from addition at the moment. Addition will go ahead even when coords in the grids are differently named (or have other non-value attributes differ) as long as the value (the coord points) are the same: then the (left and right) coords are considered equal. Multiplication treats them as different coords in this case.

  """

  global ID

  def __repr__(self):
    return self.name

  def __init__(self,name,value,grid,units = '?',direction = None, strict_v = strict_vector,long_name='?',metadata={} , squeezed_dims =Gr( ()  )):
    """
    Initialise a Field. 
    Inputs: 

    name:	the name of the Field (e.g. temperature). Displayed in console
    value:	the numpy array containing the Field data
    grid:	the grid Gr object associated with the data
    units:	data units (if known)
    direction:	scalar or, if vector Field component, axis direction (e.g. X)
    strict_v:	if True (default), addition of directional fields leads to vector fields.
    long_name:	Description of Field, corresponds to long_name Netcdf metadata.

    These inputs become attributes of the created Field object.

    """


    if not(direction):
      direction = ID()

    if isinstance(value,np.ndarray):
      if isinstance(grid,Gr):
        shape = grid.shape()

        # Check for shape consistency between grid and value to avoid creation of inconsistent Field objects:
        if shape == value.shape:
          self.name = name
          self.value = value
          self.grid = grid
          self.shape = shape
          self.units = units 
          self.direction = direction 
          self.strict_v = strict_v
          self.long_name = long_name
          self.metadata = metadata
          self.nbytes = value.nbytes

          self.squeezed_dims = squeezed_dims

        else:

          raise Exception('Error in Field creation %s using grid %s: value array argument must have same shape as grid argument! Gr shape %s while Field shape %s ' %(name,grid,str(shape),str(value.shape) ) )
          return
      else:
        raise Exception('Error in Field creation %s: argument grid %s must be a Gr object!' % (name, grid))
        return
    else:
      raise Exception('Error in Field creation %s: argument value must be an ndarray!' % name )
      return

  def __and__(self,other):
    """
    self&other
    Tests whether fields contain equal values. At the moment, if the value contains nan, this function will return false.
    """
    if (self.name == other.name) and np.array_equal(self.value,other.value) and self.grid&other.grid:
      return True
    else:
      return False


  def cdf_insert(self,file_handle, insert_dual = True, force_squeeze = False, miss_default = 9.96921e+36):
    """
    Netcdf insert method of Field .

    Writes Field to already opened file referred to with file_handle argument, along with its Coord objects.

    """


    # handle the squeeze dimensions
    if not force_squeeze and len(self.squeezed_dims) > 0:
      return unsqueeze(self).cdf_insert(file_handle = file_handle, insert_dual = insert_dual)    

    # insert the coords in own grid
    for crd in self.grid:
      if not crd.name in file_handle.variables:
        crd.cdf_insert(file_handle)
    
        if insert_dual and (crd.dual != crd) and (not crd.dual.name in file_handle.variables):
          crd.dual.cdf_insert(file_handle)         

    # This could bloat memory. Redo in a new way.
    value = copy.deepcopy(self[:])

    miss_val = miss_default
    if 'FillValue' in self.metadata:
      miss_val = self.metadata['FillValue']
    elif 'missing_value' in self.metadata:
      miss_val = self.metadata['missing_value']

    try:
      value[  np.isnan( value  ) == True   ] = miss_val
      
    except:
      try:
        value[  np.isnan( value  ) == True   ] = miss_default

      except:
        warnings.warn('Could not set missing value for Field %s.'%self.name)




    # Create the actual variable corresponding to Field.value
    var_cdf = file_handle.createVariable(self.name, value.dtype.char, tuple( [crd.name for crd in self.grid] )   )
    var_cdf[:] = value


    for k in self.metadata:
      setattr(var_cdf,k, self.metadata[k]) 

    return file_handle


  def write(self, path = None, name = None , history = 'Created from Spacegrids ' , insert_dual = True, force_squeeze = False ):

    """
    Write method of Field .

    Creates Netcdf file and writes Field to it, along with its Coord objects.

    Fields are unsqueezed before saving, along Coord objects of single length to be saved as well (override with force_squeeze = True).

    If path and name are not specified, the file will be located in the working directory.
    If only name is specified, the file will be in the wd under that name
    If path is specified, the wd is replaced by the path in the above 2 scenarios.

    insert_dual determines whether the edges of a the Coord objects are saved as well (the default).
    

    """

    if name is None:
      name = self.name

    if not name.split('.')[-1] in ['nc','cdf']:
      name = name +'.nc'
    if not path is None:
      name = os.path.join( path , name ) 
   
    

    print 'Writing Field to file %s'%name
    try:
      file_handle = netcdf_file(name , 'w')
    except IOError:
      print 'Cannot write ', name
    else:

      file_handle = self.cdf_insert(file_handle, insert_dual = insert_dual , force_squeeze = force_squeeze)

      file_handle.history = history + '%s'%str(datetime.datetime.now())
 
#    var_cdf.units = self.units

      file_handle.close()

  def cat(self,other,ax = None, name_suffix = '_cat'):
    """
    Concatenate with another Field along axis ax. If ax is None, concatenation takes place along the first encountered common axis with non-equal values.
    Grids must be orient along same axes and in same axis order.

    """

    # if no Ax object is given, an Ax is chosen where the grid Coord elements are not array equal.

#    if len(self.grid) != len(other.grid):
#      raise Exception('Error: provide grids of equal dimension.')

    if isinstance(ax,Coord):
      ax = ax.axis

    self_axis = self.grid.axis()

#    if not reduce(lambda x,y:x and y, [e^other[i] for i,e in enumerate(other)]):
    if self_axis != other.grid.axis():
      raise Exception('Error: provide fields defined on the same grid directions.')


    if ax is None:
      
      i_ax = (self.grid.array_equal(other.grid)).index(False)
      ax = self_axis[i_ax]
    
    cat_coord_self = ax*self.grid
    
    if cat_coord_self is None:
      # in this case concat is done along an axis not in the self grid
 
      raise Exception('Axis not in grid.')


    if (self.grid/ax).shape() != (other.grid/ax).shape():

      raise Exception('Field concat error %s and %s. Provide pieces of right dimensions. (now %s and %s)'%(self.name,other.name, str((self.grid/ax).shape())  , str( (other.grid/ax).shape())   ) )

      # obtain the index of the axis in the grid along which to concatenate.
      # why do we need eq_index here instead of index? because it can be an Ax object.
    ax_index = self.grid.eq_index(ax)

    # combine the two halves as dictionaries of slices of what is to be the new Coord first
 
    # pick the Coord specified by the ax argument by multiplying the grids: 
    left_coord = (ax*self.grid)
    right_coord = (ax*other.grid)
    # e here is a point in the relevant Coord: 
    Dleft = {e:self[ax,i] for i, e in enumerate( left_coord[:] ) }
    Dright = {e:other[ax,i] for i, e in  enumerate( right_coord[:] ) }

    # if one or both coords have no strings attribute set, don't give the new Coord a string attribute either.
    if (left_coord.strings is not None) and (right_coord.strings is not None):
      stringsleft = {e:left_coord.strings[i] for i, e in enumerate( left_coord[:] ) }    
      stringsright = {e:right_coord.strings[i] for i, e in enumerate( right_coord[:] ) }    
      stringscomb = dict(stringsleft.items() + stringsright.items() )
    else:

      stringscomb = None

    Dcomb = dict(Dleft.items() + Dright.items() )

      # each unravelled piece needs to have the right shape for np concatenation.
    piece_shape = list(self.shape)
    piece_shape[ax_index] = 1

      # use combined keys to construct ordered values of new concatenated Coord object.
    cat_coord_value = np.array(Dcomb.keys())
    cat_coord_value.sort()

      # create the new concatenated Coord object using the combined ordered sequence of values.
    if stringscomb is not None:
      new_strings = [stringscomb[k] for k in cat_coord_value]

      new_coord = cat_coord_self.copy(name = cat_coord_self.name +name_suffix,   value = cat_coord_value, strings = new_strings)

    else:
      new_coord = cat_coord_self.copy(name = cat_coord_self.name +name_suffix,   value = cat_coord_value, strings = None)



    new_coord|cat_coord_self
      # construct combined Field values. Reshape is needed for np.concatenate function.
    values = [Dcomb[k][:].reshape(piece_shape) for k in cat_coord_value]
   
    new_value = np.concatenate(values,axis=ax_index)

      # construct the grid of the combined object by replacing the old partial Coord with the new combined Coord in the self grid. Recall that replacement is done with left multiplication.
    new_grid = new_coord*self.grid
       
#      new_value = new_value.reshape(new_grid.shape())

    return self.copy(value = new_value,grid = new_grid )



  def roll(shift, crd):
    """
    Call roll function on self.

    Args:
      shift: (int) number of index points to roll by
      crd: (Coord) coord to shift on

    Returns:
      Field: containing the shifted value and grid.
    """
    return roll(self, shift = shift,coord = crd)


# belongs to Field 
  def __xor__(self,other):
    """
    Tests the equivalence of the grids of two fields.
    """
    return self.grid^other.grid


  def __pow__(self,n):
    if isinstance(n,int):
      return Field(name = self.name + '**' + str(n),value = self.value**n,grid = self.grid,units = self.units + '^' + str(n))
    else:
      print 'Power error: provide integer.'

  def __neg__(self):
    return Field(name = self.name,value = -self.value,grid = self.grid, units = self.units, direction = self.direction)

  

  def __add__(self,other):
    """
    Field addition F + G. Proceeds only when fields are defined on the same grid. To add fields defined on different grids, use something like F + G(F.grid) or other, depending on the space spanned by the grids.
    If the strict_v attribute of F is set to True (a default), and the direction attributes of F,G differ and are not scalar, addition leads to the formation of a vector Field F*G = (F,G).
    

    """

    L = self.value
    
    if isinstance(other,Field):
      R = other.value
      if L.shape == R.shape:
        if (self.grid == other.grid):
          if not self.grid.same(other.grid):
            warnings.warn('grids contain same data points but different other attributes (e.g. name). Proceeding.')

          if self.strict_v:
            if self.direction == other.direction:
              return Field(name=self.name, value = L+R,grid = self.grid, units = self.units, direction = self.direction)
            else: 
              return self*other

          else:

            return Field(name=self.name, value = L+R,grid = self.grid, units = self.units, direction = self.direction)
        else:
          raise Exception('Field grid error in %s + %s with Field %s: Field grids must be equal. Try F + G(F.grid).' % (self,other,self) )
          
      else:  
        raise Exception('Field shape error in %s + %s with Field %s: shapes must match. Try F + G(F.grid).' % (self,other,self) )
       

    elif isinstance(other,int):
      return self+float(other)
    elif isinstance(other,float):
      return Field(name = self.name,value = self.value + other,grid = self.grid, units = self.units, direction = self.direction)

    else:
        raise Exception('Field type error %s + %s with Field %s: right factor must be Field, int or float.' % (self,other,self) )
        

# --> belongs to Field 

  def __sub__(self,other):
    L = self.value
    
    if isinstance(other,Field):
      R = other.value
      if L.shape == R.shape:
        if (self.grid == other.grid):

          if not self.grid.same(other.grid):
            warnings.warn('grids contain same data points but different other attributes (e.g. name). Proceeding.')

          if self.strict_v:
            if self.direction == other.direction:

# should these Field creation statements be replaced with self.copy?

              return Field(name=self.name, value = L - R,grid = self.grid, units = self.units, direction = self.direction)
            else: 
              return self*other

          else:

            return Field(name=self.name, value = L - R,grid = self.grid, units = self.units, direction = self.direction)


        else:
          raise Exception('Field grid error in %s-%s with Field %s: grids must be equal. Try F - G(F.grid) or F(G.grid) - G.' % (self,other,self) )
          
      else:  
        raise Exception('Field shape error in %s-%s with Field %s: shapes must match. Try F - G(F.grid) or F(G.grid) - G.' % (self,other,self)  )
      

    elif isinstance(other,int):
      return self - float(other)
    elif isinstance(other,float):
      return Field(name = self.name,value = self.value - other,grid = self.grid, units = self.units, direction = self.direction)

    else:
        raise Exception('Field type error in %s - %s with Field %s: right factor must be Field, int or float.' % (self,other,self)  )


  def __setitem__(self,k,v):
    self.value[k] = v
    return



 

  def __getitem__(self,L):

    #getitem of Field .
    # returns a numpy array containing the sliced content of self if argument consists only of slice objects.
    # If argument is of form: (crd0,1,crd2,1:) etc for crd0,crd1 Coord objects, slicing will take place along each Coord using the slice object or integer following each crd argument as the slice object. A new Field will be returned and new associated Coord objects and a corresponding grid will be produced for the return Field.

    # The argument may also contain Ax objects X,Y,Z,T. In this case, the argument will be converted to the corresponding Coord object from the Field grid self.grid via multiplication.
    
    if isinstance(L,tuple):
      # In this case, the argument is expected to be multiple slice objects only or slice objects interspersed with Coord objects.
 
      crds = []		# holds Coord objects along which to slice
      slices = []	# holds slice objects
      
      for i in L:
        if isinstance(i,Coord):
          if i not in self.grid:
            raise Exception('Slice Coord argument %s not in Field %s grid %s.'%(i,self,self.grid))

          crds.append(i)
        elif isinstance(i,Ax):
       
          if i*self.grid is None:
            raise Exception('Slice axis argument %s not in Field %s grid %s.' % (i,self,self.grid))
          else:
            crds.append(i*self.grid) 
        elif isinstance(i,int)  | isinstance(i,slice):
          slices.append(i)
        else:
          raise Exception('Non-integer slice axis argument %s for Field %s not recognised as Ax or Coord object. The Ax/ Coord object might be stale. ' % (i, self) )


      if len(crds) == 0:
        # No Coord objects recorded
        if len(slices) == 0:
          warnings.warn( '(severe): no slices!', RuntimeWarning )
        return self.value[L]
      elif len(crds) == len(slices):        
        # In this case, we can associate a slice object (or int) to each Coord object in the argument.
        # The order then determines which slice object corresponds to which Coord object.
        # The task is now to slice the Field value appropriately and to create the associated Coord objects.

        if len(crds) == 1:
          crd = crds[0]
          if not isinstance(crd,Coord):
            raise Exception('Slice axis not valid. Value crd is: %s '  % crd)

          slc = slices[0]
          new_value = self.slice(sl_coord = crd, slice_obj = slc)

          # Create the new sliced Coord object. By default, this yields an equivalent Coord to the original.

          if isinstance(slc,int):
             # Simple slice case at a certain point along an axis.
            
             return self.copy(value = new_value, grid = self.grid/crd)

          elif isinstance(slc,slice):
            # a subset along an axis is taken. New Coord object(s) with the correct value needs to be created.
            new_crd_value = crd[slc]         
            new_crd = crd.copy(name = crd.name + '_sliced', value = new_crd_value)

            new_crd|crd

            if crd.dual == crd:
              new_crd.dual = new_crd
            else:

              if (slc.stop != None):
                slc_dual = slice(slc.start, slc.stop+1, slc.step)
              else:
                slc_dual = slc

              new_crd_dual = crd.dual.copy(name = crd.dual.name + '_sliced', value = crd.dual[slc_dual])
              new_crd_dual|crd.dual
              new_crd.dual = new_crd_dual
              new_crd_dual.dual = new_crd

            return self.copy(value = new_value, grid = new_crd*self.grid)

          else:
            # Input is neither slice object nor int
            raise Exception('Field slice error in Field %s arg %s : use slice objects only or Coord objects and slice objects.' % (self,L)  )                


        else:
          
          # length of Coord-slice argument pairs is greater than 1.
          # go through arguments recursively.
          F = self
          for e in zip(crds,slices):
            F = F[e[0],e[1]]

          return F
      else:
       raise Exception('Field slice error in Field %s arg %s : use slice objects only or Coord objects and slice objects.' % (self,L)  )         

    else:
      return self.value[L]


  def __call__(self,grid, method = 'linear'):

    return self.regrid(grid = grid, method = method)


  def regrid(self,grid, method = 'linear'):

# this method is very important. 
# If Field T is naturally defined on grid yt*xt, then T(zt*yt*xt) yields a Field with value a 3D array b such that b[k,:,:] = T(yt*xt) for all possible k.

    value = (self.grid(grid, method = method))(self.value)

    if isinstance(value,list):
# in this case the grid argument is a subspace of self.grid so that the grid of the elements is self.grid/grid due to the way self.grid(grid) has been constructed (see call method for grid objects).
      result = []
      for i,e in enumerate(value):
        result.append(Field(name = 'slice_'+str(grid.reverse() )+'_'+str(i) ,value = e, grid =self.grid/grid))
      return result
     
    else:
# the element is probably a numpy array. If not, Field init will throw an error.
      return self.copy(value = value, grid =grid)

  def __mul__(self,other):
    """
    multiplies two Field T1,T2. If T1 is defined on gr1 and T2 on gr2, then T1*T2 is defined on gr1*gr2


    """
    if isinstance(other,int):
      return self*float(other)
    
    elif isinstance(other,float):
      return Field(name = self.name ,value = self.value*other,grid = self.grid, units = self.units ,direction = self.direction)

    elif isinstance(other,Gr):
      # fields commute with Gr objects
      return other*self

    elif isinstance(other,AxGr):
      # fields commute with AxGr objects
      return other*self

    elif isinstance(other,Coord) | isinstance(other,Ax):
#      print 'Warning (benign): converting right multiplicant to Gr from Coord object.'
      return self*(other**2)
    

    elif isinstance(other,Field):
      # both multiplicants are fields
  
      if (self.direction == other.direction) | (other.direction == ID()) | (self.direction == ID()):
        # in this case, at least one of the multiplicants is a scalar (interacting with all directions), or both multiplicants are along the same direction.

    # Note that this multiplication yields precedence for the order of the left multiplicant (self). E.g. (zt*yt*xt)*(xt*yt) = zt*yt*xt
        common_gr = self.grid*other.grid
 
    # This multiplication inflates the values of self and other (arrays) onto the common grid. 
    # In case the grids contain Coord elements that are equivalent but not equal, grid multiplication dictates that common_gr will contain the elements of the left multiplicant (i.e. again a precedence for the left multiplicant). This implies that the right Field will then be interpolated on the left latice
    
        
        if other.name == '':
          new_name = self.name
        elif self.name == '':
          new_name = other.name
        else:
          new_name = self.name +'_times_'+other.name
 
        new_direction = self.direction*other.direction
        if isinstance(new_direction,AxGr):
          new_direction = new_direction[0]

        return Field(name = new_name ,value = self.grid(common_gr)(self.value)*other.grid(common_gr)(other.value),grid = common_gr, units = self.units + other.units, direction = new_direction)

      else:
        if self.grid != other.grid:
          if self.grid&other.grid:

            # If multiplicants are defined on grids that have the same values but are different objects, a duplicate grid is discovered and housekeeping is done. Duplicate grids commonly arise from earlier slicing.
            print 'Duplicate grids. FYI: replacing right gr.'
            del other.grid
            other.grid = self.grid
          else:
            # if grids are different and not duplicates, the resulting vectorfield is likely to be ill defined. Creation proceeds nonetheless, but with a warning.
            warnings.warn( '(severe) VField components defined on different grids.', RuntimeWarning)

        return VField((self,other))

    elif isinstance(other,VField):
      # the right multiplicant is a vector Field.

      if self.direction == ID():
        new_vfield = [self*e for e in other]
        return VField(new_vfield)       
        
     
      elif self.direction in other.direction():
        i = other.direction().index(self.direction)

        new_vfield = list(other)
        new_vfield[i] = self*new_vfield[i]
        return VField(new_vfield)       

      else:
        return VField([self] + list(other) )


    else:
      raise Exception('Field error in %s*%s with Field %s. Provide Field,Gr or Coord objects or int or double for right multiplicant. Hint: common mistake is when multiplying a Field F and a Coord c, and c appears to be in F.grid, c may be stale: check whether they are identical. If not, update c from Exper Coord stack. ' % (self,other,self) )
     

# --> belongs to  Field.
  def __div__(self,other):
    """
    divides two Field T1,T2. If T1 is defined on gr1 and T2 on gr2, then T1*T2 is defined on gr1*gr2


    """
    if isinstance(other,int):
      return self/float(other)
    
    elif isinstance(other,float):
      return Field(name =self.name ,value = self.value/other,grid = self.grid, direction = self.direction)

    elif isinstance(other,Gr):
      return other.mean(self)
    elif isinstance(other,AxGr) | isinstance(other,Ax):
      return self/(other*self.grid)
 
    elif isinstance(other,Coord):
#      print 'Warning: (benign) converting right multiplicant to Gr from Coord object.'
      return self/(other**2)

    elif isinstance(other,Field):

      new_name = self.name
      if other.name == '':
        new_name = self.name
      elif self.name == '':
        new_name = other.name


    # Note that this multiplication yields precedence for the order of the left multiplicant (self). E.g. (zt*yt*xt)*(xt*yt) = zt*yt*xt
      common_gr = self.grid*other.grid
 
    # This multiplication inflates the values of self and other (arrays) onto the common grid. 
    # In case the grids contain Coord elements that are equivalent but not equal, grid multiplication dictates that common_gr will contain the elements of the left multiplicant (i.e. again a precedence for the left multiplicant). This implies that the right Field will then be interpolated on the left latice

      return Field(name = new_name ,value = (self.grid(common_gr)(self.value))/(other.grid(common_gr)(other.value)),grid = common_gr, direction = self.direction)


    else:
      raise Exception('Field error in %s/%s with Field %s. Provide Field,gr or Coord objects or int or double for denominator. (Or check staleness of objects.)' % (self,other,self) )
     

  def vcumsum(self,coord, upward=True):
    """
    Apply vcumsum method of coord on Field.
    """

    return coord.vcumsum(self,upward=upward) 


# IS THIS SUM METHOD BEING CALLED??? IF NOT, REMOVE AND REPLACE WITH COORD BASED METHODS:
# --> method belongs to Field.
  def sum(self,grid=None):

    """
    Computes sum of Field over grid using masked array (nan is not counted). Outputs a float if grid is entire grid of the Field, and a Field on remaining grid (self.grid/grid) if grid argument is a subgrid.
    """

    if not(grid) or self.grid.perm(grid):
# in this case no grid argument is given, or the full grid is given (up to a permutation).
      R = ma.masked_array(self.value,np.isnan(self.value))
      return ma.sum(R)
    else:
# in this case, it is assumed the user wants to take sums along a certain set of axes, where that grid object is a subspace of self.grid  

# obtain the dual vectorspace axes of grid argument, due to the way the call method of grid objects works.   
      F = self(self.grid/grid)
      
# we assume that F is now a list of fields.
# each element has to be summed.
    
      result = []
      for e in F:
        result.append(e.sum())

      new_grid = self.grid/grid 
           
      return Field(name = self.name, value = (np.array(result)).reshape(new_grid.shape()) ,grid = new_grid)



  def ones(self, nan_val = np.nan):

    """
    Returns Field containing domain of this Field: values are 1 in grid locations where Field is defined, nan otherwise.
    """

    new_fld = self.grid.ones()
    new_fld.value[np.isnan(self[:])] = nan_val

    return new_fld

# --> method belongs to Field.
  def dV(self):

    return self.ones()*self.grid.vol()

  def vol(self, grid = None):
    """
    Compute volume (area/ length) of non-nan grid cells. 
    """
    return (self.dV()).sum(grid)
    
  def vsum(self,grid = None):
    return (self.dV()*self).sum(grid)

  def mean(self,grid = None):
    return (self.dV()*self).sum(grid)/(self.dV()).sum(grid)

  



  def slice(self, sl_coord = None,slice_obj = slice(1,None,None)):
    """
    Slice along Coord (e.g. xt) using slice_obj as slice, e.g. slice(1,None,None).
    """
  
    if isinstance(sl_coord,Coord):
     
      sl = slice(*(None,))
      
      L = []
      for e in self.grid:
        
        if sl_coord is e:
          L.append(slice_obj)
        else:
          L.append(sl)    
     
      return self.value[L]

  def draw(self, colorbar = True,**kwargs):

    if len(self.grid) == 1:
      h= plot(self,**kwargs)
      cb = None

    elif len(self.grid) == 2:
   
      h = contourf(self,**kwargs)
      cb = plt.colorbar()
      cb.set_label(self.units)

    elif len(self.grid) == 3:
      for e in self.grid:
        if hasattr(e,'axis'):
          if e.axis.name == 'Z':
            break       
      if e.axis.name != 'Z':
        e = self.grid[0]

      h = contourf(e(self))
      cb = plt.colorbar()
      cb.set_label(self.units)

    
    return h, cb




# ------------------ end Field  definition ----------------



class VField(tuple):

  """
  vector Field. A tuple of fields with extra rules. Allows multiplication.
  """

  def __mul__(self,other):

    if isinstance(other,Field):
      
        if other.direction == ID():
          # scalar Field multiplication works on individual member fields.
          return VField([e*other for e in self])

        elif other.direction in self.direction():
          i = self.direction().index(other.direction)
          new_vfield = list(self)
          new_vfield[i] = new_vfield[i]*other
          return VField(new_vfield)
        else:
          
          new_vfield = list(self)
          if new_vfield[-1].grid != other.grid:
            warnings.warn('VField components defined on different grids.',RuntimeWarning) 

        return VField(new_vfield + [other])

    elif isinstance(other,VField):
        if len(other) > 1:
          return VField( reduce(lambda x,y: x*y, [self*other[0]]+list(other[1:])) )
        else:
          return self*other[0]

    else:
       # all other types will work on the individual fields. Error messages will be generated from individual multiplication.

       return VField([e*other for e in self])


  def __div__(self,other):

    if isinstance(other,Field):
      
        if other.direction == ID():
          # scalar Field multiplication works on individual member fields.
          return VField([e/other for e in self])

        elif other.direction in self.direction():
          i = self.direction().index(other.direction)
          new_vfield = list(self)
          new_vfield[i] = new_vfield[i]/other
          return VField(new_vfield)
        else:
          new_vfield = list(self)
        return VField(new_vfield + list(other**-1))


    elif isinstance(other,VField):
        if len(other) > 1:
          return VField( reduce(lambda x,y: x*y, [self*other[0]]+list(other[1:])) )
        else:
          return self*other[0]

    else:
       # all other types will work on the individual fields. Error messages will be generated from individual multiplication.

       return VField([e/other for e in self])

  def __neg__(self):
    return VField( [-e for e in self] )

  def __sub__(self,other):
    return self + (-other)

  def __add__(self,other):

    if isinstance(other, VField):

      L_l = []
      L_r = []

      for lft in self:
        for rgt in other:
          if lft.direction == rgt.direction:
            L_l.append(lft)
            L_r.append(rgt)

      if len(L_l) == len(L_r):
        L =[]    
        for i, it in enumerate(L_l):
          sum_fld = it + L_r[i]
          # direction must be assigned, as summing fields does not retain direction.
          sum_fld.direction = it.direction
          L.append(sum_fld)
        return VField(L)
      else:
        raise Exception('Error in VField addition %s + %s. Provide equal length' % (self,other))

    elif isinstance(other,Field):
      # sum a Field to a VField. the Field is added to all members.  
#      if other.direction == ID():

    
      L = []
      for it in self:
        sum_fld = it + other
        sum_fld.direction = it.direction
        L.append(sum_fld)
      return VField(L)

    return


  def vcumsum(self,coord, upward=True):
    """
    Apply vcumsum method of coord on Field elements of self.
    """
    return VField([coord.vcumsum(e,upward=upward)  for e in self] )


  def vsum(self,coord):
    """
    Apply vsum method of coord on Field elements of self.
    """
    return VField([coord.vsum(e,upward=upward)  for e in self] )










  def innersum(self):

    return reduce(lambda x,y: x+y, self)


  def copy(self):

    return VField([ e.copy() for e in self ])    

  def direction(self):
    """ 
    Method that returns a tuple of the directions of the tuple components of this vector Field by examining these components.
    """
    return reduce(lambda x,y: x*y, [e.direction for e in self])


  def draw(self, **kwargs):

    if len(self.direction()) == 2:
      if len(self[0].grid) == 2:

        # insert quiver plot here.
        quiver(self)
      elif len(self[0].grid) == 3:

        pass
    else:

      print "Refused. Only plotting 2D fields."



# ------------ end of VField  --------------------







# Field related functions:




def concatenate(fields, ax=None, name_suffix='_cat', new_coord_name = 'gamma', new_coord= None, strings = None ):
  """
  Joins a sequence of fields together.

  concatenate((a1,a2,...),ax=None)
  Parameters
  ----------
  a1, a2,.... : sequence of Field objects
       
  The Field value ndarrays must have the same shape, except in the dimension
  corresponding to `ax` (the one with unequal Coord point values, by default). axis : Ax object, optional       
  The axis along which the arrays will be joined.  Default is the first one with unequal Coord point values.

  a new Coord is created if none of the grid elements point in the direction of the ax argument. Then, new_coord_name is used. 
  The above behaviour is overridden if the new_coord argument is given. This is a Coord object that will be used to construct one Field from the fields list argument. The list elements become slices (at single Coord values) and the new_coord values are the corresponding coordinates.
  """

  
  if fields == []:
    raise ValueError('Provide list of fields.')

  if new_coord is not None:
    if len(fields) != len(new_coord):
      raise ValueError('Provide fields and new_coord arguments of equal length if providing new_coord argument.')

    # EXIT POINT
    return fields[0].copy( name = fields[0].name +name_suffix, value = np.array( [ F[:] for F in fields ] ) , grid = new_coord*fields[0].grid )



  if ax and (ax*fields[0].grid is None):


    # the axis is not in the grid of the first Field
    expanded_fields = []

    if strings is not None:    
      for i, F in enumerate(fields):
        new_coord = Coord(name = new_coord_name,value = np.array([i]), direction = ax.name , axis = ax , strings = [strings[i],] )
        expanded_fields.append( F(new_coord*F.grid) )
    else:
      for i, F in enumerate(fields):
        new_coord = Coord(name = new_coord_name,value = np.array([i]), direction = ax.name , axis = ax  )
        expanded_fields.append( F(new_coord*F.grid) )

    fields = expanded_fields    
    name_suffix = ''

  # EXIT POINT
  return reduce(lambda x,y:x.cat(y,ax=ax, name_suffix = name_suffix), fields)

def squeeze(F, hard = False):
  """
  Equivalent to Numpy squeeze method. Remove dimensions and associated coords in grid of length 1. Reversible operation as squeezed dimensions are recorded. Setting argument hard to True yields an irreversible squeeze where the squeezed dims are not recorded (and cannot be unsqueezed later). 
  """

  dims = list(F.grid)
  body = F.value
  
  squeezed_dims = []

  if hard:
    # irreversible squeeze
    # In this case, the squeezed dims are not recorded for later inflation
    for i,dim in enumerate(dims):
      if body.shape[i] == 1:
        dims.remove(dims[i])

  else:
    # reversible squeeze
    # In this case, the squeezed dims are ecorded for later inflation
    for i,dim in enumerate(dims):
      if body.shape[i] == 1:
        squeezed_dims.append(dims[i])
        dims.remove(dims[i])
      
   
  body = np.squeeze(body)
  return F.copy(value=body,grid = Gr(dims) , squeezed_dims =  Gr(squeezed_dims) )

def unsqueeze(F ):
  """
  Opposite of squeeze. Uses the grid stored in squeezed_dims Field attribute to restore the unit-length dimensions (coords) of the Field. 
  
  """

  gr_unsqueezed = F.squeezed_dims*F.grid

  return F.copy( value = F.value.reshape(gr_unsqueezed.shape() ) , grid = gr_unsqueezed, squeezed_dims =  Gr( () )  )

 




def nugget(path = None, name = None,fields = [] , history = 'Created from Spacegrids '  ):

    """
    Write.

    Creates Netcdf file and writes all loaded Field to it, along with their Coord objects.

    """

    if name is None:
      name = 'nugget'

    if not name.split('.')[-1] in ['nc','cdf']:
      name = name +'.nc'
    if not path is None:
      name = os.path.join( path , name ) 
   
    

    print 'Writing Field to file %s'%name

    try:
      file_handle = netcdf_file(name , 'w')

    except IOError:
      print 'Cannot open %s'%name

    else:

      for fld in fields:

        file_handle = fld.cdf_insert(file_handle)

      file_handle.history = history + '%s'%str(datetime.datetime.now())
 
#    var_cdf.units = self.units

      file_handle.close()




def roll(F,shift=1,coord=None,axis=None,mask=False,keepgrid = False, nan_val = np.nan):

  """
  Function that rolls a Field similar to np.roll on numpy arrays (sg roll actually calls np.roll). Axis can be picked via coord name. If mask is True, the elements that rolled from the other side of the array are set to nan (appropriate for non re-entrant domains). The rolled coord element of the grid belonging to Field F is replaced by a new Coord object reflecting the roll operation. To disable this Coord replacement, use argument keepgrid = True

  

  NOTE: axis here means np array index.

  """

  if isinstance(coord,Ax):
    coord = coord*F.grid

  if not(axis):
    if coord in F.grid:
      
      axis = F.grid.index(coord)
    else:
      print 'coord not in Field grid'
      return 

# avoid deepcopy for fields
# Fr is the rolled Field.

  if keepgrid is True:
    # keep the original grid of Field F
    newgr = F.grid 
  elif keepgrid is False:
    # replace the grid with one with rolled coord
    newgr = coord.roll(shift = shift)*F.grid 
  else:
    raise Exception('Argument error in roll of Field %s. Provide True or False for keepgrid argument. ') % F


  Fr = F.copy(value = np.roll(F.value,shift=shift,axis=axis), grid = newgr )
  
  if mask:
    # handle the areas in the Field that need to be set to nan
    sl = slice(*(None,))
    
    if shift > 0:
     
      sl_exposed = slice(0,shift,None)

    elif shift < 0:
# note that shift is negative here, indicating last elements of array.
      sl_exposed = slice(shift,None,None)


    L = []
    for e in F.grid:
      L.append(sl)
    L[axis] = sl_exposed
    Fr.value[L] = nan_val

  return Fr

def ones(grid):

  return Field('ones',np.ones(grid.shape()),grid)



def finer_Field(F,factor =5.):

  """
  This is a more UVic specific function to prepare a Field containing the outline of the continents for horizontal plots.
  """
  
  return F(finer_grid(grid = F.grid,factor = factor),method ='nearest')






# ------------- some Coord related functions ----------------------



# used in function cdfsniff
cdf_axes = {'X':XCoord,'Y':YCoord,'Z':Coord,'T':Coord,'scalar':Coord}



def make_dual(crd,name = None,guess_append = True,append_last=True, zero_boundary = False):
  """
  Create a dual Coord by appending one entry, of which the width is guessed based on the adjacent cell width.
  """

  if name is None:
    name = crd.name

  if guess_append:
    # Guesses according to CSIRO  model conventions
    if (crd.direction == 'X') | (crd.direction == 'Y'):
      append_last = False 
    else:
      append_last = True      

  if len(crd[:]) > 1:
    if append_last:
      value = np.append(crd[:] , np.array(2*crd[-1] - crd[-2]))

    else:
  
      value = np.append(np.array(2*crd[0] - crd[1]),crd[:])

  else:
      value = crd[:]
 
  return crd.copy(name = name +'_edges',value= value, long_name = crd.long_name + ' as edges')

def find_set_dual(cstack, force = None):
  """
  This function tries to find duals among a list cstack (argument) of Coord objects.

  Checks if duals have been defined before. If one such Coord is found, function is aborted (it is assumed it is not needed then). Override with argument force = True.

  """

  if force is None:
    # Check if duals have been defined before. If one such Coord is found, function is aborted (it is assumed it is not needed then).
    for c in cstack:
      if c.dual != c:
     
        return cstack

  # create grid, and therefore tuple, of all axis objects associated with Coord objects in list cstack.
  axes_available = reduce(lambda x,y: x*y, [c.axis for c in cstack])


  for a in axes_available:
    # L is list of all Coord objects that have same axis.
    L = [c for c in cstack if c.axis == a]
    
    # if this has two elements, these are considered duals.
    if len(L) == 2:
  
    # only unique pairs with the same axis will be interpreted as duals.    
      if guess_grid_type(L[0]) == 'ts_grid':
        crd_dual = make_dual(L[1],name = L[0].name) 
        L[0].dual = crd_dual
      elif guess_grid_type(L[1]) == 'ts_grid':
        crd_dual = make_dual(L[0],name = L[1].name) 
        L[1].dual = crd_dual

    elif len(L) == 1:
       
        crd_dual = make_dual(L[0],name = L[0].name) 
        if len(crd_dual) > 1:
          L[0].dual = crd_dual
        else:
          # in this case, the Coord is made to be self-dual. this could happen for a time Coord with only 1 time slice. make_dual will return a length 1 dual for a Coord of length 1. Note that for UVic data this else clause is NOT needed to make the Coord self dual. 
          L[0].dual = L[0]
        
  return cstack
      




def find_equal_axes(lstack,rstack):

  """
  Expects two lists of Coord objects and determines which Coord objects are equal. This is needed when different Coord objects have identical attributes.

  """

  for lc in lstack:
    for i,rc in enumerate(rstack):
      if (lc.axis == rc.axis):
        # use Coord equality method & (__and__):
        if lc&rc:
          # if all 3 attributes are equal values, replace right stack element with left stack element
          rstack[i] = lc
        else:
          # in this case the Coord elements only have the same axis attribute, and are merely equivalent.
          if not rstack[i]^lc:
            rstack[i] | lc


#  return rstack




# -------- io related Coord functions --------------------


def cdfsniff(path_parent, file_extensions = cdf_file_extensions, verbose = False):
  """
  This sg function looks inside the path_parent path (path to directory containing the Netcdf files, provided as argument) for Netcdf files and extracts Coord objects from the dim data using sg.cdfsniff_helper.

  Returns all Coord objects that contain different data, to be used in the Coord stack cstack.
  """

  if os.path.isfile(path_parent):
    # In this case, a file path is provided. This occurs when experiment object correspond to (Netcdf) files instead of directories containing Netcdf files.
    return rem_equivs(cdfsniff_helper( path_parent , verbose = verbose ))

  # all files within path_parent
  fnames = os.listdir(path_parent)

  # cstack will contain all Coord objects constructed from dims in Netcdf 
  cstack = []

  # prepare glob patterns to look for Netcdf files
  globfpaths = [os.path.join(path_parent , e) for e in file_extensions]

  cdf_filepaths = reduce(lambda x,y: x + y, [glob.glob(e) for e in globfpaths  ]  )

  # construct combined cstack out of individual Netcdf files via cdfsniff_helper: 
  for cdf_filepath in cdf_filepaths:

    cstack = cstack + cdfsniff_helper( cdf_filepath , verbose = verbose )

  # remove equivalent Coord objects (containing the same data) and return    
  return rem_equivs(cstack)


def cdfsniff_helper(filepath, verbose = False):
  """
  Takes inventory of coords in netcdf file. 

  Input:
  filepath	total file path the specific Netcdf file.

  Output:
  A list of spacegrids Coord objects.

  Directions and therefore types of Coord objects (e.g. XCoord) are guessed from description and naming of Netcdf vars.


  """

# axis to the possible axes encountered in netcdf: X,Y,Z
  global cdf_axes

  file = netcdf_file(filepath,'r')

  coord_stack = []
  dimensions = file.dimensions
 
  dimensions = [e for e in dimensions if e in file.variables]

  for dim_name in dimensions:

    # guess which direction the Coord is pointing in, based on netcdf descriptions. The netcdf .axis attribute is included!
    # leave directional_names wild card: no filter on general name (e.g. velocity).

    # maybe rely more on this dictionary in future:
    metadata = {k:file.variables[dim_name].__dict__[k] for k in file.variables[dim_name].__dict__.keys() if k not in ['data','dimensions','_shape','_size']  }

    coord_name = dim_name
    
    direction = guess_direction(file.variables[dim_name],  name_atts = ['axis','long_name','standard_name'] , x_dir_names = coord_dir_names['x_dir_names'], y_dir_names = coord_dir_names['y_dir_names'], z_dir_names = coord_dir_names['z_dir_names'],t_dir_names = coord_dir_names['t_dir_names'],directional_names = '*')

    if direction == 'scalar':
      # double check that this Coord has no direction by looking at dim_name itself.
      if dim_name in cdf_axes:
        direction = dim_name
        coord_name = dim_name + '_crd'
        print 'OK. Inferring direction from dimension name %s itself. Renaming Coord to %s. '%(dim_name,coord_name)

      else:
        warnings.warn('No direction inferred for %s. Guessed direction is scalar.'%dim_name)


    if hasattr(file.variables[dim_name],'units'):
      units = copy.deepcopy(file.variables[dim_name].units)
    else:
      units = '?'
      warnings.warn('(mild): no units assigned to %s in cdfsniff'%dim_name )

    if hasattr(file.variables[dim_name],'long_name'):
      long_name = copy.deepcopy(file.variables[dim_name].long_name )
    else:
      long_name = '?'
      if verbose:

        warnings.warn('No long_name for %s. Assigning %s'%(dim_name, coord_name) )
      long_name = coord_name

    # look only at the keys
    if hasattr(file.variables[dim_name] , 'axis' ):
      # If the netcdf vatiable has an axis attribute, as in UVic, we look for edges to determine the dual and assign an axis attribute. 

      # only edges Coord objects do not have an axis attribute, so edges
      # Coord objects need to be created simultaneously with their dual.
      # note that time coords have axis attributes but not edges (self-dual).

      # Get the netcdf name of the dual variable (the edges, or bounds). Failure signal is None.
      dual_var_name = get_att(file.variables[dim_name], edge_names, fail_val = None)

      if dual_var_name != None:
        # if edges are defined, we create Coord and its dual in pairs
        if file.variables[dim_name].axis in cdf_axes:

          # convert the netcdf name of the dual to an actual cdf variable
          if dual_var_name in file.variables:
          
            dual_var = copy.deepcopy(file.variables[dual_var_name][:] )
          else:

# THIS IS A TEMPORARY FUDGE IN CASE A FILE HINTS AT COORD EDGES BUT DOESN'T STORE THEM:
             dual_var = copy.deepcopy(file.variables[dim_name][:] )           
       
       
          # using call method of Coord object in cdf_axes global

          this_coord = cdf_axes[direction](coord_name, copy.deepcopy( file.variables[dim_name][:] ), axis = copy.deepcopy(file.variables[dim_name].axis ),direction = direction, units = units, long_name = long_name , metadata = metadata)  

          #this_coord = cdf_axes[file.variables[dim_name].axis](dim_name, file.variables[dim_name][:], axis = file.variables[dim_name].axis, units = units)  
          dual_coord = cdf_axes[direction](dual_var_name,prep_dual_array(dual_var),dual = this_coord, axis = copy.deepcopy(file.variables[dim_name].axis ), direction = direction, units = units, long_name = long_name, metadata = metadata)

          this_coord.dual = dual_coord

          coord_stack.append(this_coord)
          coord_stack.append(dual_coord)

        else:
          warnings.warn('Unknown axis.',RuntimeWarning)

      else:
        # this is the case of self-dual objects such as time, so only 1 object needs to be made
        if file.variables[dim_name].axis in cdf_axes:
          coord_stack.append(cdf_axes[direction](coord_name,copy.deepcopy(file.variables[dim_name][:] ), axis = copy.deepcopy(file.variables[dim_name].axis ),direction = direction, units = units, long_name = long_name , metadata = metadata ))

    else:
    # In this case, no axis attribute has been detected.
      if verbose:
        print 'No Netcdf axis attribute detected. Creating attribute from direction guess.'

      if hasattr(file.variables[dim_name] , 'edges' ):
        # if edges are defined, we create coord and its dual in pairs
        if direction in cdf_axes:
          dual_var = copy.deepcopy(file.variables[file.variables[dim_name].edges][:] )

          # using call method of Coord object in cdf_axes global

          this_coord = cdf_axes[direction](coord_name, copy.deepcopy(file.variables[dim_name] [:] ), axis = file.direction,direction = direction, units = units, long_name = long_name, metadata = metadata)  

          #this_coord = cdf_axes[file.variables[dim_name].axis](dim_name, file.variables[dim_name][:], axis = file.variables[dim_name].axis, units = units)  
          dual_coord = cdf_axes[direction]( copy.deepcopy( file.variables[dim_name].edges ),prep_dual_array(dual_var),dual = this_coord, axis = direction, direction = direction, units = units, long_name = long_name, metadata = metadata)


          this_coord.dual = dual_coord

          coord_stack.append(this_coord)
          coord_stack.append(dual_coord)

        else:
          warnings.warning('guessed direction not in cdf_axes, strange!')  
      else:
        # this is the case of self-dual objects such as time, so only 1 object needs to be made
        if direction in cdf_axes:
          coord_stack.append(cdf_axes[direction](coord_name, copy.deepcopy(file.variables[dim_name][:] ), axis = direction,direction = direction, units = units ,long_name =long_name, metadata = metadata))

        else:
          warnings.warning('guessed direction not in cdf_axes, strange!')  
  
  if coord_stack:
    bins = {}
    for ax in cdf_axes:
      bins[ax] = []
      for i,cc in enumerate(coord_stack):
        if ax == cc.axis:
          bins[ax].append(i)
 
    for name in bins:
      for i in range(len(bins[name])-1):
        coord_stack[bins[name][i]] | coord_stack[bins[name][i+1]]

  file.close()

  # if no duals were defined above, the following function call will detect the absence of duals and try to guess them:


  return coord_stack





def guess_helper(desc, guess_names, true_val = None, false_val = None):
  """
  Helper function for guess_direction
  """


#  if reduce(lambda x,y: x| y ,[e in desc for e in x_dir_names]):  
#    return 'X'

  denied_found = [e[1:] in desc for e in guess_names if e[0] =='!']
  if denied_found:
    deny = reduce(lambda x,y: x| y , denied_found)
  else:
    deny = False
  allowed_found = [e in desc for e in guess_names if e[0] !='!']
  if allowed_found:
    allow = reduce(lambda x,y: x| y , allowed_found)
  else:
    allow = False

  if not(deny) and allow:  
    return true_val
  else:
    return false_val

def guess_direction(cdf_var,  name_atts = ['long_name','standard_name'], x_dir_names = ['eastward','Eastward','zonal','Zonal'], y_dir_names = ['northward','Northward','meridional','Meridional'], z_dir_names = ['upward','Upward','vertical','Vertical'],t_dir_names = [],directional_names = ['velocity','stress','momentum flux','momentum_flux']):

  """
  Helper function for cdfread. Used to guess, based on keywords in the netcdf data descriptions, whether a Field is a (space-) vector Field component and in what direction it points. The directional_names argument is a list of keywords that might show up in a description that indicates a vector component: e.g. the word velocity. If this list is empty, the function will not search for those keywords (less restrictive). The name_atts argument indicates the possible name of a descriptive attribute in a netcdf file. The {x,y,z}_dir_names correspond to keywords indicating that particular direction (x,y,z).

  """

 # the keywords in the description will indicate a direction Field.
  # i.e. a vector Field component. If found, their direction attribute will be set in the appropriate direction. Otherwise, it is a scalar.

 
  # loop through the possible names the attribute can have

  for na in name_atts:
 
    if hasattr(cdf_var,na):
    
      desc = getattr(cdf_var,na)     
      
      if reduce(lambda x,y: x| y ,[e in desc for e in directional_names]) | (directional_names == '*'):  
      # directional description found:

        try_dict = {'X':x_dir_names, 'Y':y_dir_names, 'Z': z_dir_names, 'T':t_dir_names}

        for XX in try_dict:

          try_dir = guess_helper(desc, try_dict[XX], true_val = XX)
        
          if try_dir:
            return try_dir
       
 
  return 'scalar'


def cdfread(filepath,varname,coord_stack=[], ax_stack = [], verbose = True,squeeze_Field=False):
  """
  Reads data corresponding to variable name varname from netcdf file. Returns Field object. coord_stack is used to provide Field with grid object built from corresponding Coord objects according to information in netcdf.
  Input filepath is complete path pointing to file.

  """

  file = netcdf_file(filepath,'r')
 
  if varname not in file.variables:
    if verbose:
      warnings.warn('(moderate) from cdfread: var name not in file.')
    return None

  var_cdf_ob = file.variables[varname]

  # in future we are going to use this metadata instead of below attributes. For now, it is used when fields are saved.
  metadata = {k:var_cdf_ob.__dict__[k] for k in var_cdf_ob.__dict__.keys() if k not in ['data','dimensions','_shape','_size']  }

  dims = list(var_cdf_ob.dimensions)  
  body = copy.deepcopy(var_cdf_ob[:])
  if isinstance(body,np.ma.masked_array):
    # The Netcdf4 module yields masked arrays. Convert to ndarray to work well with sg.
    # fill_value is a standard attribute of masked arrays (no checks):
    mis_val = body.fill_value

    body = np.array(body)
  else:
    fvn = get_att(var_cdf_ob, fval_names,fail_val = [np.nan])
    if isinstance(fvn,list):  
      mis_val = fvn[0] 
    else:
      mis_val = fvn

#  mis_val = var_cdf_ob.missing_value[0]
  if hasattr(var_cdf_ob,'units'):
    units = var_cdf_ob.units
  else:
    units = '?'

  if hasattr(var_cdf_ob,'long_name'):
    long_name = var_cdf_ob.long_name
  else:
    long_name = varname

 
# attempts at interpreting data by obtaining the string name of the direction. these names are a convention: 'X','Y','Z'. This direction guess is for fields that could be components of a vector fields.

  direction = guess_direction(var_cdf_ob)

  Dict = {e.direction:e for e in ax_stack}
  Dict['scalar'] = ID()

  direction = Dict[direction]


  if squeeze_Field:
  # if there are dimensions of length 1, remove them.
    if 1 in body.shape:
   
      for i,dim in enumerate(dims):
        if body.shape[i] == 1:
          dims.remove(dims[i])
   
    body = np.squeeze(body)

  try: 
    body[body == mis_val] = np.nan
  except:
    
    warnings.warn('missing value not set to NaN.')
 
  grid = []

  for dim in dims:
    dim_val = file.variables[dim][:]
   
    for crd in coord_stack:
   
      if (dim == dimname(crd) ) and np.array_equal(dim_val , crd[:]):
     
        grid.append(crd)
            
  file.close()
#  print '-----'
#  print Gr(tuple(grid)).shape()
#  print body.shape
 
  return Field(varname,body,grid=Gr(tuple(grid)),units = units, direction = direction, long_name = long_name, metadata = metadata)


def dimname(crd):
  """
  Strip off _crd suffix if it had been added because dim name equalled axis name
  """
  if ('_crd' in crd.name) and (crd.name.split('_crd')[0] == crd.axis.name):
    return crd.axis.name
  else:
    return crd.name


def delta(x):
  """
  Function that calls the d method of the Coord object depending on the kind of coordinate (i.e. x or y).

  """
  if isinstance(x,XCoord):
    for oth in x.others:
      if isinstance(oth,YCoord):
        return x.d(oth)
    print 'No y found.'
    return
  else:
    return x.d()





def prep_dual_array(raw_array):

  if raw_array.ndim == 1:
    return raw_array
  elif raw_array.ndim == 2:
    if raw_array[0,0] == raw_array[1,1]:
      # assume FAMOUS-type encoding of coord edges.
        new_array = np.array([raw_array[0,1],] + list(raw_array[:,0]) )
    elif raw_array[0,1] == raw_array[1,0]:
      # assume FAMOUS-type encoding of coord edges.     
        new_array = np.array(list(raw_array[:,0]) + [raw_array[-1,1],] )
    else:
      warnings.warn('! edge/ boundary var of ndim 2 not recognised.')
 
      print raw_array
      # (un!)lucky guess:
      new_array = raw_array[:,0]    

  else:
    warnings.warn('! edge/ boundary var not recognised.')
    new_array = raw_array    
    
  return new_array


# ------------- some Ax related functions ----------------------

def guess_grid_type(crd, default = 'ts_grid'):
  """
  Function that guesses the grid type using the keywords contained in the (sg) global dictionary grid_type_names by testing for keywords (contained as lists in that dictionary).

  Returns None if no grid type is found. 
  """
  for grtn in grid_type_names:
    if reduce(lambda x,y: x|y, [e in crd.long_name for e in grid_type_names[grtn]]):
      return grtn
  # No grid type found, default is assumed.
  return default






def make_axes(cstack):
  """

  Replaces axis attribute of Coord objects if it is a string with newly created (non-repeating) corresponding axis objects.

  inputs: cstack, a list of Coord objects.
  outputs: returns a list of all unique (no repeats) Ax objects that have been created to replace the axis attribute of the elements of the Coord list (cstack) argument that were strings. None is returned when all cstack Coord elements already have Ax axis attributes.


   The Ax objects are created here.

  Returns list of Ax objects!!

  NOTE THAT THIS FUNCTION DOES 2 THINGS: IT RETURNS A LIST OF AXES AND MODIFIES THE CSTACK ARGUMENT. 

  """
  # No Coord objects will be removed from the cstack list. But cstack argument is modified!
  # created_axes will contain the newly created Ax objects!! So created_axes is NOT cstack!!
  created_axes = []
  for c in cstack:
    if hasattr(c,'axis'):
      # string Ax attribute will be replaced with corresponding Ax object attribute
      if isinstance(c.axis,str) or isinstance(c.axis,unicode):
        # str attr found --> create corresponding Ax object.
        new_ax = Ax(c.axis,direction = c.direction)
        if id_in(created_axes,new_ax):
          # however, if we already have that Ax object in created_axes, assign existing Ax object instead.
          i = id_index(created_axes,new_ax)
          if c.direction != created_axes[i].direction:
            # if direction inconsistent, proceed but issue warning.
            warnings.warn('Warning! While associating Ax object %s with coord %s. %s.direction =%s, %s.direction = %s ' %(created_axes[i], c, created_axes[i], created_axes[i].direction, c, c.direction) )
          c.axis = created_axes[i]
          
        else:
          # new Ax object not in created_axes, proceed using new Ax object.
          created_axes.append(new_ax)
          c.axis = new_ax 
        # Ax object equivalent (parallel) to Coord object.
        c.axis | c     
  # return the list of Ax objects. 

  # there might be dual coords associated with coords that still have an axis attribute that hasn't been converted from str yet. Replace them with elements from the cstack:
  for j, crd in enumerate(cstack):
    if hasattr(crd,'dual'):
      if not crd.dual is crd:
        i = id_index(cstack, crd.dual)
        if i is not None:
          
          cstack[j].dual = cstack[i]

  return created_axes        




# define the scalar axis via a lazy class. Then ID*X = X etc.
class GetId(object):
  
  def __init__(self):
    self.ret_id = None  

  def __call__(self):
    if not(self.ret_id):
      self.ret_id = Ax('scalar')


    return self.ret_id

ID = GetId()


# ---------------- Gr related functions --------------------

def finer_grid(grid, factor = 5.):

  return reduce(lambda x,y: x*y, [crd.finer(factor = factor) for crd in grid])






from plotting import *




