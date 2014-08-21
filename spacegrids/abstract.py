#encoding:utf-8

""" Abstract classes to build Coord, Gr, Ax, AxGr and Field class.

Most classes inherit from the abstract base classes Named, Associative, Directional, Membered, Valued contained in abstract.py. 

The abstract module contains the following classes:

Named
-----

Base class for most other sg classes, representing objects with copy and same methods.

Associative
-----------

Associative class that objects with the equiv method can belong to. Two objects will be equivalent if they belong to the same associative class. 

Directional
-----------

Base class for derived Coord and Ax classes, representing "direction" (e.g. "latitude" or "depth"). An abstract equivalence relationship is defined among Directional objects where two objects are equivalent when they have the same 'associative' attribute (pointing to an Associative object). This relationship is generally used to indicate whether two Directional objects have the same direction (e.g. X,Y), but could represent other relationships depending on the user.

Membered
--------

Base class for classes containing members such as a grid (Gr) object containing coordinate (Coord) members, or an AxGr object containing Ax objects, e.g. (X, Y).

Valued
------

Base class for classes that contain a ndarray value attribute. The Field class is derived from this.
"""

import numpy as np
import inspect
import copy
import warnings
from utilsg import *
from _config import *

from decorators import check_equiv, method2members, att2members

warnings.formatwarning = warning_on_one_line




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

  def __repr__(self):
    return self.name

  def copy(self, *args, **kwargs):
    """
    Copy method for Named. See __init__ for arguments.

    Most child classes should inherit this method.

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

  def json(self, types_allow = []):
    """convert self to a json friendly object.

    Usage: json.dumps(X.json())
    """
    types_basic = ['int','double','float','str']

    
    # class encoding doesn't work yet:
    return_dict = {'class':str( type(self)  )}

    
    for k in self.__dict__:

      ob_type = type(self.__dict__[k])
      if ob_type not in types_basic:

        if ob_type in types_allow:
          # this is the nested case on which to call the method recursively

          return_dict[k] = self.__dict__[k].json(types_allow= [t for t in types_allow if t != ob_type ] )

        else:
          if isinstance(self.__dict__[k] , np.ndarray ):
            return_dict[k] = self.__dict__[k].tolist()
          else:

            return_dict[k] = self.__dict__[k].__repr__()



    
    return return_dict



class Associative(Named):
  """
  Associative class that objects with the equiv method can belong to. 

  Two objects will be equivalent if they belong to the same associative class. 

  For Coord objects, this should remain consistent with the axis attribute: two Coord objects belong to the same associative class iff they have the same axis attribute. In this case, the associative class of Coord objects is effectively their direction or axis. The mechanisms for the two remain independent.

  For classes using Associative as equivalence principle:

  Their copy method should carry over the Associative object of the parent.
  Their make_equiv method should make the associate the Associative object of the argument equal to the Associative object of the calling object.
  Their __init__ method should create a new Associative class as default behaviour, and assign an argument Associative class if given.
  """

  def __init__(self,name):

    self.name = name
    self.associative = self


class Directional(Named):
  """
  Base class for derived Coord and Ax classes, representing "direction" (e.g. "latitude" or "depth").

  An abstract equivalence relationship is defined among Directional objects where two objects are equivalent when they have the same 'associative' attribute (pointing to an Associative object). This relationship is generally used to indicate whether two Directional objects have the same direction (e.g. X,Y), but could represent other relationships depending on the user.

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

  def __init__(self,name='scalar',direction ='scalar',long_name= '', associative = None ):  
    """
    Initialisation of Directional object. 

    Args:
      name: (str) name of Object
      direction: (str) name of direction in which object points
      long_name: (str) longer description (e.g. for display or in Netcdf)
      associative: (Directional or Associative) object that this new object is equivalent to, or its Associative

    Returns:
      Directional object

    Raises:
      ValueError if associative has no associative attribute
    """
  
# choosing the name ID creates an identity object. ID*b = b for all Coord elements b.
# could implement the identity Field in __call__

# Metric could be a class. Objects of this class could be constructed by a method of the Coord class (Coord objects then spawn metric objects).
    self.equivs = [self]
    self.name = name
    self.direction = direction
    self.long_name = long_name

    if associative is None:
      self.associative = Associative(self.name+'_assoc')
    else:
      if hasattr(associative, 'associative'):
        self.associative = associative.associative
      else:
        raise ValueError('provide object with associative attribute for associative.')

  def __neg__(self):
    """
    To be overriden for now. Could introduce +-1 here.
    """
    pass

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
    Register equivalence of two Directional objects. 

    Args:
      other: (Directional)

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

    other.associative = self.associative  

    return


  def is_equiv(self,other, checks = False):
    """
    Test for equivalence (under make_equiv) between Directional objects. e.g. xt is equivalent to xu

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
    
#    if (other in self.equivs) | (self in other.equivs):
    if self.associative is other.associative: 

      return True     
    else:
      # Warnings helpful for debugging and spotting potential problems:

      if checks is True:

        try:
          if ( self.same(other) ):
            warnings.warn('Warning (severe): %s.is_equiv(%s) is False, but %s.same(%s) is True! ' % (self,other,self,other) )
          elif ( self.weaksame(other) ):
            warnings.warn('Warning (severe): %s.is_equiv(%s) is False, but %s.weaksame(%s) is True! ' % (self,other,self,other) )
        except: 
          warnings.warn('Consistency check failed.')

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
  Base class for classes containing members such as a grid (Gr) object containing coordinate (Coord) members, or an AxGr object containing Ax objects, e.g. (X, Y).

  This class is intended for multiple inheritance with classes that provide container functionality. For example Tuple. The elements inside the tuples are then referred to as "members", hence the name of this class.  Methods relate to general operations with these members. Note that more specific member-related methods are relegated to derived classes.

  __init__ method of joint-inheritance class must take a container of elements. 
  """


  def same(self,other):
    """
    Member-wise same comparison.
    """

    if len(self) == len(other):
      for i,c in enumerate(self):
        if not(c.same(other[i]) ):    
          return False

      return True
    else:
      return False



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


  def call_on_members(self, method, *args, **kwargs):
    """
    Call method on all members and construct new Membered object.
    """

    return self.__class__( [ getattr(member, method)(*args, **kwargs) for member in self ] )


  def get_from_members(self, att_name):
    """
    Call method on all members and construct new Membered object.
    """

    return self.__class__( [ getattr(member, att_name) for member in self ] )


  def __and__(self,other):
    """
    Shorthand to member-wise weaksame comparison.
    """

    return self.weaksame(other)

  @method2members
  def __neg__(self):
    pass

#    """
#    Call __neg__ on members and return corresponding Membered object.
#    """
#    return self.__class__( [-member for member in self] )

  def reverse(self):
    """
    Reverse the order of the members.

    Examples:

    >>> coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) )
    >>> coord2 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]) )
    >>> (coord1*coord2).reverse()
    (test2, test1)
    """

    return self.__class__([ self[len(self) -i -1 ] for i in range(len(self))  ])


  @method2members
  def copy(self,*args,**kwargs):
    """
    Member-wise copy method.
    """
    pass


  def strict_equiv(self, other):
    """
    Tests whether two Membered objects have equivalent Members at each position.
    This is a stricter test than Membered equivalence testing via gr1.is_equiv(gr2), which only tests whether both Membered objects describe the same linear space (elements equivalent up to a permutation).

    """
    if len(self) == len(other):
      
      RO = True      
      for i,it in enumerate(self):
        RO *= (it.is_equiv(other[i] ) )
      return bool(RO)
    else:
      return False


  def __xor__(self, other):
    """
    Shorthand for is_equiv.
    """

    if len(self) == len(other):
      
      return self.is_equiv(other)


  def is_equiv(self, other):
    """
    Checks member-wise equivalence between Membered objects up to a permutation. 

    For grids,  objects are equivalent if they define the same physical subspace, based on the equivalence definition for Coord classes. In other words, checks whether the individual Coord elements of the two grid (Gr object) arguments are equivalent up to a permutation. A stricter version of this test is strict_equiv, which allows no permutation.
 
    Args:
      other: (Membered) the object to compare with

    Returns:
      True if all (self) members are equivalent to a member of other and vice versa and both have equal length. False otherwise
    """

    if len(self) == len(other):
      
      if self.eq_perm(other):
        return True
      return False

    else:
      return False



  def eq_in(self, member):
    """ Determines whether argument is equivalent to any of the constituent members.

    
    Args:
      member: (Membered) object to be checked.

    Returns:
      True when member is equivalent to one of the member objects, False otherwise.

    See also:
    eq_in method of Coord 
    """


    for i in self:
      if member.is_equiv(i): return True
    return False

  def eq_index(self,member):
    """
    Returns index of argument in members.
    """

    for i,v in enumerate(self):
      if member.is_equiv(v): return i
    return -1

  def rearrange(self,permutation):
    """
    Rearranges the order of the members of this object via permutation arrgument. 

    Args:
      permutation: (List or Tuple) permutation to rearrange by

    Returns:
      object of same type as self with member rearranged.
  
    Examples: 

    >>> g1 = latitude*depth
    >>> g1.rearrange( (1,0) ) 
    (depth, latitude)

    See also Gr.perm method    
    """
    return self.__class__((self[i] for i in permutation))
    
    


  def perm(self, other,verbose = False):     
    """
    yields permutation of axes going from self to other.

    E.g. for grids gr1 and gr2, g2 = g1.rearrange( g1.perm(g2) )

    Returns None if no permutation exists.

    See also rearrange.
    """  

    return find_perm(self,other,verbose = verbose)

  def eq_perm(self, other, verbose = True):      
    """
    Yields permutation of members going from self to other, where equivalent members are treated as identical. 

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





  def json(self, types_allow = []):
    """convert self to a json friendly object.

    Usage: json.dumps(X.json())
    """
    types_basic = ['int','double','float','str']


    members = [member.json(types_allow=types_allow) for member in self]

    return_dict = {'class':str(self.__class__),'members':members}

    for k in self.__dict__:

      ob_type = type(self.__dict__[k])
      if ob_type not in types_basic:

        if ob_type in types_allow:
          # this is the nested case on which to call the method recursively

          return_dict[k] = self.__dict__[k].json(types_allow= [t for t in types_allow if t != ob_type ] )

        else:
          if isinstance(self.__dict__[k] , np.ndarray ):
            return_dict[k] = self.__dict__[k].tolist()
          else:

            return_dict[k] = self.__dict__[k].__repr__()


    return return_dict


class Valued(Named):
  """
  Base class for classes that contain a ndarray value attribute. 

  This class derives its name from the presence of an attribute named "value" that contains a Numpy ndarray. Methods relate to this attribute.
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


  def get_value(self,i):
    return self.value[i]


  def set_value(self, value):
    self.value = value


  def __getitem__(self,i):
    """Obtain item from value atttribute.
    """
    return self.get_value(i)


  def __setitem__(self,value):

    return self.set_value(value)

  def sliced(self,slice_obj = None ,suffix = '_sliced'):
    """Create new sliced Valued object with sliced value.

    The slice argument must match the value dimensions, there are no checks.

    Args:
      slice_obj: (slice objects or tuple of) to slice value with
      suffix: (str) suffix to use for sliced Valued object

    Returns:
      Valued object containing sliced value
    """
    return self.copy(name = affix(self.name, suffix) , value = self.value[slice_obj]  )

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



  def __pow__(self,n):
    """
    Repeated multiplication of object with itself.
    """
    return reduce(lambda x,y: x*y, n*[self])

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



