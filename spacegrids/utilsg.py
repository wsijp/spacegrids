#encoding:utf-8

""" Useful general functions and utilities.
"""

import numpy as np
import types
import math

import fnmatch
import itertools

from _config import *

class Report():
  """
  Reporting class containing useful method for building a string containing messages.

  Attributes
    
    value : `str` The reporting string being built
 
  Examples:
 .. doctest::

  >>> import spacegrids as sg
  >>> REP = sg.Report()
  >>> REP.echoln('---')
  >>> REP.echo('a test')
  >>> REP.line()
  >>> D ={'a':1,'b':2,'c':3,'d':5}
  >>> REP.table(D)
  >>> REP
  ---
  a test
  ----------

  a               1              c               3              
  b               2              d               5              
  """

  def __init__(self,value = ''):
    """ Initializes Report instance.

    Args:
         value: string containing initial string value of Report

    Returns:
         A Report object with attribute value containing a str.
    """

    self.value = value

  def __repr__(self):
    return self.value

  def line(self,char = '-',times = 10):
    """
    Insert a line of --- characters on a new line into Report value.

    Args:
         char: string containing character used for line (default '-')
         times: the number of times char is repeated
    Returns:
         None. 
    """

    if not self.value[-1] == '\n':
      self.echoln()

    self.echoln(char*times)

  def block(self,L,delim = '\ ', width = 12, cols = 4 ):
    """Creates string to display list of strings as a deliminated rows x cols block.

    Args:
         L: list of strings
         delim: delimiter character passed to block method
         width: width in characters of each cell (int)
         cols: number of columns of block (int)
    Returns:
         None
    """


    Ls = _split_list(L , size = len(L)/cols + 1)
    Ls = _transpose_list(Ls)
 

    for LL in Ls:
      self.echo(what = delim.join([ ("%-{}s".format(str(width)))%it for it in  LL ] ) )
      self.echoln()

  def echoln(self,what = '', delim = ' ', maxlen = 10, width = 12, cols = 4):
    """
    Insert a string on a new line into Report value.

    Args:
         what: string to insert (default '') or list/ tuple of strings
         delim: delimiter character passed to block method
         maxlen: max length (int) passed to block method
         width: int passed to block method
         cols: int passed to block method
    Returns:
         A Report object with updated value attribute.    

    This is the echo method with a '\n' appended.

    **See Also**
    echo method
    """

    self.echo(what = what , delim = delim, maxlen = maxlen, width = width,cols = cols)
    self.echo('\n') 


  def echo(self,what='' , delim = ' ', maxlen = 10, width = 12, cols = 4):
    """
    Insert a string into Report value.

    Args:
         what: string to insert (default '') or list/ tuple of strings
         delim: delimiter character passed to block method
         maxlen: max length (int) passed to block method
         width: int passed to block method
         cols: int passed to block method
    Returns:
         None

    If 'what' is a list of strings, the block method will be used to embed a table of the strings into the Report value.

    **See Also**
    echoln method
    """


    if isinstance(what,str) or isinstance(what,unicode):
      self.value = self.value + what
    elif isinstance(what,list) or isinstance(what,tuple):
      if len(what) > maxlen:
        # too many strings to display, truncate.
        self.block(L = what[:maxlen/2] + ['...',] +  what[-maxlen/2:] , delim = delim, width = width, cols = cols)

      else:
        # display entire list
        self.block(L = what, delim = delim, width = width, cols = cols)
     
  def table(self,D,cols = 4):
    """ Create a table based on dictionary D.


    Args:
         D: dictionary of name:value pairs to be included in table
         cols: number of columns of table (int)

    Returns:
         None
    """

    for i, k in enumerate(D.keys()):
      if (2*i%cols == 0):
        self.echoln()

      self.echo('%-15s %-15s' % (k , D[k]))

    

  def __add__(self,other):
    """
    Creates new Report object with str values concatenated.
    """

    return Report(value = self.value + other.value)



def plural(n):
  """Determine whether to use plural of a word describing a qty n of something. 

     Used in stdout messages that mention quantities of something.

    Args:
         n: int or float, the quantity of the thing that is reported on.

    Returns:
         A string that is '' or 's', indicating a grammatical plural or singular.
  """

  if int(n) == 1:
    return ''
  else:
    return 's'


def _split_list(L, size = 10):
  """
  Very general function splitting a list L into a list of sublists of equal length (and a remainder). Used in text formatting functionality.
  """

  length = len(L)
  full_blocks = length/size

  new_L = [ L[i*size:(i+1)*size]   for i in range(full_blocks) ]
  
  if full_blocks * size != length: 
    new_L.append(L[size*full_blocks:])  

  return new_L

def _transpose_list(L):
  """
  Transpose a list of lists. Used in text formatting functionality.

  Args:
      L: list of lists to be transposed.

  Returns:
      A list of lists

  """
  if len(L) > 0:

    new_L = [[] for e in L[0] ] 
 
    for LL in L:
      for i, e in enumerate(LL):
        new_L[i].append(e)

    return new_L
  else:

    return L



def merge(A1,A2, sort = True):
  """
  _merge two ndarrays or iterables and order them. 

  Args:
       A1: first ndarray
       A2: second ndarray to be combined with first array
       sort: flag. sort called on result if True.

  Returns:
       new_list: sorted list of combined values.

  Used to _merge Coord values.
  """

  L1 = list(A1)
  L2 = list(A2)

  newlist = []

  while L1 and L2:
    if L1[0] == L2[0]:
      newlist.append(L1.pop(0))
      L2.pop(0)
    elif L1[0] < L2[0]:
      newlist.append(L1.pop(0))
    else:
      newlist.append(L2.pop(0))

  if L1:
    # if there are remaining non-popped elements of L1, append them (at end)
    newlist.extend(L1)
  if L2:
    # if there are remaining non-popped elements of L2, append them also (at end)
    newlist.extend(L2)

  if sort:
    newlist.sort()

  return np.array(newlist)


def sublist(L,pick_vals):
  """
  Picks a sublist from a list L of strings (e.g. names) using the filematch function. 

  Args:
       L: list of strings from which to pick
       pick_vals: pattern or list of patterns used in fnmatch.fnmatch
 
  Returns:
       sL: a list of matching values

  For instance, a wildcard will pick all elements.
  """

  if not isinstance(pick_vals,types.ListType):
    pick_vals = [pick_vals]
  
  sL = []
  
  for val in pick_vals:
    for el in L:
      if fnmatch.fnmatch(el,val):
        sL.append(el)

  return sL	




def add_alias(L):
  """Create alias attribute based on name attribute, and identify objects that have identical name attributes to yield corresponding numbered alias attributes to tell them apart.

  Args:
       L: list of Coord objects (or any objects with a name attribute).
  
  Returns: 
       A list -- the list L with alias attributes added

  Takes list L of objects with name attribute. Assigns alias attribute.
  Some of these names might be the same. To yield unique identifiers, a unique
  alias attribute can be assigned to the objects, to be used in practice instead
  of the names (as the names might need to correspond to the netcdf names). E.g. ('latitude','latitude','latitude','longitude','depth') yields aliases: ('latitude','latitude2','latitude3','longitude','depth')

  This is useful when Coord names are used interactively in the namespace of
  e.g. ipython, and allows the user to tell the Coord objects that have the 
  same name apart.  Coord and Gr __repr__ methods will display alias attribute 
  instead of name if the alias attribute is available.

  Examples:
 
  For example, we could apply add_alias to a list of 4 test Coord objects, two of which have the same name.
   
.. doctest::

  >>> import spacegrids as sg
  >>> import numpy as np
  >>> coord1 = sg.Coord(name = 'test',direction ='X',value =np.array([1,2,3]) )

  >>> coord2 = sg.Coord(name = 'test',direction ='Y',value =np.array([1,2,3,4]) )

  >>> coord3 = sg.Coord(name = 'test3',direction ='X',value =np.array([5,1,2,3,4]))

  >>> coord4 = sg.Coord(name = 'test4',direction ='X',value =np.array([5,1,2,3,4]))

  >>> L = sg.add_alias([coord1, coord2, coord3, coord4])  
  
  >>> print [it.alias for it in L]
  ['test', 'test2', 'test3', 'test4']

  """

  # construct list of names of the objects. Can be non-unique.
  L_names =[e.name for e in L]

  # Create dictionary of name vs the number of times it occurs in L
  counts = dict((name,L_names.count(name)) for name in L_names)

  for base_name in counts:
    # for each name that occurs, create list L_name of corresponding objects
    L_name = [e for e in L if e.name == base_name ]

    for i,e in enumerate(L_name):
      if i == 0:
        # The first element takes the name of the element for alias
        e.alias = e.name 
      else:
        # Any possible further elements by the same name are suffixed with a number
        e.alias = e.name +str(i+1)

  return L


def _complete(comp_sets):
  """
  Assigns attribute others to every element of a list of lists that is a list of all elements in the own group except the element itself. 

  Args:
       comp_sets: a list or tuple of tuples (or lists).
  
  Returns: 
       None

  E.g. 
  ((xt,yt,zt),(xu,yu,zw))
  It then assigns an attribute others to every element xt,yt,...  that is a tuple containing all other elements in the subset except the element itself.

  i.e. xt.others = (yt,zt)
  """

  for set in comp_sets:
    
    for co in set:
      co.others = tuple([e for e in set if e != co])

  return

def _flat_list(L):
  """Documentation to come.
  """
  W = list(itertools.chain(*L))
  while isinstance(W[0],list):
    W = list(itertools.chain(*W))

  return W


def _rav_index(L,sh):
  """Documentation to come.
  """
  lsh = list(sh)
  lsh[-1] = 1
  if isinstance(L,int):
    return L
  elif isinstance(L,tuple) or isinstance(L,list):
    if len(L) == len(lsh):
   
      result = 0
      coef = 1
      for e in sh:
        coef *= e

      for i in range(len(L)):
        coef = coef/sh[i]

        result += coef*L[i]

      return result
        

    else:
      print 'provide equal length'
      return
  else:
    print 'provide list or tuple'
    return 




def find_perm(left,right, verbose = True):
  """
  Find a permutation to obtain tuple (iterable) argument "right" from "left".

  Args:
       left: a tuple (or list)
       right: a tuple (or list) 
       
  Returns: 
       None or perm, a permutation such that [left[i] for i in perm] = right

  If left and right are of equal length and contain the same elements up to their order, a permutation will be returned.
  """

  if len(left) == len(right):
    perm = []
   
    for r in right:
      if r in left:
        perm.append(left.index(r))
      else:
        if verbose:
          print 'Inputs not permutable: %s'%r.name
        return

  else:
    if verbose:
      print "message from find_perm: inputs must be of equal length."
    return 
  return tuple(perm)



def round_order(value, order = None):
  """
  Round value up to order of magnitude order.

  Args:
       value: float
       order:  (int) the order to which to round value (e.g. order=2 is 100)
       
  Returns: 
       value rounded to the order 'order'

  Examples:

.. doctest::

  >>> import spacegrids as sg
  >>> sg.round_order(1140.,2.)
  1100.0
  >>> sg.round_order(1140.,1.)
  1140.0
  """

  if order is None:
    order = math.log10(abs(value))//1.
  
  return (10**order)*round(value*10**-order)


def order_mag(val):
  """
  Find order of magnitude of value.

  Args:
       val: float
       
  Returns: 
       order of magnitude of val

  Examples:

.. doctest::

  >>> import spacegrids as sg
  >>> sg.order_mag(1140.)
  3.0
  """
  
  if val == 0.:

    return 0.
  else:
    return math.log10(abs(val))//1.




def _which_att(obj,att_list,fail_val = None):
  """
  Returns first element of list of names if it matches an attribute name of object, None otherwise.

  Args:
       obj: an object of which to examine the attributes
       att_list: an iterable containing strings
       fail_val: value to return if att_list elements not in object attribute names
       
  Returns: 
       attribute name that gives firs match in list, fail_val otherwise

  **See also**
  get_att
  """

  for attname in att_list:

    if hasattr(obj,attname):
      return attname
  
  return fail_val


def get_att(obj, att_list,fail_val = None):
  """
  Returns attribute value corresponding to first element of a list of name strings that matches any attribute name of object, None otherwise.

  Args:
       obj: an object of which to examine the attributes
       att_list: an iterable containing strings
       fail_val: value to return if att_list elements not in object attribute names
       
  Returns: 
       attribute value that gives firs match in list, fail_val otherwise

  **See also**
  _which_att
  """

  this_att = _which_att(obj, att_list)
  if this_att is None:
    return fail_val
  else:
    return getattr(obj, this_att)



def id_index(L,e):
  """
  Same as index method of lists, but with equality defined as &.

  Args:
       L: list of objects that implement &-equality
       e: an object of that same type.
       
  Returns: 
       None or integer representing the index of e in L.

  Used on both Ax and Coord objects in code. In Coord case, the & method is identical to the weaksame method.

  **See also**
  id_in
  rem_equivs
  __and__ methods of Ax and Coord class
  List method 'index'
  
  """

  for i,it in enumerate(L):
    if (it&e):
      return i

# not raising error as normally it would. 
#  raise ValueError('%s is not in list' % e)

# return None instead:
  return

def id_in(L,e):
  """
  Same as in function for lists, but with equality defined as &.

  Args:
       L: list of objects that implement &-equality
       e: an object of that same type.
       
  Returns: 
       True or False depending on whether there is an element e2 of L such that e&e2 is True


  Used on Ax objects in code.

  **See also**
  id_index
  rem_equivs
  List method 'index'
  """

  for it in L:
    if (it&e):
      return True
  return False

def rem_equivs(L):
  """
  Removes all all elements from list that are equivalent under &-relationship

  Args:
       L: list of objects that implement &-relationship
      
  Returns: 
       A list of elements that are unique with respect to the &-relationship
      

  Uses function id_in.

  **see also**
  id_in
  id_index

  """

  L_new = []

  for e in L:
    if not(id_in(L_new, e)):
      L_new.append(e)

  return L_new



def read(fname):
  """
  Open fname (full path) for simple reading and return data.
  """

  data = None

  with open(fname,'r') as f:

    data = f.read();

  return data


def _parse_control_line(line):
  """
  Interpret the data within one named block of control file, including the name of the block.

  Args:
       line: string containing read data occurring inside a named block	
       

  Returns: 
       (block_name, pairs) where pairs is a list of the (name, value) pairs interpreted from the control file and block_name is the name of the block.  

  Helper function to function parse_control_file
  """

  if line.replace(' ','') == '':
    return

  block_name = line.split(' ')[0]
  line = line.replace(block_name,'', 1)
  line = line.replace(' ','')
  L = line.split(',')
  pairs = []
  for ass in L:
 
    sides = ass.split('=')
    if len(sides) == 1:
      # there is = sign in string sides. It must be concatenated with the previous one
      if len(pairs) > 0 and len(pairs[-1]) == 2:

        if isinstance(pairs[-1][1], float):
          pairs[-1][1] = [ pairs[-1][1], float(sides[0]) ]
        elif isinstance(pairs[-1][1], list):
          if isinstance( pairs[-1][1][-1], float  ):
            pairs[-1][1].append( float(sides[0]) )
          elif isinstance( pairs[-1][1][-1], str  ) or isinstance( pairs[-1][1][-1], unicode  ):
            pairs[-1][1].append( sides[0] )     

          else:
            raise Exception('Unknown datatype')
        else:
            raise Exception('Error')

    else:

      try:
        value = float(sides[1])
      except:
        value = sides[1]

      pairs.append( [ sides[0] ,  value ]  )

  pairs = [ (e[0] , e[1] ) for e in pairs  ]

  return (block_name,pairs)

def parse_control_file(fname,rem_chars = ['\n','/']):

  """
  Read and parse a UVic/ MOM2 - style control file to yield list(s) of (name, value) pairs.

  Args:
       fname: (str) full path to configuration file	
       rem_chars: (list of str) characters to ignore in file (defaults recommended)

  Returns: 
       La, L. Here, La is list of (name, value) pairs of all parameters found in the control file, excluding names of parameter groups (generally preceded with & in file). If the same parameter name occurs more than once in the file, it is overwritten in the list (unique names are expected). L is a more detailed list. List of (parameter group name, list of (name, value) ) pairs, where the lists of (name, value) pairs belong to each parameter group name.

  Data is organized in named blocks. The names of these blocks are included in the second return value L, but not L_all. L_all contains the raw (name, value) pairs
  """

  data = read(fname) 

  if data:
    # Remove the chars to be ignored:
    for ec in rem_chars:
      data = data.replace(ec,'')
   
    # identify lines not by \n but by &, indicating new block
    lines = data.split('&')

    # obtain a list of (block name, list of (name,value) pair) pairs. 
    L = []
    for line in lines:
      pair = _parse_control_line(line)
      if pair:
        L.append(pair)

    # obtain the first return value L_all by concatenating the lists of all (name, value) pairs for the list of (block name, list) pairs: 
    L_all = reduce(lambda x,y: x+y, [it[1] for it in L] )  

    return L_all, L


def simple_glob(L, name_filter):

  """
  Very simple wildcard expansion on list of strings.
 
  Args:
       L: list of strings	
       name_filter: the wildcard glob pattern (e.g. '*hello')

  Returns: 
       A list of strings that match.

  Applies very simple shell wildcard \* expansion-type matching of argument name_filter on a list of strings, argument L, without using re module.
  \* can only appear at beginning or end of name_filter (e.g. '\*foo\*').


  E.g. name_filter = '\*oo' and L=['foo','bar'] yields ['foo']
    """

  if name_filter is None:
    return L
  else:
    if name_filter[0] == '*' and name_filter[-1] != '*':
      L = [ it for it in L if it.endswith(name_filter[1:])  ] 
    elif name_filter[0] != '*' and name_filter[-1] == '*':
      L = [ it for it in L if it.startswith(name_filter[:-1])  ] 
    elif name_filter[0] == '*' and name_filter[-1] == '*':
      L = [ it for it in L if name_filter[1:-1] in it  ] 
    else:
      L = [ it for it in L if name_filter == it  ] 
  return L


def end_of_filepath(path):
  """
  Little function that finds the last element of a path, without using os module.
  e.g. /foo/bar yields bar, as does /foo/bar/

  """
# only works on unix-like systems. 
  L_tmp = path.split('/')
  L_tmp = [it for it in L_tmp if it !='']
  return L_tmp[-1]


def affix(coord_name ,affix = '', kind = 'suffix'):
  """
  Append affix to string, ignore when affix already present.

  Can replace this with a nillpotent decorator.
  """

  if affix in coord_name:
    return coord_name

  if kind == 'suffix':
    return coord_name + affix
  elif kind == 'prefix':
    return affix + coord_name 
  else:
    raise Exception('Provide suffix or prefix for kind.')


# ---------------- fill related functions ---------------------------

def _eval_node(A,node):

  return A[node[0],node[1]]

def _set_node(A,node, value = 1):

  A[node[0],node[1]] = value

def _test_node(A, node, test_value = 0, fill_value = 1):

  en = _eval_node(A,node)

  return ( en != test_value) and ( en != fill_value)

def _test_node_append(A, node,Q, test_value = 0, fill_value = 1):

  if _test_node(A,node, test_value, fill_value):
    Q.append(node)
    


def move_west(node):

  return (node[0],node[1]-1)


def move_east_maker(cyclic = False, m = 10000):

  if cyclic:

    def move_east(node):
      return ( node[0], (node[1] + 1)%m  )

  else:

    def move_east(node):
      return ( node[0], node[1] + 1  )

  return move_east


def move_north_maker(cyclic = False, m = 10000):

    if cyclic:

      def move_north(node):
        return ( (node[0] + 1)%m,node[1] )
   
    else:

      def move_north(node):
        return (node[0] + 1,node[1] )

    return move_north

def move_south(node):

  return (node[0] - 1,node[1] )


def _embed_param(shape,x_cyclic,y_cyclic):

  if x_cyclic:
    slicex = slice(None)  
  else:
    shape[1] += 2 
    slicex = slice(1,-1,None)  

  if y_cyclic:
    slicey = slice(None)  
  else:
    shape[0] += 2 
    slicey = slice(1,-1,None)  

  return shape, slicex, slicey

def _embed(mask,x_cyclic = False, y_cyclic = False):

  shape = list(mask.shape)

  shape, slicex, slicey = _embed_param(shape,x_cyclic,y_cyclic)

  new_mask = np.zeros(shape,np.byte)

  new_mask[slicey,slicex] = mask

  return new_mask

def _de_embed(mask,x_cyclic = False, y_cyclic = False):

  shape = list(mask.shape)
  shape, slicex, slicey = _embed_param(shape,x_cyclic,y_cyclic)

  return mask[slicey,slicex]

def _slice_mask(mask,xmin=0,xmax=10000,ymin=0,ymax=10000):
  """Pick a rectangular subsection of the mask.
  """

  slicex_min = slice(None,xmin,None)
  slicex_max = slice(xmax,None,None)

  slicey_min = slice(None,ymin,None)
  slicey_max = slice(ymax,None,None)

  mask[:,slicex_min] = 0
  mask[:,slicex_max] = 0      

  mask[slicey_min,:] = 0
  mask[slicey_max,:] = 0      


def floodfill(A, node = (0,0),boundary_value = np.nan, xmin=0,xmax=10000,ymin=0,ymax=10000,x_cyclic = False, y_cyclic = False, mask_val = 2):
  """
  Fill array (e.g. ocean) from node up to boundary defined by boundary_value (e.g. land) using floodfill.

  Creates mask to fill array (e.g. ocean) in area contained within boundary_value (e.g. land), containing node.

  Args:
    A: (ndarray) 2 dimensional array to use fill on
    node: (2 tuple of int) coordinates of starting point: (y,x) values.
    boundary_value: (float or nan) value of the boundary of the filled domain
    xmin: (int) set mask to boundary value up to this x-index
    xmax: (int) set mask to boundary value from this x-index
    ymin: (int) set mask to boundary value up to this y-index
    ymax: (int) set mask to boundary value from this y-index
    x_cyclic: (Boolean) indicates whether x-domain cyclic if True
    y_cyclic: (Boolean) indicates whether y-domain cyclic if True
    mask_val: (int) value of the masked out region: nodes outside fill in mask

  Returns:
    2 dimension ndarray int mask of filled values
  """

  shape = A.shape

  xmin = max(0,xmin)
  ymin = max(0,ymin)
  xmax = min(shape[1],xmax)
  ymax = min(shape[0],ymax)

  mask = mask_val*np.ones(A.shape,np.byte)

  if np.isnan(boundary_value):
    mask[np.isnan(A)] = 0
  else:
    mask[A == boundary_value] = 0

  _slice_mask(mask,xmin,xmax,ymin,ymax)
   
  mask = _embed(mask,x_cyclic, y_cyclic )

  # Create these functions depending on cyclicity (re-entrance) for speed:
  move_east = move_east_maker(x_cyclic,shape[1])
  move_north = move_north_maker(y_cyclic,shape[0])

  Q = []

  if not _test_node(mask,node):

    warnings.warn('node selected inside boundary value or fill value area. Returning non-filled mask')
    return _de_embed(mask,x_cyclic = x_cyclic, y_cyclic = y_cyclic)

  Q.append(node)

  while Q:
    N = Q.pop(0)

    if _test_node(mask,N):
      w = N
      e = move_east(N)

      while _test_node(mask,w):
        _set_node(mask,w)
        n = move_north(w)
        _test_node_append(mask, n,Q)          
        s = move_south(w)
        _test_node_append(mask, s,Q)      
        w = move_west(w)        
        
      while _test_node(mask,e):
        _set_node(mask,e)
        n = move_north(e)
        _test_node_append(mask, n,Q)          
        s = move_south(e)
        _test_node_append(mask, s,Q)      
        e = move_east(e)        

  return _de_embed(mask,x_cyclic = x_cyclic, y_cyclic = y_cyclic)


# ------------- general time series related functions ----------------


def _str_html_row(row, cell_tag='td'):
  """
  Helper function to construct a HTML table row from a list of string values.
  row: list of row values
  cell_tag = the HTML tag for a cell. Usually 'td' or 'th'

  Used by _add_html_row_str
  """

  open_tag = '<' + cell_tag + '>'
  close_tag = '</' + cell_tag + '>'


  row_str = open_tag
  for r in row[:-1]:
    row_str = row_str + str(r) + close_tag + open_tag

  row_str = row_str + str(row[-1]) + close_tag

  return row_str



def _float_html_row(row, cell_tag='td', decimal = "%.3f"):
  """
  Helper function to construct a HTML table row from a list of string values.
  row: list of row values
  cell_tag = the HTML tag for a cell. Usually 'td' or 'th'

  """

  open_tag = '<' + cell_tag + '>'
  close_tag = '</' + cell_tag + '>'


  row_str = open_tag
  for r in row[:-1]:
    row_str = row_str + decimal%r + close_tag + open_tag

  row_str = row_str + decimal%row[-1] + close_tag

  return row_str

def _add_html_row(row, t_str = '', row_tag = 'tr', cell_tag = 'td' , decimal = "%.3f"):

  open_tag = '<' + row_tag + '>'
  close_tag = '</' + row_tag + '>'

  row_str = _float_html_row(row, cell_tag = cell_tag, decimal = decimal)

  return t_str + open_tag +'\n' + row_str + '\n' + close_tag   

def _add_html_row_str(row, t_str = '', row_tag = 'tr', cell_tag = 'td' ):

  open_tag = '<' + row_tag + '>'
  close_tag = '</' + row_tag + '>'

  row_str = _str_html_row(row, cell_tag = cell_tag)

  return t_str + open_tag +'\n' + row_str + '\n' + close_tag   



def _add_html_rows(T, t_str = '', row_tag = 'tr', cell_tag = 'td', format = "%.3f"):

  """
  Convert block of rows to HTML. Does not yet add table tags.
  format argument determines whether data is string (value "string" passed) or float (formatting string passed). 
  """


  if format == "string":

    for row in T:
      t_str = _add_html_row_str(row = row, t_str = t_str , row_tag = row_tag, cell_tag = cell_tag)

  else:
    # should check format string here with a regexp to avoid strange errors downstream.

    for row in T:
      t_str = _add_html_row(row = row, t_str = t_str , row_tag = row_tag, cell_tag = cell_tag, decimal = format)

  return t_str

def _tryfloat(value, nanval):

  try:
    ret_val = float(value)
  except:
    ret_val = nanval

  return ret_val

def tohtml(T, H = None, D = None,topleft = '', table_tag = 'table', format = "%.3f"):

  """
  Converts lists containing data, headers and dates to a HTML table.

 
  Args:
       T: data
       H: header
       D: dates
       topleft: string to fill topleft cell of table (e.g. a name)
       table_tag: HTML tag to indicate table
       format: string indicating the text formatting (see default).

  Returns: 
       A string containing the HTML table only.
  """


  if H is not None:
    if D is not None:

      T_str  = _add_html_row_str(  row = [topleft,] + H  , cell_tag = 'th' )
      new_T = []
      for line, row in enumerate(T):
         new_T.append( [D[line], ] + row  )
      T_str = _add_html_rows(t_str = T_str, T = new_T, format = format)

    else:
      T_str  = _add_html_row_str( row =  H  , cell_tag = 'th' )
  
      T_str = _add_html_rows(t_str = T_str, T = T, format = format)

  else:

    if D is not None:

      
      new_T = []
      for line, row in enumerate(T):
         new_T.append( [D[line], ] + row  )
      T_str = _add_html_rows(t_str = '', T = new_T)

    else:
      T_str = _add_html_rows(t_str = '', T = T)

  table_open_tag = '<' + table_tag + '>\n'
  table_close_tag = '</' + table_tag + '>'

  return table_open_tag + T_str + table_close_tag
    

def fromcsv(raw_string,row_header=0,col_dates=0,data_start_row =1,data_start_col = 1,hor_desc_col=None, separator =',',row_delin = '\n', nanval = -999999):

  """
  Convert csv string to list of lists of floats T, list of header strings H and (optional) list pf dates D.

  """

  if hor_desc_col is None:

    hor_desc_col = col_dates

  if data_start_col is None:
    data_start_col = col_dates + 1 


  ROWS = raw_string.split(row_delin)

  H = [h.strip('\r') for h in ROWS[row_header].split(separator)[data_start_col:] ]
  

  D = []
  T = []
  for row in ROWS[data_start_row:]:
    cols = row.split(separator)
    row_list = [_tryfloat(e,nanval) for e in cols[data_start_col:]]

    if row_list:
      T.append(row_list)
    dte_str = str(cols[col_dates]  )     
    if dte_str not in  ('\r', ''):     
      D.append(dte_str)

  return T,H,D


def tocsv(T,H=None,D=None, separator =','):  

  """ 
  Converts lists to csv string.
  T list of row lists containing cols of floats.
  optional:
  H list of head strings
  D list of dates. Must be of length len(T)

  """

  if H is None: 
    T_str = ''  
  else:
    T_str = str(H[0])
    for hcol in H[1:]:
      T_str = T_str + separator + str(hcol)
    T_str = T_str + '\n'

  if D is None:

    for row in T:
      row_str =  str(row[0])
      for col in row[1:]:
        row_str = row_str +separator + str(col)
      T_str = T_str + row_str + '\n'   

    return T_str

  else:
    
    for line, row in enumerate(T):
      row_str = str(D[line]) + separator + str(row[0])
      for col in row[1:]:
        row_str = row_str + separator + str(col)
      T_str = T_str + row_str + '\n'   

    return T_str



