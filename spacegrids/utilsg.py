#encoding:utf-8

""" useful general functions
"""

import numpy as np
import types
import math

import fnmatch
import itertools

from config import *

class Report():
  """
  Reporting class containing useful method for building a string containing messages.
  """

  def __init__(self,value = ''):

    self.value = value

  def __repr__(self):
    return self.value

  def line(self,char = '-',times = 10):
    if not self.value[-1] == '\n':
      self.echoln()

    self.echoln(char*times)

  def block(self,L,delim = '\ ', width = 12, cols = 4 ):

    Ls = split_list(L , size = len(L)/cols + 1)
    Ls = transpose_list(Ls)
 

    for LL in Ls:
      self.echo(what = delim.join([ ("%-{}s".format(str(width)))%it for it in  LL ] ) )
      self.echoln()

  def echoln(self,what = '', delim = ' ', maxlen = 10, width = 12, cols = 4):

    self.echo(what = what , delim = delim, maxlen = maxlen, width = width,cols = cols)
    self.echo('\n') 


  def echo(self,what='' , delim = ' ', maxlen = 10, width = 12, cols = 4):
    """
    Doc here
    """


    if isinstance(what,str) or isinstance(what,unicode):
      self.value = self.value + what
    elif isinstance(what,list) or isinstance(what,tuple):
      if len(what) > maxlen:

        self.block(L = what[:maxlen/2] + ['...',] +  what[-maxlen/2:] , delim = delim, width = width, cols = cols)

      else:
        self.block(L = what, delim = delim, width = width, cols = cols)
     
  def table(self,D,cols = 4, numspace =2):
    """ Doc here

    """
    for i, k in enumerate(D.keys()):
      if (2*i%cols == 0):
        self.echoln()

      self.echo('%-15s %-15s' % (k , D[k]))

    

  def __add__(self,other):
    return Report(value = self.value + other.value)

 
def print_box(L, cols = 5, numspace = 2):
  """
  Used for formatted output to stdout
  """

  for i, it in enumerate(L):
    if (i%cols == 0):
      print '\n', 
    print '%-15s' % (it) ,

def print_table(D, cols = 4, numspace = 2):
  """
  Used for formatted output to stdout
  """

  for i, k in enumerate(D.keys()):
    if (2*i%cols == 0):
      print '\n', 
    print '%-15s %-15s' % (k , D[k]),



def split_list(L, size = 10):
  """
  Very general function splitting a list L into a list of sublists of equal length (and a remainder).
  """

  length = len(L)
  full_blocks = length/size

  new_L = [ L[i*size:(i+1)*size]   for i in range(full_blocks) ]
  
  if full_blocks * size != length: 
    new_L.append(L[size*full_blocks:])  

  return new_L

def transpose_list(L):
  """
  Transpose a list of lists.

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
  Merge two ndarrays or iterables and order them. 

  Args:
       A1: first array
       A2: second array to be combined with first array
       sort: flag. sort called on result if True.

  Returns:
       new_list: sorted list of combined values.

  Used to merge Coord values.
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


def mark_sublist(Lbig,Lsub, indicator ='*'):
  """
  Marks all strings in a list of strings that also occur in a sublist by appending an indicator (\*) at the end of those strings.

  Args:
       Lbig:	The larger list, of which the contained strings will be marked.
       Lsub:	The list of strings that need to be marked
       indicator:	string containing the marker.
       
  Returns: 
       cLbig -- the list with marked strings
  """

  cLbig = []
  for it in Lbig:
    if it in Lsub:
      it =   it +indicator
  
    cLbig.append(it)
  
  return cLbig    


def add_alias(L):
  """Create alias attribute based on name attribute, and identify objects that have identical name attributes to yield corresponding numbered alias attributes to tell them apart.

  Args:
       L:	list of Coord objects (or any objects with a name attribute).
  
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

  **Examples**
 
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


def complete(comp_sets):
  """
  Expects a list or tuple of tuples (or lists).
  E.g. 
  ((xt,yt,zt),(xu,yu,zw))
  It then assignes an attribute others to every element xt,yt,...  that is a tuple containing all other elements in the subset except the element itself.

  i.e. xt.others = (yt,zt)
  

  """
  for set in comp_sets:
    
    for co in set:
      co.others = tuple([e for e in set if e != co])

  return

def flat_list(L):
  W = list(itertools.chain(*L))
  while isinstance(W[0],list):
    W = list(itertools.chain(*W))

  return W


def rav_index(L,sh):
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
  Inputs: a,b tuples or lists of equal length containing the same elements, possibly permuted.

  Output: a tuple of indices perm such that [left[i] for i in perm] = right

  """

  if len(left) == len(right):
    perm = []
   
    for r in right:
      if r in left:
        perm.append(left.index(r))
      else:
        if verbose:
          print 'inputs not permutable: ' + r.name
        return

  else:
    if verbose:
      print "message from find_perm: inputs must be of equal length."
    return 
  return tuple(perm)



def round_order(value, order = None):

  if order is None:
    order = math.log10(abs(value))//1.
  
  return (10**order)*round(value*10**-order)


def order_mag(val):
  
  if val == 0.:

    return 0.
  else:
    return math.log10(abs(val))//1.





def dist_angle(angles):

   new_angles = []
   for a in angles:
     if a> 180.:
       new_angles.append(360. - a)
     else:
       new_angles.append(a)
  
   return np.array(new_angles)

def which_att(obj,att_list,fail_val = None):

  for attname in att_list:

    if hasattr(obj,attname):
      return attname
  
  return fail_val


def get_att(obj, att_list,fail_val = None):

  this_att = which_att(obj, att_list)
  if this_att is None:
    return fail_val
  else:
    return getattr(obj, this_att)



def id_index(L,e):

  for i,it in enumerate(L):
    if (it&e):
      return i

# not raising error as normal in would. 
#  raise ValueError('%s is not in list' % e)

# return None instead
  return

def id_in(L,e):

  for it in L:
    if (it&e):
      return True
  return False

def rem_equivs(L):

  L_new = []

  for e in L:
    if not(id_in(L_new, e)):
      L_new.append(e)

  return L_new



def read(fname):

  data = None

  with open(fname,'r') as f:

    data = f.read();

  return data


def interp_line(line):

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
  Read and parse a UVic/ MOM2 - style control file. Returns list(s) of (name, value) pairs found in the control file.

  Input: fname full path to configuration file	
  Output: La, L
  La	list of (name, value) pairs of all parameters found in the control file, excluding names of parameter groups (generally preceded with & in file). If the same parameter name occurs more than once in the file, it is overwritten in the list (unique names are expected).
  L	more detailed list. List of (parameter group name, list of (name, value) ) pairs, where the lists of (name, value) pairs belong to each parameter group name.


  """

  data = read(fname) 

  if data:
  
    for ec in rem_chars:
      data = data.replace(ec,'')
   
    lines = data.split('&')

    L = []
    for line in lines:
      pair = interp_line(line)
      if pair:
        L.append(pair)

    L_all = reduce(lambda x,y: x+y, [it[1] for it in L] )  

    return L_all, L


def simple_glob(L, name_filter):

    """
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
  utility that finds the last element of a path without using os module.
  e.g. /foo/bar yields bar, as does /foo/bar/

  """
# only works on unix-like systems. 
  L_tmp = path.split('/')
  L_tmp = [it for it in L_tmp if it !='']
  return L_tmp[-1]





# ------------- general time series related functions ----------------


def str_html_row(row, cell_tag='td'):
  """
  Helper function to construct a HTML table row from a list of string values.
  row: list of row values
  cell_tag = the HTML tag for a cell. Usually 'td' or 'th'

  """

  open_tag = '<' + cell_tag + '>'
  close_tag = '</' + cell_tag + '>'


  row_str = open_tag
  for r in row[:-1]:
    row_str = row_str + str(r) + close_tag + open_tag

  row_str = row_str + str(row[-1]) + close_tag

  return row_str



def float_html_row(row, cell_tag='td', decimal = "%.3f"):
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

def add_html_row(row, t_str = '', row_tag = 'tr', cell_tag = 'td' , decimal = "%.3f"):

  open_tag = '<' + row_tag + '>'
  close_tag = '</' + row_tag + '>'

  row_str = float_html_row(row, cell_tag = cell_tag, decimal = decimal)

  return t_str + open_tag +'\n' + row_str + '\n' + close_tag   

def add_html_row_str(row, t_str = '', row_tag = 'tr', cell_tag = 'td' ):

  open_tag = '<' + row_tag + '>'
  close_tag = '</' + row_tag + '>'

  row_str = str_html_row(row, cell_tag = cell_tag)

  return t_str + open_tag +'\n' + row_str + '\n' + close_tag   



def add_html_rows(T, t_str = '', row_tag = 'tr', cell_tag = 'td', format = "%.3f"):

  """
  Convert block of rows to HTML. Does not yet add table tags.
  format argument determines whether data is string (value "string" passed) or float (formatting string passed). 
  """


  if format == "string":

    for row in T:
      t_str = add_html_row_str(row = row, t_str = t_str , row_tag = row_tag, cell_tag = cell_tag)

  else:
    # should check format string here with a regexp to avoid strange errors downstream.

    for row in T:
      t_str = add_html_row(row = row, t_str = t_str , row_tag = row_tag, cell_tag = cell_tag, decimal = format)

  return t_str

def tryfloat(value, nanval):

  try:
    ret_val = float(value)
  except:
    ret_val = nanval

  return ret_val

def tohtml(T, H = None, D = None,topleft = '', table_tag = 'table', format = "%.3f"):

  """
  T data
  H header
  D dates

  """


  if H is not None:
    if D is not None:

      T_str  = add_html_row_str(  row = [topleft,] + H  , cell_tag = 'th' )
      new_T = []
      for line, row in enumerate(T):
         new_T.append( [D[line], ] + row  )
      T_str = add_html_rows(t_str = T_str, T = new_T, format = format)

    else:
      T_str  = add_html_row_str( row =  H  , cell_tag = 'th' )
  
      T_str = add_html_rows(t_str = T_str, T = T, format = format)

  else:

    if D is not None:

      
      new_T = []
      for line, row in enumerate(T):
         new_T.append( [D[line], ] + row  )
      T_str = add_html_rows(t_str = '', T = new_T)

    else:
      T_str = add_html_rows(t_str = '', T = T)

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
    row_list = [tryfloat(e,nanval) for e in cols[data_start_col:]]

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



