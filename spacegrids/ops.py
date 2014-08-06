#encoding:utf-8

""" Mathematical operators based on Field and Gr operations.

The mathematical operators defined here, such as Der, provide basic building blocks for the user to define more detailed Operators. Arithmatic operations such as addition and multiplication among these operators yield new operators. These new operators can be called on Field objects later, performing the intended operation. 

Example:

  >>> ddX = sg.Der(X) # create Operator that would take X-derivative on Field args
  >>> ddY = sg.Der(Y) # similar
  >>> ddXddY = ddX*ddY  # create composite operator
  >>> ddXddY(F) # differentiation in X direction following differentiation in Y
"""

from _config import *

from fieldcls import *

class Operator(object):
  """
  Operator base class that is called on Field objects. 

  Has addition and multiplication-related operations, and empty ___call___. 
  """

  def __call__(self,F):
    """
    To be overriden in subclasses.
    """
    return None 

  def __pow__(self,n):
    """Power: repeated multiplication.
    """
    return reduce(lambda x,y: x*y, n*[self])

  def __neg__(self):
    """Negative: multiplication by -1.
    """
    return FloatMul(-1.)*self

  def __sub__(self,other):
    return self + (-other)

  def __add__(self,other):
    return Add(self,other)

  def __mul__(self,other):

    if isinstance(other,Operator):
      return Mul(self,other)

    elif isinstance(other,Field) | isinstance(other,VField):
      return self(other)

class FloatMul(Operator):
  """
  Operator that is initalized with a float, and multiplies a Field argument with that float on calling.
  """

  def __init__(self, value):
    """
    Initialize float multiplication by assigning value attribute.
    """
    self.value = value

  def __call__(self,F):
    """
    Args:
      F: (Field) to act on.
    
    Returns:
      Field that is arg Field times value attribute, a float.
    """
    return F*self.value

class Add(Operator):
  """
  Operator for addition among Operator objects. 

  Stores left and right Operator __init__ arguments, that will both act on Field __call__ argument and then added there.
  """
  def __init__(self,left,right):
    """
    Takes "left" and "right" Operator arguments, and stores them in attributes.
    """
    self.left = left
    self.right = right

  def __call__(self,F):
    """
    Args:
      F: (Field) to act on.
    
    Returns:
      Field: sum of "left" and "right" operator called on arg Field, or None if something went wrong.
    """

    right_result = self.right(F)
    if right_result:
      left_result = self.left(F)
      if left_result:
        combined_value = left_result + right_result
        return combined_value
    return None

class Mul(Operator):
  """
  Operator for multiplication among Operator objects.


  Stores left and right Operator __init__ arguments, that will both act on Field __call__ argument and then multiplied there.
  """


  def __init__(self,left,right):
    """
    Takes "left" and "right" Operator arguments, and stores them in attributes.
    """

    self.left = left
    self.right = right

  def __call__(self,F):
    """
    Args:
      F: (Field) to act on.
    
    Returns:
      Field: product of "left" and "right" operator called on arg Field, or None if something went wrong.
    """

    right_result = self.right(F)
    if right_result:
      return self.left(right_result)
   
       
    return None

def make_lin_form(Xi):
    """ pass
    """  
    return lambda :lambda x : Xi(x)

class Lazy(Operator):
  """
  Operator created with a function argument. The function yields an operator, and is called only in __call__, where an operator attribute is assigned to self with that return value. From then on, that operator will act on the __call__ Field argument.
  """

  def __init__(self,op_fctn):
    self.operator = None
    self.fctn = op_fctn 

  def __call__(self, F):

    if not(self.operator):
      self.operator = self.fctn()
    
    return self.operator(F)

class Pick(Operator):
  """Constructed with argument Xi, picks the VField component with direction matching that Xi and returns it.
  """
  def __init__(self,Xi):
    """
    Initialize Pick Operator.

    Args:
      Xi: (str) containing direction (e.g. 'X') compatible with direction attributes of VField members.

    Raises:
      TypeError if argument not Field or VField.
    """
    self.ax = Xi
     
  def __call__(self,UU):
    """
    Picks matching VField.direction to self.ax in VField.

    Projection of vector field to an axis.

    Args:
      UU: (VField) to pick member from

    Returns:
      Field or None depending on success.
    """

    if isinstance(UU,Field):
      if UU.direction == self.ax:
        return UU
      else:
        return None

    elif isinstance(UU,VField):
      for U in UU:
        if U.direction == self.ax:
          return U
      return None
    else:
      raise TypeError('Error in pick operator %s on %s, argument must be Field or VField. ' % (self.ax, U.name))

class SetDirection(Operator):
  """Operator class to set the direction on a Field to self.ax.
  """  

  def __init__(self,Xi):
    """Initialize SetDirection class by assigning ax attribute (ax).

    Args:
      ax: (Ax) axis that calling this class will asign to Field __call__ argument 
    """
    self.ax = Xi
    return 

  def __call__(self,F):
     """ Assign direction to copy of Field.
     """
     return F.copy(direction = self.ax)


class Innprod(Operator):
  """
  Inner product on two (V)Fields. 

  Both fields must be in same linear space: their direction must match up to a permutation. Can be extended later to include non-trivial inner products.
  """
  
  def __call__(self,left,right):
    """
    Produce inner product of two (V)Fields. 

    The product of two Field objects yields normal field multiplication. 
    Multiplying two VField objects requires multiplicants to be of equal length and have the same direction attributes up to a permutation (e.g. X*Y vs Y*X is ok).

    Args:
      left: (Fiel or VField). left multiplicant
      right: (Fiel or VField). right multiplicant

    Returns: 
      A Field containing the sum of the matching products.

    Raises:
      ValueError if VFields are not in same linear space (directions don't match up to permutation). 
    """


    if isinstance(left,Field) and isinstance(right,Field):
      return left*right

    elif len(left) == len(right):

      if find_perm(left.direction(),right.direction()):
        return (left*right).innersum()
      else:
        raise ValueError('Error in inner product %s * %s, must be in same space. ' % (left,right))

    else:
      raise ValueError('Error in inner product %s * %s, must be of equal length ' % (left,right))

class Der(Operator):
  """Directional derivative Operator.

  Attributes:
    ax: (Ax) axis in which direction __call__ will take the derivative. 

  Example:

    >>> ddX = sg.Der(X) # create Operator that would take X-derivative on Field args
    >>> ddY = sg.Der(Y) # similar
    >>> ddXddY = ddX*ddY  # create composite operator
    >>> ddXddY(F) # differentiation in X direction following differentiation in Y
  """

  def __init__(self,Xi):
    """Initialize Der Operator by assigning Ax object from Xi arg. to ax attribute.
    """
    self.ax = Xi

  def __call__(self,vF): 
    """Call der method of (V)Field argument using self.ax as ax argument.

    This will then call the Ax.der method on VF, which in turn calls the Gr.der method. That will take care of any Coord-dependencies.
    """

    return vF.der(self.ax)
 

class Integ(Operator):
  """Compute total integral in ax direction. 

  Attributes:
    ax: (Ax) axis direction in which to integrate.
   
  See also:
    vsum method.
  """

  def __init__(self, Xi):
    """Initialize the Integ Operator by assigning ax attribute.

    Args:
      ax: (Ax) axis to integrate over.
    """

    self.ax = Xi

  def __call__(self,vF):
    """Initialize the Integ Operator by assigning ax attribute.

    Args:
      vF: (Field or VField) to integrate.

    Returns:
      vsum method called on vF.

    Raises:
      TypeError 
    """

    if isinstance(vF,Field):
      return self.ax.vsum(vF)
    elif isinstance(vF,VField):
      return VField([self.ax.vsum(e) for e in vF])
    else:
      raise TypeError('Error in %s - integration of %s, argument must be Field or VField ' % (self.ax,vF))

class Mean(Operator):
  """Calculate mean of (V)Field object using stored Ax object.
  """
  def __init__(self, Xi):
    """Initialize Mean Operator.
  
    Args:
      Xi: (Ax or Gr) axis or grid along which to take mean when calling __call__
    """
    self.ax = Xi

  def __call__(self,vF):
    """Initialize Mean Operator.
  
    Args:
      vF: (Field or VField) to take mean of 

    Returns:
      Mean of Field.

    Raises:
      TypeError 
    """


    if isinstance(vF,Field):
      return vF/self.ax
    elif isinstance(vF,VField):
      return VField([self*e for e in vF])
    else:
      raise Exception('Error in taking %s - primitive of %s, argument must be Field or VField ' % (self.ax,vF))    


  def __mul__(self,other):
    """Multiplication of Mean object with other Operator.
    """   
    if isinstance(other,Mean):
      return Mean(self.ax*other.ax)
    else:
      return self(other)

class Prim(Operator):
  """Take primitive (antiderivative) of (V)Field along axis.

  Calls vcumsum.
  """
  def __init__(self, Xi):
    """Initialize Prim Operator by specifying ax attribute.

    Args:
      Xi: (Ax) along which to take primitive.
    """
    self.ax = Xi

  def __call__(self,vF):
    """Call vF.vcumsum on self.Xi.
    """
    return vF.vcumsum(self.Xi)


class FieldMul(Operator):
  """Operator storing a Field and multiplying (V)Field argument with that Field when called.
  """
  def __init__(self, F):
    """Initialize FieldMul Operator by storing Field in F attribute.
    """

    self.F = F

  def __call__(self,vF):
    """Multiply (V)Field with self.F and return product (V)Field.
    """
    return self.F*vF

class Sum(Operator):
  """Operator adding the vector components of a VField and returning the resulting Field.
  """
  def __call__(self, vF):
    """
    Args:
      vF: (Field or VField) to take inner sum of.

    Returns:
      Field containing sum.

    Raises:
      TypeError if argument other than (V)Field
    """
    if isinstance(vF,Field):
      return vF
    elif isinstance(vF,VField):
      result = reduce(lambda x,y: x+y, vF)
      result.direction = ID()
      return result

    else:

      raise TypeError('Error in Sum %s, argument must be Field or VField ' % vF)

class Slice(Operator):
  """Slice Operator storing tuple of Coord vs slice object pairs to slice (V)Field arguments with on __call__.

  Attributes:
    tup: (tuple of Coord vs slice objects) to slice Field with
  """

  def __init__(self,tup):
    """ Initialize Slice Operator.
    
    Args:
      tup (tuple of Coord vs slice obj pairs).
    """

    self.tup = tup

  def __call__(self, vF):
   """ Slice (V)Field argument with self.tup attribute.

   Args:
     vF: (Field or VField)
   
   Returns:
     Sliced (V)Field

   Raises:
     TypeError
   """

   if isinstance(vF,Field):
     return vF[self.tup]
   elif isinstance(vF,VField):
     return VField([e[self.tup] for e in vF])
   else:

      raise TypeError('Error in Slice %s, argument must be Field or VField ' % vF)

class Identity(Operator):
  """Operator representing the identity. 

  No change upon __call__
  """
  def __call__(self,vF):
    """Return argument as is.
    """
    return vF

# define an instance of the identity:
nop = Identity();

class If(Operator):
  """If Operator as in If then.

  The state of the system is defined by the Field(s) vF. Therefore, the condition can be internal to the call function and must be defined in terms of vF (e.g. 'len(vF.grid) == 3'). 
  """
  def __init__(self, Cond_Op, True_Op, Else_Op = nop):
    """Initialize If Operator by storing args in attributes. 

    Args:
      Cond_Op: (Operator): to apply to (V)Field argument
      True_Op: (Operator): to apply if this yields a different result
      Else_Op: (Operator): to apply if (V)Field remains unchanged.
    """

    self.Cond_Op = Cond_Op
    self.True_Op = True_Op
    self.Else_Op = Else_Op

  def __call__(self,vF):
   """Perform if-then depending on (V)Field argument.

   If Cond_Op applied to vF yields a different result (from vF), apply True_Op. Apply Else_Op otherwise.

   Raises: TypeError
   """
   if isinstance(vF,Field) |  isinstance(vF,VField):
     
     result = self.Cond_Op(vF)

     if result != vF:
       return self.True_Op(vF)
     else:
       return self.Else_Op(vF)

   else:

      raise TypeError('Error in If %s, argument must be Field or VField ' % vF)


class While(Operator):
  """Similar to If Operator
  """

  def __init__(self, True_Op,Cond_Op):
    """Similar to If Operator
    """

    self.Cond_Op = Cond_Op
    self.IF = If(Cond_Op = Cond_Op,True_Op = True_Op)
  
  def __call__(self,vF):
    """Similar to If Operator
    """

    if isinstance(vF,Field) |  isinstance(vF,VField):
     
      prev_result = vF
      result = self.IF*prev_result

      while result != prev_result:
        prev_result = result
        result = self.IF*result

      return result       

    else:

      raise Exception('Error in While %s, argument must be Field or VField ' % vF)

class Try(Operator):
  """Operator storing an Operator that is called with except on (V)Field argument in __call__ 

  If successful, return result, otherwise return original (V)Field argument.

  Similar to If Operator
  """


  def __init__(self, Op):
    """Initialize Try object by assigning Op Operator argument to attribute Op.
    """
    self.Op = Op

  def __call__(self,vF):
    """Try calling self.Op on vF (V)Field argument and return result.
 
    If this fails, return original vF
    """

    try:
      result = self.Op*vF
    except:
      result = vF

    return result

# define an instance of Innprod:
dot = Innprod()
