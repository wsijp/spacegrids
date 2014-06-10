#encoding:utf-8

""" io related
"""

from config import *

from fieldcls import *

class Operator():

  def __call__(self,F):
    return None # subclasses will override this

  def __pow__(self,n):

    return reduce(lambda x,y: x*y, n*[self])

  def __neg__(self):
    return FloatMul(-1.)*self

  def __sub__(self,other):
    return self + (-other)

  def __add__(self,other):
    return Add(self,other)

  def __mul__(self,other):

    if isinstance(other,Operator):
      return Mul(self,other)

    elif isinstance(other,field) | isinstance(other,vfield):
      return self(other)

class FloatMul(Operator):

  def __init__(self, value):
    self.value = value

  def __call__(self,F):
    return F*self.value

class Add(Operator):

  def __init__(self,left,right):
    self.left = left
    self.right = right

  def __call__(self,F):

    right_result = self.right(F)
    if right_result:
      left_result = self.left(F)
      if left_result:
        combined_value = left_result + right_result
        return combined_value
    return None

class Mul(Operator):

  def __init__(self,left,right):
    self.left = left
    self.right = right

  def __call__(self,F):

    right_result = self.right(F)
    if right_result:
      return self.left(right_result)
   
       
    return None

def make_lin_form(Xi):
  
    return lambda :lambda x : Xi(x)

class Lazy(Operator):

  def __init__(self,op_fctn):
    self.operator = None
    self.fctn = op_fctn 

  def __call__(self, F):

    if not(self.operator):
      self.operator = self.fctn()
    
    return self.operator(F)

class Pick(Operator):

  def __init__(self,Xi):

    self.ax = Xi
     
  def __call__(self,UU):

    if isinstance(UU,field):
      if UU.direction == self.ax:
        return UU
      else:
        return None

    elif isinstance(UU,vfield):
      for U in UU:
        if U.direction == self.ax:
          return U
      return None
    else:
      raise Exception('Error in pick operator %s on %s, argument must be field or sfield. ' % (self.ax, U.name))

class Set_Direction(Operator):
  
  def __init__(self,Xi):
    self.ax = Xi
    return 

  def __call__(self,F):
    
     return F.copy(direction = self.ax)


class Innprod(Operator):
  """
  Inner product on vfields, yields a field. Both fields must be in same linear space: their direction must match up to a permutation. Can be extended later to include non-trivial inner products.
  """
  
  def __call__(self,left,right):

    if isinstance(left,field) and isinstance(right,field):
      return left*right

    elif len(left) == len(right):

      if find_perm(left.direction(),right.direction()):
        return (left*right).sum()
      else:
        raise Exception('Error in inner product %s * %s, must be in same space. ' % (left,right))

    else:
      raise Exception('Error in inner product %s * %s, must be of equal length ' % (left,right))

class Der(Operator):

  def __init__(self,Xi):
    self.ax = Xi

  def __call__(self,vF): 
    if isinstance(vF,field):
      return self.ax^vF
    elif isinstance(vF,vfield):
      return vfield([ self.ax^e for e in vF ])
    else:
      raise Exception('Error in %s - derivative of %s, argument must be field or vfield ' % (self.ax,vF))


class Integ(Operator):

  def __init__(self, Xi):

    self.ax = Xi

  def __call__(self,vF):

    if isinstance(vF,field):
      return self.ax*vF
    elif isinstance(vF,vfield):
      return vfield([self.ax*e for e in vF])
    else:
      raise Exception('Error in %s - integration of %s, argument must be field or vfield ' % (self.ax,vF))

class Mean(Operator):

  def __init__(self, Xi):
    # Argument can be a grid too
    self.ax = Xi

  def __call__(self,vF):


    if isinstance(vF,field):
      return vF/self.ax
    elif isinstance(vF,vfield):
      return vfield([self*e for e in vF])
    else:
      raise Exception('Error in taking %s - primitive of %s, argument must be field or vfield ' % (self.ax,vF))    


  def __mul__(self,other):
     
    if isinstance(other,Mean):
      return Mean(self.ax*other.ax)
    else:
      return self(other)

class Prim(Operator):

  def __init__(self, Xi):

    self.ax = Xi

  def __call__(self,vF):

    if isinstance(vF,field):
      return self.ax|vF
    elif isinstance(vF,vfield):
      return vfield([self.ax|e for e in vF])
    else:
      raise Exception('Error in taking %s - primitive of %s, argument must be field or vfield ' % (self.ax,vF))


class Field_Mul(Operator):

  def __init__(self, F):

    self.F = F

  def __call__(self,vF):

    return self.F*vF

class Sum(Operator):

  def __call__(self, vF):

    if isinstance(vF,field):
      return vF
    elif isinstance(vF,vfield):
      result = reduce(lambda x,y: x+y, vF)
      result.direction = ID()
      return result

    else:

      raise Exception('Error in Sum %s, argument must be field or vfield ' % vF)

class Slice(Operator):

  def __init__(self,tup):

    # tup argument contains a tuple containing coord-slice obj pairs.

    self.tup = tup

  def __call__(self, vF):

   if isinstance(vF,field):
     return vF[self.tup]
   elif isinstance(vF,vfield):
     return vfield([self*e for e in vF])
   else:

      raise Exception('Error in Slice %s, argument must be field or vfield ' % vF)

class Identity(Operator):

  def __call__(self,vF):

    return vF

nop = Identity();

class If(Operator):

  def __init__(self, Cond_Op, True_Op, Else_Op = nop):

    # The state of the system is defined by the field(s) vF. Therefore, the condition can be internal to the call function and must be defined in terms of vF (e.g. 'len(vF.gr) == 3'). 
    self.Cond_Op = Cond_Op
    self.True_Op = True_Op
    self.Else_Op = Else_Op

  def __call__(self,vF):

   if isinstance(vF,field) |  isinstance(vF,vfield):
     
     result = self.Cond_Op(vF)

     if result != vF:
       return self.True_Op(vF)
     else:
       return self.Else_Op(vF)

   else:

      raise Exception('Error in If %s, argument must be field or vfield ' % vF)


class While(Operator):

  def __init__(self, True_Op,Cond_Op):

    self.Cond_Op = Cond_Op
    self.IF = If(Cond_Op = Cond_Op,True_Op = True_Op)
  
  def __call__(self,vF):

   if isinstance(vF,field) |  isinstance(vF,vfield):
     
     prev_result = vF
     result = self.IF*prev_result

     while result != prev_result:
       prev_result = result
       result = self.IF*result

     return result       

   else:

      raise Exception('Error in While %s, argument must be field or vfield ' % vF)

class Try(Operator):

  def __init__(self, Op):

    self.Op = Op

  def __call__(self,vF):

    try:
      result = self.Op*vF
    except:
      result = vF

    return result


dot = Innprod()
