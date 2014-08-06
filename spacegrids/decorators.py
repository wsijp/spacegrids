# ---- decorators -----



def check_equiv(func):
  """
  Decorator to perform preliminary equivalence check between self and other.

  Args:
    func: function.

  Returns:
    Function: the decorator.
  """

  def checker(caller,other, *args, **kwargs): 

    
    if not(caller.is_equiv(other)):
      warnings.warn("Ordered equivalence wrong, aborting with None!")
      return

    return func(caller,other, *args, **kwargs)

  return checker




def method2members(func):
  """
  Decorator to transfer method call from Membered object to its members.

  Args:
    func: function.

  Returns:
    Function: the decorator.
  """

  def wrap(caller, *args, **kwargs): 

    method_name = func.__name__
    
    return caller.call_on_members(method_name, *args, **kwargs)

  return wrap


def att2members(func):
  """
  Decorator to transfer attribute retrieval from Membered object to its members.

  Args:
    func: function.

  Returns:
    Function: the decorator.
  """

  def wrap(caller, *args, **kwargs): 

    att_name = func.__name__
    
    return caller.get_from_members(att_name, *args, **kwargs)

  return wrap






