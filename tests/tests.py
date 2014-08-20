#encoding:utf-8

""" unit tests.
Data required for these tests is the example project my_project mentioned in the tutorial.

"""

import inspect
import copy
import unittest
import os
import numpy as np
import spacegrids as sg


class TestValuedClass(unittest.TestCase):


  def setUp(self):
    pass

  def tearDown(self):
    pass


  def test_slice_method(self):

    K = sg.Valued('K',np.array([1.,2.,3.,4.]))
    R=K.sliced(slice(1,None,None))

    self.assertEqual( np.array_equal( R.value, np.array([ 2.,  3.,  4.]) ), True )    


# test the info function

class TestInfo(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__

    self.fixture = sg.info_dict()

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture

  def test_type(self):
    
    self.assertEqual(type(self.fixture),dict)

  def test_type2(self):
    D = self.fixture
    if len(D) > 0:
      self.assertEqual(type(D.keys()[0]),str)

  def test_paths_in_D_exist(self):
    D = self.fixture
    for path in D.values():
      self.assertEqual(os.path.exists(path), True)  

class Test_project_helpers(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info_dict()
    self.fixture = D

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture



  def test_isexpdir_on_project_dir(self):

    D = self.fixture
    self.assertEqual(set(sg.isexpdir(os.path.join(D['my_project']))),  set(['DPO', 'DPC','Lev.cdf'] ) ) 


  def test_isexpdir_on_exper_dir(self):

    D = self.fixture
    self.assertEqual(sg.isexpdir(os.path.join(D['my_project'], 'DPO')),  ['time_mean.nc'] ) 


class TestCoordsOnTheirOwn(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__

    def provide_axis(cstack):
      for i, c in enumerate(cstack):
        cstack[i].axis = cstack[i].direction 

      return cstack
  
    # Note that some coord values are deliberately unordered.   

    # Coords ---    
    coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]), strings = ['one','two','three'] , metadata = {'hi':5} )
    coord2 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]), metadata = {'hi':7})
    coord3 = sg.fieldcls.Coord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':3})
    # identical in main attributes to previous set (in order):
    coord4 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]), metadata = {'hi':8})
    coord5 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1,2,3, 4]), metadata = {'hi':10})
    coord6 = sg.fieldcls.Coord(name = 'test',direction ='X',value =np.array([5,1,2,3, 4]), metadata = {'hi':12})

    # providing coord1 and coord2 with duals. coord3 is self-dual

    coord1_edges = sg.fieldcls.Coord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), strings = ['a','b','c','d'] , dual = coord1 , metadata = {'hi':25} )
    coord2_edges = sg.fieldcls.Coord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = coord2, metadata = {'hi':77})

    # identical in main attributes to previous set (in order):
    coord4_edges = sg.fieldcls.Coord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = coord4 , metadata = {'hi':25} )
    coord5_edges = sg.fieldcls.Coord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = coord5, metadata = {'hi':77})


    # YCoords ---

    ycoord1 = sg.fieldcls.YCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) , metadata = {'hi':5} )
    ycoord2 = sg.fieldcls.YCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]), metadata = {'hi':7})
    ycoord3 = sg.fieldcls.YCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':3})
    # identical in main attributes to previous set (in order):
    ycoord4 = sg.fieldcls.YCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]), metadata = {'hi':8})
    ycoord5 = sg.fieldcls.YCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3., 4.]), metadata = {'hi':10})
    ycoord6 = sg.fieldcls.YCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3., 4.]), metadata = {'hi':12})


    ycoord1_edges = sg.fieldcls.YCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = ycoord1 , metadata = {'hi':25} )
    ycoord2_edges = sg.fieldcls.YCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = ycoord2, metadata = {'hi':77})

    # identical in main attributes to previous set (in order):
    ycoord4_edges = sg.fieldcls.YCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = ycoord4 , metadata = {'hi':25} )
    ycoord5_edges = sg.fieldcls.YCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = ycoord5, metadata = {'hi':77})




    # XCoords ---

    xcoord1 = sg.fieldcls.XCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) , metadata = {'hi':5} )
    xcoord2 = sg.fieldcls.XCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]), metadata = {'hi':7})
    xcoord3 = sg.fieldcls.XCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':3})
    # identical in main attributes to previous set (in order):
    xcoord4 = sg.fieldcls.XCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]), metadata = {'hi':8})
    xcoord5 = sg.fieldcls.XCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3., 4.]), metadata = {'hi':10})
    xcoord6 = sg.fieldcls.XCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3., 4.]), metadata = {'hi':12})

    xcoord1_edges = sg.fieldcls.XCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = xcoord1 , metadata = {'hi':25} )
    xcoord2_edges = sg.fieldcls.XCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = xcoord2, metadata = {'hi':77})

    # identical in main attributes to previous set (in order):
    xcoord4_edges = sg.fieldcls.XCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = xcoord4 , metadata = {'hi':25} )
    xcoord5_edges = sg.fieldcls.XCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = xcoord5, metadata = {'hi':77})




    # we are testing for Coord, YCoord and XCoord 
    cstack1 = provide_axis([coord1,coord2,coord3,coord1_edges,coord2_edges])
    cstack2 = provide_axis([coord4,coord5,coord6,coord4_edges,coord5_edges])

    ycstack1 = provide_axis([ycoord1,ycoord2,ycoord3,ycoord1_edges,ycoord2_edges])
    ycstack2 = provide_axis([ycoord4,ycoord5,ycoord6,ycoord4_edges,ycoord5_edges])

    xcstack1 = provide_axis([xcoord1,xcoord2,xcoord3,xcoord1_edges,xcoord2_edges])
    xcstack2 = provide_axis([xcoord4,xcoord5,xcoord6,xcoord4_edges,xcoord5_edges])


    self.fixture = [cstack1, cstack2, ycstack1, ycstack2,xcstack1, xcstack2,]

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture

  def test_init_method(self):
    """Test the __init__ method of Coord
    """
    
    self.assertRaises(ValueError, sg.fieldcls.Coord, **{'name' : 'test1','direction' :'X','value': np.array([1.,2.,3.]) , 'metadata': {'hi':5}, 'strings': ['foo','bar'] })

  def test_get_item_method(self):
    """Test the __getitem__ method of Coord class for success, failure and raised error.
    """
    coord1 = self.fixture[0][0]

    self.assertEqual(coord1[1], 2.)


  def test_coord_array_equal_method(self):
    """Test the array_equal method of Coord class for success, failure and raised error.
    """
    coord1 = self.fixture[0][0]
    coord2 = self.fixture[0][1]
    coord4 = self.fixture[1][0]

    self.assertEqual(coord1.array_equal(coord2), False)
    self.assertEqual(coord1.array_equal(coord4), True)
    self.assertRaises(TypeError, coord1.array_equal, 5)

  def test_coord_init_attributes_assigned(self):
    """
    Test whether all passed are assigned to attributes as intended. This is easy to forget when adding new arguments.
    """

    pass




  def test_coord_sliced_method(self):
    """Tests whether slicing works"""

    coord1 = self.fixture[0][0] # this one has string property set
    coord2 = self.fixture[0][1] # this one doesn't
    coord4 = self.fixture[1][0]

    K=coord1(coord1*coord2)
    R = coord1.coord_shift(K,1)

    self.assertEqual( (R[1,1:3]).shape, (2,)   )
    self.assertEqual( (isinstance(R[1,1:3]), sg.Field), True   )
    self.assertEqual( (isinstance(R[1,1:2]), sg.Field), False   ) # float
    self.assertEqual( (isinstance(R[1,2]), sg.Field), False   ) # float


  def test_coord_sliced_method(self):
    """Tests whether slicing works"""

    coord1 = self.fixture[0][0] # this one has string property set
    coord2 = self.fixture[0][1] # this one doesn't
    coord4 = self.fixture[1][0]

    self.assertEqual(coord1.sliced(slice(None,None,None) ) is coord1, True   )

    slice_obj = slice(1,None,None)

    coord1_sliced = coord1.sliced( slice_obj )
    coord2_sliced = coord2.sliced( slice_obj )


    self.assertEqual(np.array_equal(coord1_sliced.value,  coord1.value[slice_obj] ) , True)

    self.assertEqual(np.array_equal(coord1_sliced.strings,  coord1.strings[slice_obj] ) , True)


    # the dual Coord should also be sliced and be properly assigned:
    self.assertEqual(len(coord1_sliced.dual.value) , len(coord1.dual.value) -1 )

    # the dual should remain one longer
    self.assertEqual(len(coord1_sliced.dual.value) , len(coord1_sliced.value) + 1 )


    self.assertEqual(coord1_sliced.dual.dual, coord1_sliced)

    self.assertEqual(coord2_sliced.strings,  None)

    # test integer slice:
    self.assertEqual(len(coord1.sliced(2) )  , 1)



    # Now make coord1 self dual to test for self-dual Coord object:

    coord1.give_dual()

    # as an aside, test whether give_dual worked:

    self.assertEqual(coord1.dual is coord1, True)

    # ok, slice again:
    coord1_sliced = coord1.sliced( slice_obj )

    # the sliced coord should remain self-dual:
    self.assertEqual(coord1_sliced.dual is coord1_sliced, True)


  # ---------- test block for Coord class ------



  def test_coord_mult_with_AxGr(self):
    """
    Test copy method of Coord. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
 
    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')
    Z = sg.fieldcls.Ax('Z')


    coord1.give_axis(X)
    coord2.give_axis(Y)
    coord3.give_axis(Z)

    coord3.direction = 'Z'


    self.assertEqual( (X*Y)*(coord1*coord2*coord3) , coord1*coord2  )
    self.assertEqual( (Y*X)*(coord1*coord2*coord3) , coord2*coord1  )





  def test_copy_method_yields_not_same_for_case_name(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is the same as the original (although a different object in memory) and differs in that specific attribute.

    """

    cstack1 = self.fixture[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(name = 'joep')

    self.assertEqual(coord3_copy.name, 'joep'  )

  def test_copy_method_yields_not_same_for_case_dual(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is the same as the original (although a different object in memory) and differs in that specific attribute.

 
    """

    cstack1 = self.fixture[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
    Z = sg.Ax('Z')
   
    coord3_copy = coord3.copy(dual = coord2)

    test_args = {'name':'joep', 'value':np.array([1.,2.,3.]),'dual':coord2,'axis':Z,'direction':'Z','units':'cm','long_name':'this is a coordinate in the x direction','metadata':{'hi':0},'strings':['five','one','two','three','four']}
 

    for ta in test_args:
      value = test_args[ta]
      coord3_copy = coord3.copy(**{ta:value})

      coord_att = getattr(coord3_copy,ta)
      if isinstance(coord_att,np.ndarray):
        self.assertEqual(np.array_equal(coord_att, value), True  )
      else:
        self.assertEqual(coord_att, value  )



  def test_coord_copy_self_dual(self):
    """
    Test copy method of Coord. 
    """

    cstack1 = self.fixture[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
 
    # coord1 and coord2 have non-self duals, coord3 is self-dual. 

    copy_coord3 = coord3.copy()

    # test whether coord3 remains self-dual under operation:
    self.assertEqual(copy_coord3.dual is copy_coord3.dual , True )



  def test_coord_copy_other_dual(self):
    """
    Test copy method of Coord. 
    """

    cstack1 = self.fixture[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
 
    # coord1 and coord2 have non-self duals, coord3 is self-dual. 

    copy_coord2 = coord2.copy()

 
    self.assertEqual(copy_coord2.dual is coord2.dual , True )




  def test_coord_neg_self_dual(self):
    """
    Test __neg__ method of Coord. 
    """

    cstack1 = self.fixture[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
 
    # coord1 and coord2 have non-self duals, coord3 is self-dual. 

    minus_coord3 = -coord3

    # test whether coord3 remains self-dual under operation:
    self.assertEqual(minus_coord3.dual, minus_coord3 )



  def test_coord_neg_other_dual(self):
    """
    Test __neg__ method of Coord. 
    """

    cstack1 = self.fixture[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
 
    # coord1 and coord2 have non-self duals:

    minus_coord2 = -coord2

    self.assertEqual( np.array_equal(minus_coord2.dual.value, -coord2.dual.value), True )





  def test_same_method_yields_same(self):
    """
    Test whether making a copy with no arguments passed to .copy method yields a Coord object that is the same (with respect to .same method) as the original (although a different object in memory). Also tested for other Coord objects from fixture and for hybrid axis attributes (one str, one Ax).
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord4 = self.fixture[1][0]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy()

    self.assertEqual(coord3.same(coord3_copy),True  )
    self.assertEqual(coord1.same(coord4),True  )

    coord4.axis = sg.fieldcls.Ax(coord4.axis)
    self.assertEqual(coord1.same(coord4),True  )
    self.assertEqual(coord4.same(coord1),True  )
    self.assertEqual(coord1.same(coord3),False  )

  def test_same_method_yields_not_same_for_case_array(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[0]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(value = np.array([5,6,7]))

    self.assertEqual(coord3.same(coord3_copy), False  )

  def test_same_method_yields_not_same_for_case_name(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[0]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(name = 'joep')

    self.assertEqual(coord3.same(coord3_copy), False  )

  def test_same_method_yields_not_same_for_case_axis(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[0]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(axis = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )

  def test_same_method_yields_not_same_for_case_direction(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[0]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(direction = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )


  def test_same_method_yields_not_same_for_case_direction(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[0]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(direction = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )



  def test_cast_method_2D_grid(self):
    """
    Test Coord cast method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]     
    coord3 = cstack1[2]

    F = coord1.cast(coord1*coord2)
  
    self.assertEqual(F.shape, (3,4)  )
    self.assertEqual( np.array_equal( F.value[1,:], np.array([2,2,2,2]) ), True  )
    self.assertEqual( np.array_equal( F.value[:,1], np.array([1.,2.,3.]) ), True  )


  def test_copy_equiv_method(self):
    """
    Test whether Coord.copy yields a new self-equivalent Coord object.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    K = coord1.copy(name='ho')  

    self.assertEqual(K.is_equiv(K),  True  )


  def test_make_equiv_method(self):
    """
    Test Coord make_equiv method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]     
    coord3 = cstack1[2]

    coord1.make_equiv(coord2)
  
    self.assertEqual(coord2.associative, coord1.associative  )

  def test_is_equiv_method_false(self):
    """
    Test Coord is_equiv method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]     
    coord3 = cstack1[2]
  
    self.assertEqual(coord1.is_equiv(coord2),  False  )
 

  def test_is_equiv_method_true(self):
    """
    Test Coord is_equiv method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]     
    coord3 = cstack1[2]

    coord1.make_equiv(coord2)
  
    self.assertEqual(coord1.is_equiv(coord2),  True  )
    self.assertEqual(coord2.is_equiv(coord1),  True  )


  def test_eq_in_method_false(self):
    """
    Test Coord eq_in method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]     
    coord3 = cstack1[2]
 
    self.assertEqual(coord1.eq_in(coord2*coord3),  None  )

  def test_eq_in_method_true(self):
    """
    Test Coord eq_in method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]     
    coord3 = cstack1[2]

    coord1.make_equiv(coord2)
  
    self.assertEqual(coord1.eq_in(coord2*coord3),  coord2  )


  def test_pow_method(self):
    """
    Test Coord __pow__ method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
  
    self.assertEqual(coord1**2,  sg.Gr((coord1,))  )

  def test_mul_method_non_equiv(self):
    """
    Test Coord __mul__ method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
  
    self.assertEqual((coord1*coord2).shape(),  (3,4)  )

  def test_mul_method_equiv(self):
    """
    Test Coord __mul__ method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
      
    coord1.make_equiv(coord2)

    self.assertEqual((coord1*coord2).shape(),  (3,)  )

  def test_roll_function_non_masked(self):
    """Test the sg roll function
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    K = coord1(coord1*coord2)
    R= sg.roll(K,coord=coord1,mask = False)

    # test whether Field.roll method compatible with sg.roll.
    self.assertEqual( np.array_equal(R.value, K.roll(shift=1 , crd=coord1).value),True)
    
    self.assertEqual( np.array_equal( R.value[0,:],  np.array([3.,3.,3.,3.]) ), True  )
    self.assertEqual( np.array_equal( R.value[1,:],  np.array([1.,1.,1.,1.]) ), True  )
    # first coord in R.grid is replaced:
    self.assertEqual( R.grid[0] is coord1, False  )
    # second coord in R.grid is not replaced:
    self.assertEqual( R.grid[1] is coord2, True  )

    # test whether coord in R.grid is properly rolled:
    self.assertEqual( np.array_equal(R.grid[0].value , np.array( [3., 1., 2.] )  ) , True  )

  def test_roll_function_non_masked_keepgrid(self):
    """Test the sg roll function
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    K = coord1(coord1*coord2)
    R= sg.roll(K,coord=coord1,mask = False, keepgrid = True)
    
    # first coord in R.grid is not replaced:
    self.assertEqual( R.grid[0] is coord1, True  )
    # second coord in R.grid is not replaced:
    self.assertEqual( R.grid[1] is coord2, True  )


  def test_roll_function_masked(self):
    """Test the sg roll function
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    K = coord1(coord1*coord2)
    R= sg.roll(K,coord=coord1,mask = True)
    
    self.assertEqual( np.isnan( R.value[0,:] ).all() , True  )
    self.assertEqual( np.array_equal( R.value[1,:],  np.array([1.,1.,1.,1.]) ), True  )


  def test_coord_shift_method(self):
    """
    Test Coord coord_shift method.

    Need to check the nan's that show up as numbers in the exposed area.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
     
    K = coord1(coord1*coord2) 
    R = coord1.coord_shift(K,1)

    self.assertEqual( np.isnan( R.value[0,:] ).all() , True  )
    self.assertEqual( np.array_equal( R.value[1,:],  np.array([1.,1.,1.,1.]) ), True  )

  def test_trans_method(self):
    """
    Test Coord trans method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
     
    K = coord1(coord1*coord2) 
    R = coord1.trans(K)

    self.assertEqual( np.array_equal( R.value[1,:],  np.array([1.,1.,1.,1.]) ), True  )
    self.assertEqual( np.array_equal( R.value[2,:],  np.array([1.,1.,1.,1.]) ), True  )

    self.assertEqual( R.grid[0] is coord1   , True  )


  def test_sum_method(self):
    """
    Test Coord sum method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
     
    K = coord1(coord1*coord2) 
    R = coord1.sum(K)

    self.assertEqual( np.array_equal(R.value,  6*np.array([1.,1.,1.,1.]) ), True  )

    self.assertEqual( coord1.sum(sg.ones(coord1**2)) , 3.0 )

    # should error if coord1 not in gr of argument field:
    self.assertRaises( ValueError, coord1.sum, sg.ones(coord2*coord3) )     

  def test_roll_method(self):
    """
    Test Coord roll method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    K = coord1.roll(1)

    self.assertEqual( np.array_equal(coord1,  np.array([1.,2.,3.]) ), True  )
    self.assertEqual( np.array_equal(K.value,  np.array([3.,1.,2.]) ), True  )

  def test_flip_method(self):
    """
    Test Coord flip method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    K = coord1(coord1*coord2) 
    R = coord1.flip(K)

    self.assertEqual( np.array_equal(R.value[:,1],  np.array([3.,2.,1.]) ), True  )

  def test_flip_method_transpose_of_previous(self):
    """
    Test Coord flip method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    # order of coord product reversed with respect to previous test:
    K = coord1(coord2*coord1) 
    R = coord1.flip(K)

    self.assertEqual( np.array_equal(R.value[1,:],  np.array([3.,2.,1.]) ), True  )


  def test_cumsum_method(self):
    """
    Test Coord cumsum method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2]

    # order of coord product reversed with respect to previous test:

    ONES = sg.ones(coord1*coord2)

    R = coord1.cumsum(ONES  )

    self.assertEqual(R.grid,ONES.grid)
    self.assertEqual( np.array_equal(R.value[0,:],  np.array([3.,3.,3.,3.]) ), True  )

    R = coord1.cumsum(ONES , upward = True )

    self.assertEqual( np.array_equal(R.value[-1,:],  np.array([3.,3.,3.,3.]) ), True  )

    # should error if coord1 not in gr of argument field:
    self.assertRaises( ValueError, coord1.cumsum, sg.ones(coord2*coord3) )     


  def test_der_method(self):
    """
    Test Coord der method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
  
    # make up Ax to use for coord1:
    W = sg.fieldcls.Ax('W', direction='W')
    coord1.give_axis(W)

    W2 = sg.fieldcls.Ax('W2', direction='W2')
    coord1.give_axis(W2)

    K = coord1(coord1*coord2)
    R=coord1.der(K)
    R_with_Ax_method = W.der(K)
    R_with_Field_method = K.der(W)

    value_R = copy.deepcopy(R.value)
    value_R_wam = copy.deepcopy(R_with_Ax_method.value)
    value_R_wfm = copy.deepcopy(R_with_Field_method.value)

    value_R[np.isnan(value_R)] = 0. # to be able to do np.array_equal
    value_R_wam[np.isnan(value_R_wam)] = 0.
    value_R_wfm[np.isnan(value_R_wfm)] = 0.

    # test whether Ax.der calls coord1.der properly:
    self.assertEqual(np.array_equal(value_R, value_R_wam), True)

    # test whether Field.der calls Ax.der properly:
    self.assertEqual(np.array_equal(value_R, value_R_wfm), True)

    self.assertEqual( np.array_equal(  R.value[1,:], np.array([1.,1., 1.,   1.]) ), True )

    self.assertEqual( (np.isnan(  R.value )[0,:]).all(), True )

  def test_delta_dist_method(self):
    """
    Test Coord delta_dist method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
 
    self.assertEqual( np.array_equal(  coord1.delta_dist().value[1:], np.array([ 1.,   1.]) ), True )

    self.assertEqual( np.isnan(  coord1.delta_dist()[0] ), True )


  def test_d_method(self):
    """
    Test Coord d method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
 
    self.assertEqual( np.array_equal(  coord1.d().value, np.array([ 1., 1.,   1.]) ), True )

  def test_vol_method(self):
    """
    Test Coord vol method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2] 

    self.assertEqual( np.array_equal(  coord1.vol(coord1*coord2).value, np.array([ 1., 1.,   1.]) ), True )

    self.assertRaises(ValueError, coord1.vol , coord2*coord3 )

# -------- test block for YCoord class ---------------

  def testY_copy_method_yields_not_same_for_case_name(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is the same as the original (although a different object in memory) and differs in that specific attribute.

    """

    cstack1 = self.fixture[2]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(name = 'joep')

    self.assertEqual(coord3_copy.name, 'joep'  )

  def testY_copy_method_yields_not_same_for_case_dual(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is the same as the original (although a different object in memory) and differs in that specific attribute.

 
    """

    cstack1 = self.fixture[2]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
    Z = sg.Ax('Z')
   
    coord3_copy = coord3.copy(dual = coord2)

    test_args = {'name':'joep', 'value':np.array([1.,2.,3.]),'dual':coord2,'axis':Z,'direction':'Z','units':'cm','long_name':'this is a coordinate in the x direction','metadata':{'hi':0},'strings':['five','one','two','three','four']}
 

    for ta in test_args:
      value = test_args[ta]
      coord3_copy = coord3.copy(**{ta:value})

      coord_att = getattr(coord3_copy,ta)
      if isinstance(coord_att,np.ndarray):
        self.assertEqual(np.array_equal(coord_att, value), True  )
      else:
        self.assertEqual(coord_att, value  )





  def test_Ysame_method_yields_same(self):
    """
    Test whether making a copy with no arguments passed to .copy method yields a Coord object that is the same (with respect to .same method) as the original (although a different object in memory).
    """

    cstack1 = self.fixture[2]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy()

    self.assertEqual(coord3.same(coord3_copy),True  )

  def testY_same_method_yields_not_same_for_case_array(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[2]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(value = np.array([5,6,7]))

    self.assertEqual(coord3.same(coord3_copy), False  )

  def testY_same_method_yields_not_same_for_case_name(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[2]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(name = 'joep')

    self.assertEqual(coord3.same(coord3_copy), False  )

  def testY_same_method_yields_not_same_for_case_axis(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[2]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(axis = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )

  def testY_same_method_yields_not_same_for_case_direction(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[2]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(direction = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )


  def testY_same_method_yields_not_same_for_case_direction(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[2]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(direction = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )








# -------- test block for XCoord class ---------------

  def testX_copy_method_yields_not_same_for_case_name(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is the same as the original (although a different object in memory) and differs in that specific attribute.

    """

    cstack1 = self.fixture[4]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(name = 'joep')

    self.assertEqual(coord3_copy.name, 'joep'  )

  def testX_copy_method_yields_not_same_for_case_dual(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is the same as the original (although a different object in memory) and differs in that specific attribute.

 
    """

    cstack1 = self.fixture[4]
    coord2 = cstack1[1]
    coord3 = cstack1[2]
    Z = sg.Ax('Z')
   
    coord3_copy = coord3.copy(dual = coord2)

    test_args = {'name':'joep', 'value':np.array([1.,2.,3.]),'dual':coord2,'axis':Z,'direction':'Z','units':'cm','long_name':'this is a coordinate in the x direction','metadata':{'hi':0},'strings':['five','one','two','three','four']}
 

    for ta in test_args:
      value = test_args[ta]
      coord3_copy = coord3.copy(**{ta:value})

      coord_att = getattr(coord3_copy,ta)
      if isinstance(coord_att,np.ndarray):
        self.assertEqual(np.array_equal(coord_att, value), True  )
      else:
        self.assertEqual(coord_att, value  )



  def testX_roll_method(self):
    """
    Test XCoord roll method.
    """

    cstack1 = self.fixture[4]

    coord1 = cstack1[0]
 

    # Check the shift is re-entrant:
    self.assertEqual( np.array_equal(coord1.roll(1).value,  np.array([-357., 1., 2.])) , True  )
    self.assertEqual( np.array_equal(coord1.roll(-1).value,  np.array([-358., -357., 1.])) , True  )

  def testX_coord_shift_method(self):
    """
    Test XCoord coord_shift method.

    Need to check no nan's show up in the exposed area.
    """

    cstack1 = self.fixture[4]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
     
    K = coord1(coord1*coord2) 
    R = coord1.coord_shift(K,1)

    # Check the shift is re-entrant:
    self.assertEqual( np.array_equal( R.value[0,:],  np.array([3.,3.,3.,3.]) ), True  )
    self.assertEqual( np.array_equal( R.value[1,:],  np.array([1.,1.,1.,1.]) ), True  )

    # check sum is preserved:
    self.assertEqual( np.sum(K.value), np.sum(R.value)  )



  def testX_delta_dist_method(self):
    """ Test the XCoord delta_dist method 
    """

    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    K = xcoord1.delta_dist(ycoord1)

    # This must be a 2D Field:
    self.assertEqual(K.shape, (7,4))

    # test some values
    self.assertAlmostEqual(K[1,0], 5002986.3008417469, places =3  )
    self.assertAlmostEqual(K[3,0], 10005972.601683492, places =3  )

    # kind of "checksum"
    self.assertAlmostEqual(np.sum(K.value), 149371192.51449975, places =3  )


    # distances between all consecutive points around circle must be same
    # test whether constant in this direction, and no nan:
    K2 = K - K[1,0]
    idx = K2.value == 0

    self.assertEqual( idx[1,:].all(), True )


# I want more tests in this area. Also more tests to indicate that actual calculations (not just code logic) are correct.

  def testX_der_method(self):
    """ Test the XCoord der method 
    """

    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    xcoord2 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    K = xcoord1.delta_dist(ycoord1)

    R=xcoord1.der(K,ycoord1)

    # This must be a 2D Field:
    self.assertEqual(R.shape, (7,4))

    Idx = R.value == 0.

    # result must be all zero
    self.assertEqual(Idx.all(),True)
  
    # note that xcoord2 is identical to xcoord1, but not the same object in memory, hence error:

    self.assertRaises(ValueError, xcoord2.der, **{'F':K,'y_coord':ycoord1})



  def testX_dist_method(self):
    """ Test the XCoord dist method 
    """

    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    K = xcoord1.dist(ycoord1)

    # This must be a 2D Field:
    self.assertEqual(K.shape, (7,4))

    # test some value:
    self.assertAlmostEqual(K[2,2], 17330852.925257955 , places =3  )

    # kind of "checksum"
    self.assertAlmostEqual(np.sum(K.value), 224056788.77174962 , places =3  )


    R = xcoord1.der(K,ycoord1)
    dR = R.value[:,1:] - 1.

    # result must be all be ~zero
    self.assertEqual(np.max(dR) < 1e-15,True)

    dR = R.value[:,:1] -3.

    # result must be all be ~zero
    self.assertEqual(np.max(dR) < 1e-15,True)

  
    # note that xcoord2 is identical to xcoord1, but not the same object in memory, hence error:
    xcoord2 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    self.assertRaises(ValueError, xcoord2.der, **{'F':K,'y_coord':ycoord1})



  def testX_d_method(self):
    """ Test the XCoord d method 
    """

    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    xcoord1_edges = sg.fieldcls.XCoord(name = 'testx_edges',direction ='X',value = np.arange(0.,360+45.,90.) -45. , dual = xcoord1  )
#    ycoord1_edges = sg.fieldcls.YCoord(name = 'testy_edges',direction ='Y',value = np.arange(-90.+y_step/2,90.,y_step) , dual = ycoord1  )

    K = xcoord1.d(ycoord1)

    # This must be a 2D Field:
    self.assertEqual(K.shape, (7,4))

    self.assertAlmostEqual(np.sum(K.value), 149371192.51449975 , places =3  )

  def testX_vol_method(self):
    """ Test the XCoord vol method.

        Might want to extend this with the introduction of new derived Coord classes.
    """

    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    xcoord1_edges = sg.fieldcls.XCoord(name = 'testx_edges',direction ='X',value = np.arange(0.,360+45.,90.) -45. , dual = xcoord1  )
#    ycoord1_edges = sg.fieldcls.YCoord(name = 'testy_edges',direction ='Y',value = np.arange(-90.+y_step/2,90.,y_step) , dual = ycoord1  )

    # a YCoord must be in the grid:
    self.assertRaises(RuntimeError, xcoord1.vol, xcoord1**2 )

    # Now a YCoord is in the grid:
    K = xcoord1.vol(ycoord1*xcoord1)

    # This must be a 2D Field:
    self.assertEqual(K.shape, (7,4))

    self.assertEqual( np.array_equal( K.value, xcoord1.d(ycoord1).value ) , True)

    xcoord2 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )

    # identically valued coord2 is a different object to those in the grid, hence no go:
    self.assertEqual(xcoord2.vol(ycoord1*xcoord1) , None )



  def test_Xsame_method_yields_same(self):
    """
    Test whether making a copy with no arguments passed to .copy method yields a Coord object that is the same (with respect to .same method) as the original (although a different object in memory).
    """

    cstack1 = self.fixture[4]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy()

    self.assertEqual(coord3.same(coord3_copy),True  )

  def testX_same_method_yields_not_same_for_case_array(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[4]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(value = np.array([5,6,7]))

    self.assertEqual(coord3.same(coord3_copy), False  )

  def testX_same_method_yields_not_same_for_case_name(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[4]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(name = 'joep')

    self.assertEqual(coord3.same(coord3_copy), False  )

  def testX_same_method_yields_not_same_for_case_axis(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[4]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(axis = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )

  def testX_same_method_yields_not_same_for_case_direction(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[4]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(direction = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )


  def testX_same_method_yields_not_same_for_case_direction(self):
    """
    Test whether making a copy with 1 argument passed to .copy method yields a Coord object that is NOT the same (with respect to .same method) as the original (and a different object in memory).

    Note that in general, the .same method tests for:

    self.array_equal(other)
    self.name == other.name
    self.axis == other.axis
    self.direction == other.direction 

    """

    cstack1 = self.fixture[4]
    coord3 = cstack1[2]
   
    coord3_copy = coord3.copy(direction = 'Z')

    self.assertEqual(coord3.same(coord3_copy), False  )












# -----------------------

# ------- further general Coord tests --------

  def test_sort(self):
    cstack1 = self.fixture[0]
    coord3 = cstack1[2]
    coord3.sort()
    value = copy.deepcopy(coord3.value)
    value.sort()
    self.assertEqual( np.array_equal(coord3.value , value ) ,True )

  def test_equality_relation_weaksame(self):
    """"
    Does the &-relationship yield equality?
    """
    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    # First two coord objects should have same content

    self.assertEqual(cstack1[0].weaksame(cstack2[0]), True)

  def test_inequality_relation_weaksame(self):
    """"
    Does the &-relationship yield inequality?
    """

    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    # These two coord objects are not the same

    self.assertEqual(cstack1[0].weaksame(cstack2[1]), False)

  
  def test_equality_relation_weaksame_grid(self):
    """"
    Does the weaksame-relationship yield equality for multiple-member object?
    """
    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    # First two coord objects should have same content

    self.assertEqual( (cstack1[0]*cstack1[1]).weaksame(cstack2[0]*cstack2[1]), True)



  def test_inequality_relation_weaksame_grid(self):
    """"
    Does the weaksame-relationship yield inequality for multiple-member object?
    """
    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    # First two coord objects should have same content

    self.assertEqual( (cstack1[0]*cstack1[1]).weaksame(cstack2[1]*cstack2[0]), False)

    self.assertEqual( (cstack1[0]*cstack1[1]).weaksame(cstack2[1]*cstack2[2]), False)





# ----- some make_axes related tests:


  def test_equality_relation_find_equal_axes(self):
    """"
    Does the function find_equal_axes recognise equivalent coord objects in the two cstacks and replace the elements of the 2nd stack accordingly?
    """

    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    # this should remove all redundant coord objects with respect to &-equality
    sg.find_equal_axes(cstack1,cstack2)

    self.assertEqual(cstack1,cstack2)

  def test_make_axes_function_type_output(self):
    """
    The output should be a list of Ax objects ([X,Y] expected, see below)
    """

    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    self.assertEqual(isinstance(sg.make_axes(cstack1 + cstack2)[0],sg.Ax ) , True   )


  def test_make_axes_function_output_expected(self):
    """
    The test coords contain only X and Y direction Ax objects
    """

    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    self.assertEqual(str( sg.make_axes(cstack1 + cstack2) ) , '[X, Y]'  )


  def test_make_axes_function_no_output_expected(self):
    """
    Calling make_axes twice should not yield further output
    """

    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]
    sg.make_axes(cstack1 + cstack2) 
    self.assertEqual(str( sg.make_axes(cstack1 + cstack2) ) , '[]'  )


class TestAxAndAxGr(unittest.TestCase):

# ----- for Ax and AxGr objects (not using fixture)

  def test_equality_relation_weaksame_grid(self):
    """"
    Does the weaksame-relationship yield equality for multiple-member object?
    """
    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')

    X2 = sg.fieldcls.Ax('X')
    Y2 = sg.fieldcls.Ax('Y')


    # First two coord objects should have same content

    self.assertEqual( (X*Y).weaksame(Y*X), False)
    self.assertEqual( (Y*X).weaksame(Y2*X2), True)

# ----- for Ax and AxGr objects (not using fixture)

  def test_copy_AxGr(self):
    """"
    Test copy method of AxGr
    """
    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')

    ag_copy = (X*Y).copy()


    self.assertEqual( ag_copy.__repr__(), '(X,Y,)')


  def test_eq_in_AxGr(self):
    """"
    Test eq_in method of AxGr
    """
    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')
    Z = sg.fieldcls.Ax('Z')
  
    self.assertEqual( (X*Y).eq_in(X), True )
    self.assertEqual( (X*Y).eq_in(Z), False )

    X2 = sg.fieldcls.Ax('X2')

    self.assertEqual( (X*Y).eq_in(X2), False )
    X2.make_equiv(X)
    self.assertEqual( (X*Y).eq_in(X2), True )


  def test_eq_index_AxGr(self):
    """"
    Test eq_index method of AxGr
    """
    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')
    Z = sg.fieldcls.Ax('Z')
  
    self.assertEqual( (X*Y).eq_index(Y), 1 )
    self.assertEqual( (X*Y).eq_index(Z), -1 )

    X2 = sg.fieldcls.Ax('X2')

    self.assertEqual( (X*Y).eq_index(X2), -1 )
    X2.make_equiv(X)
    self.assertEqual( (X*Y).eq_index(X2), 0 )


  def test_eq_perm_AxGr(self):
    """"
    Test eq_perm method of AxGr
    """
    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')
    Z = sg.fieldcls.Ax('Z')
  
    self.assertEqual( (X*Y).eq_perm(Y*X), (1,0) )
    self.assertEqual( (X*Y).eq_perm(Y*Z), None )

    X2 = sg.fieldcls.Ax('X2')

    self.assertEqual( (X*Y).eq_perm(Y*X2), None )
    X2.make_equiv(X)
    self.assertEqual( (X*Y).eq_perm(Y*X2), (1,0) )




  def test_ax_div_mult(self):

    # set up some independent axes to test on:
    a1 = sg.fieldcls.Ax(name='a1')
    a2 = sg.fieldcls.Ax(name='a2')
    a3 = sg.fieldcls.Ax(name='a3')
    a4 = sg.fieldcls.Ax(name='a4')

    self.assertEqual(len(a1*a2*a3),3)
    self.assertEqual((a1*a2*a3)/a2 , a1*a3  )
    self.assertEqual((a1*a3)/a2 , a1*a3  )
  

class TestGr(unittest.TestCase):



  def setUp(self):
    print 'Setting up %s'%type(self).__name__

    def provide_axis(cstack):
      for i, c in enumerate(cstack):
        cstack[i].axis = cstack[i].direction 

      return cstack
  
    # Note that some coord values are deliberately unordered.   

    # Coords ---    
    coord1 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) , metadata = {'hi':5} )
    coord2 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]), metadata = {'hi':7})
    coord3 = sg.fieldcls.Coord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':3})
    # identical in main attributes to previous set (in order):
    coord4 = sg.fieldcls.Coord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]), metadata = {'hi':8})
    coord5 = sg.fieldcls.Coord(name = 'test2',direction ='Y',value =np.array([1,2,3, 4]), metadata = {'hi':10})
    coord6 = sg.fieldcls.Coord(name = 'test',direction ='X',value =np.array([5,1,2,3, 4]), metadata = {'hi':12})

    # providing coord1 and coord2 with duals. coord3 is self-dual

    coord1_edges = sg.fieldcls.Coord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = coord1 , metadata = {'hi':25} )
    coord2_edges = sg.fieldcls.Coord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = coord2, metadata = {'hi':77})

    # identical in main attributes to previous set (in order):
    coord4_edges = sg.fieldcls.Coord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = coord4 , metadata = {'hi':25} )
    coord5_edges = sg.fieldcls.Coord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = coord5, metadata = {'hi':77})


    # YCoords ---

    ycoord1 = sg.fieldcls.YCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) , metadata = {'hi':5} )
    ycoord2 = sg.fieldcls.YCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]), metadata = {'hi':7})
    ycoord3 = sg.fieldcls.YCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':3})
    # identical in main attributes to previous set (in order):
    ycoord4 = sg.fieldcls.YCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]), metadata = {'hi':8})
    ycoord5 = sg.fieldcls.YCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3., 4.]), metadata = {'hi':10})
    ycoord6 = sg.fieldcls.YCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3., 4.]), metadata = {'hi':12})


    ycoord1_edges = sg.fieldcls.YCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = ycoord1 , metadata = {'hi':25} )
    ycoord2_edges = sg.fieldcls.YCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = ycoord2, metadata = {'hi':77})

    # identical in main attributes to previous set (in order):
    ycoord4_edges = sg.fieldcls.YCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = ycoord4 , metadata = {'hi':25} )
    ycoord5_edges = sg.fieldcls.YCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = ycoord5, metadata = {'hi':77})




    # XCoords ---

    xcoord1 = sg.fieldcls.XCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]) , metadata = {'hi':5} )
    xcoord2 = sg.fieldcls.XCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3.,4.]), metadata = {'hi':7})
    xcoord3 = sg.fieldcls.XCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':3})
    # identical in main attributes to previous set (in order):
    xcoord4 = sg.fieldcls.XCoord(name = 'test1',direction ='X',value =np.array([1.,2.,3.]), metadata = {'hi':8})
    xcoord5 = sg.fieldcls.XCoord(name = 'test2',direction ='Y',value =np.array([1.,2.,3., 4.]), metadata = {'hi':10})
    xcoord6 = sg.fieldcls.XCoord(name = 'test',direction ='X',value =np.array([5.,1.,2.,3., 4.]), metadata = {'hi':12})

    xcoord1_edges = sg.fieldcls.XCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = xcoord1 , metadata = {'hi':25} )
    xcoord2_edges = sg.fieldcls.XCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = xcoord2, metadata = {'hi':77})

    # identical in main attributes to previous set (in order):
    xcoord4_edges = sg.fieldcls.XCoord(name = 'test1_edges',direction ='X',value =np.array([0.5,1.5,2.5,3.5]), dual = xcoord4 , metadata = {'hi':25} )
    xcoord5_edges = sg.fieldcls.XCoord(name = 'test2_edges',direction ='Y',value =np.array([0.5,1.5,2.5,3.5,4.5]), dual = xcoord5, metadata = {'hi':77})


    # we are testing for Coord, YCoord and XCoord 
    cstack1 = provide_axis([coord1,coord2,coord3,coord1_edges,coord2_edges])
    cstack2 = provide_axis([coord4,coord5,coord6,coord4_edges,coord5_edges])

    ycstack1 = provide_axis([ycoord1,ycoord2,ycoord3,ycoord1_edges,ycoord2_edges])
    ycstack2 = provide_axis([ycoord4,ycoord5,ycoord6,ycoord4_edges,ycoord5_edges])

    xcstack1 = provide_axis([xcoord1,xcoord2,xcoord3,xcoord1_edges,xcoord2_edges])
    xcstack2 = provide_axis([xcoord4,xcoord5,xcoord6,xcoord4_edges,xcoord5_edges])

    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')

    coord1.axis = X
    coord2.axis = Y
    
    coord4.axis = X
    coord5.axis = Y



    self.fixture = [cstack1, cstack2, ycstack1, ycstack2,xcstack1, xcstack2,]

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture


  def test_copy_Gr(self):
    """"
    Test copy method of Gr
    """
 
    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    gr_copy = (coord1*coord2).copy()

    self.assertEqual( gr_copy.__repr__(), '(test1, test2)')

  def test_Gr_array_equal_method(self):
    """"
    Test array_equal method of Gr
    """

    cstack1 = self.fixture[0]
    cstack2 = self.fixture[1]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    coord4 = cstack2[0]
    coord5 = cstack2[1]

    self.assertEqual( (coord1*coord2).array_equal(coord4*coord5), [True ,True] )


  def test_Gr_axis_method(self):
    """"
    Test axis method of Gr
    """
 
    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]

    self.assertEqual( (coord1*coord2).axis().__repr__(), '(X,Y,)')


  def test_Gr_axis_method(self):
    """
    Test axis method of Coord. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]
 
    X = sg.fieldcls.Ax('X')
    Y = sg.fieldcls.Ax('Y')


    coord1.give_axis(X)
    coord2.give_axis(Y)

    self.assertEqual( (coord1*coord2).axis() , X*Y  )
 


  def test_Gr_reverse_method(self):
    """
    Test reverse method of Coord. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]
 
    self.assertEqual( (coord1*coord2).reverse() , coord2*coord1  )




  def test_Gr_is_equiv_method(self):
    """
    Test is_equiv method. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]

    cstack2 = self.fixture[1]
    coord4 = cstack2[0]
    coord5 = cstack2[1]

 
    self.assertEqual( (coord1*coord2).is_equiv(coord5*coord4) , False  )

    coord1.make_equiv(coord4)
    coord2.make_equiv(coord5)
 
    self.assertEqual( (coord1*coord2).is_equiv(coord5*coord4) , True  )



  def test_Gr_eq_in_method(self):
    """
    Test Gr.eq_in method. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]

    cstack2 = self.fixture[1]
    coord4 = cstack2[0]
    coord5 = cstack2[1]

 
    self.assertEqual( (coord1*coord2).eq_in(coord4) , False  )

    coord1.make_equiv(coord4)
#    coord2.make_equiv(coord5)
 
    self.assertEqual( (coord1*coord2).eq_in(coord4) , True  )


  def test_Gr_rearrange_method(self):
    """
    Test Gr.rearrange method. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]

    self.assertEqual( (coord1*coord2).rearrange([1,0]) , coord2*coord1  )



  def test_Gr_perm_method(self):
    """
    Test perm method of Gr. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]

    cstack2 = self.fixture[1]
    coord4 = cstack2[0]
    coord5 = cstack2[1]


    self.assertEqual( (coord1*coord2).perm(coord2*coord1) , (1,0)  )
 
    self.assertEqual( (coord1*coord2).perm(coord5*coord4) is None , True  )

    coord1.make_equiv(coord4)
    coord2.make_equiv(coord5)
 
    self.assertEqual( (coord1*coord2).perm(coord5*coord4) , None  )


  def test_Gr_eq_perm_method(self):
    """
    Test eq_perm method of Gr. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]

    cstack2 = self.fixture[1]
    coord4 = cstack2[0]
    coord5 = cstack2[1]


    self.assertEqual( (coord1*coord2).eq_perm(coord2*coord1) , (1,0)  )
 
    self.assertEqual( (coord1*coord2).eq_perm(coord5*coord4) is None , True  )

    coord1.make_equiv(coord4)
    coord2.make_equiv(coord5)
 
    self.assertEqual( (coord1*coord2).eq_perm(coord5*coord4) , (1,0)  )


  def test_Gr_shape_method(self):
    """
    Test Gr.shape method. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]


    self.assertEqual( (coord1*coord2).shape() , (3,4)  )
 

  def test_Gr_ones_method(self):
    """
    Test Gr.ones method. 
    """

    cstack1 = self.fixture[0]
    coord1 = cstack1[0]
    coord2 = cstack1[1]

    cstack2 = self.fixture[1]
    coord4 = cstack2[0]
    coord5 = cstack2[1]

    K = (coord1*coord2).ones()

    self.assertEqual( K.value.shape , (3,4)  )
 


  def test_Gr_der_method(self):
    """
    Test Gr.der method. 
    """


    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    xcoord2 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    # construct simple test Field:
    K = xcoord1.delta_dist(ycoord1)

    # take derivative:
    R= (ycoord1*xcoord1).der(xcoord1, K)

    # This must be a 2D Field:
    self.assertEqual(R.shape, (7,4))

    Idx = R.value == 0.

    # Few tests to distinguish between possible problem causes:
    self.assertEqual(Idx.all(),True)
    self.assertEqual( np.array_equal( R.value, xcoord1.der(K,ycoord1).value ) , True )



  def test_vol_method(self):
    """
    Test Coord vol method.
    """

    cstack1 = self.fixture[0]

    coord1 = cstack1[0]
    coord2 = cstack1[1]
    coord3 = cstack1[2] 

    self.assertEqual( np.array_equal(  coord1.vol(coord1*coord2).value, np.array([ 1., 1.,   1.]) ), True )

    self.assertRaises(ValueError, coord1.vol , coord2*coord3 )




  def test__find_args_coord_method(self):
    """
    Test Coord _find_args_coord method.
    """


    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    xcoord2 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    # construct simple test Field:
    K = xcoord1.delta_dist(ycoord1)

    # take derivative:
    R= (ycoord1*xcoord1).der(xcoord1, K)

    # This must be a 2D Field:
    self.assertEqual(R.shape, (7,4))

    Idx = R.value == 0.

    # Few tests to distinguish between possible problem causes:
    self.assertEqual(Idx.all(),True)
    self.assertEqual( np.array_equal( R.value, xcoord1.der(K,ycoord1).value ) , True )


    A = (ycoord1*xcoord1)._find_args_coord({'x_coord':sg.fieldcls.XCoord,'y_coord':sg.fieldcls.YCoord,'z_coord':sg.fieldcls.Coord})

    self.assertEqual(A,[[], [ycoord1]] )


    A = (xcoord1*ycoord1)._find_args_coord({'x_coord':sg.fieldcls.XCoord,'y_coord':sg.fieldcls.YCoord,'z_coord':sg.fieldcls.Coord})

    self.assertEqual(A,[ [ycoord1] , []  ] )





  def test__call_on_members_Gr_method(self):
    """
    Test Coord _find_args_coord method.
    """


    y_step=30;

    xcoord1 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    xcoord2 = sg.fieldcls.XCoord(name = 'testx',direction ='X',value =np.arange(0.,360.,90.) )
    ycoord1 = sg.fieldcls.YCoord(name = 'testy',direction ='Y',value =np.arange(-90.,90.+y_step,y_step) )

    grid = ycoord1*xcoord1

    R= grid.call_on_members('__neg__')

    # This must be a 2D Field:
    self.assertEqual(R.shape(), (7,4))


    self.assertEqual( np.array_equal( R[0].value
 , -grid[0].value )  ,True)

    self.assertEqual( np.array_equal( R[1].value
 , -grid[1].value )  ,True)

# ------------- Test utilsg.py module --------------------


class TestUtilsg(unittest.TestCase):

  def test_id_index_id_in_rem_equivs_functions(self):
    """Test id_index, id_in and rem_equivs from sg.utilsg.
    """

    # set up some axes to test on:
    a1 = sg.fieldcls.Ax(name='a1')
    a2 = sg.fieldcls.Ax(name='a2')
    a3 = sg.fieldcls.Ax(name='a3')
    a4 = sg.fieldcls.Ax(name='a4')

    # b2 is equivalent to a2, b3 to none.
    b2 = sg.fieldcls.Ax(name='a2', direction ='Q', long_name='Q')
    b3 = sg.fieldcls.Ax(name='b3', direction ='Q', long_name='Q')
 
    # the tests
    self.assertEqual(sg.utilsg.id_in([a1,a2,a3,a4],b2  ) , True)
    self.assertEqual(sg.utilsg.id_in([a1,a2,a3,a4],b3  ) , False)   

    self.assertEqual(sg.utilsg.id_index([a1,a2,a3,a4],b2  ) , 1)
    self.assertEqual(sg.utilsg.id_index([a1,a2,a3,a4],b3  ) , None)
    self.assertEqual(sg.utilsg.rem_equivs([a1,a2,a3]+[b2,] ),  [a1, a2, a3] )

  def test_get_att_function(self):
    """Tests sg.utilsg.get_att
    """

    # define some test class with some attributes
    class Tmp(object):
      test =0
      test2=20
      test3=30

    W = Tmp()
    self.assertEqual(sg.utilsg.get_att(W,['test','test2'] ), 0 )
    self.assertEqual(sg.utilsg.get_att(W,['test2','test'] ), 20 )
    self.assertEqual(sg.utilsg.get_att(W,['test100'] ), None )


  def test_merge_function(self):
     """
     Tests whether 2 test arrays are properly merged.
     """

     # two test arrays to merge
     A = np.array([1.,2.,3.,4.])
     B = np.array([-10.,1.5,2.5,3.5,4.5,11.])

     self.assertEqual( np.array_equal(sg.utilsg.merge(A,B), np.array([-10. ,   1. ,   1.5,   2. ,   2.5,   3. ,   3.5,   4. ,   4.5,  11. ]) ), True  )


# 3 tests for very simple function sublist in utilsg.py
# --------------
  def test_sublist(self):

    self.assertEqual(sg.utilsg.sublist(['test','hi'] ,'hi' ) , ['hi'])
     
  def test_sublist_all(self):

    self.assertEqual(sg.utilsg.sublist(['test','hi'] ,'*' ) , ['test','hi'])

  def test_sublist_none(self):

    self.assertEqual(sg.utilsg.sublist(['test','hi'] ,'ho' ) , [])


# -------------

  def test_add_alias(self):
    """
    Create some test coords to test the add_alias function in utilsg.py.

    An alias attribute is assigned, which is the same as the name attribute unless the name appears more than once.
    Two names are the same in this example, and in the created alias, the second of those two names must receive a suffix "2". 
    """
    coord1 = sg.fieldcls.Coord(name = 'test',direction ='X',value =np.array([1.,2.,3.]) , metadata = {'hi':5} )
    coord2 = sg.fieldcls.Coord(name = 'test',direction ='Y',value =np.array([1.,2.,3.,4.]), metadata = {'hi':7})

    coord3 = sg.fieldcls.Coord(name = 'test3',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':3})

    coord4 = sg.fieldcls.Coord(name = 'test4',direction ='X',value =np.array([5.,1.,2.,3.,4.]), metadata = {'hi':5})

    L = sg.utilsg.add_alias([coord1, coord2, coord3, coord4])

    # test that alias is correct (same as names, but if the same name occurs >1 times, it is numbered)
    self.assertEqual([it.alias for it in L]  , ['test', 'test2', 'test3', 'test4'] )

    # test that names remain the same
    self.assertEqual([it.name for it in L]  , ['test', 'test', 'test3', 'test4'] )

  def test_find_perm_function_equal_length_permutables(self):
    """
    Test whether the permutation between two permutable lists yields the right result.
    """
    left = ['a','b','c']
    right = ['c','a','b']

    perm = sg.utilsg.find_perm(left,right)

    self.assertEqual([left[i] for i in perm] , right)

  def test_find_perm_function_non_equal_length(self):
    """
    Test whether the permutation between two non-permutable lists yields the right result.
    """
    left = ['a','b','c']
    right = ['c','a']

    perm = sg.utilsg.find_perm(left,right)

    self.assertEqual(perm, None)

  def test_find_perm_function_equal_length_non_permutables(self):
    """
    Test whether the permutation between two permutable lists yields the right result.
    """

    a=sg.fieldcls.Coord('a')

    b=sg.fieldcls.Coord('b')

    c=sg.fieldcls.Coord('c')

    left = [a,b]
    right = [b,c]

    perm = sg.utilsg.find_perm(left,right)

    self.assertEqual(perm, None)



  def test_simple_glob_function_left_wildcard(self):

    self.assertEqual(sg.utilsg.simple_glob(['foo','bar'],'*oo'  ), ['foo'] )

  def test_simple_glob_function_right_wildcard(self):

    self.assertEqual(sg.utilsg.simple_glob(['foo','bar'],'oo*'  ), ['foo'] )

  def test_simple_glob_function_right_wildcard(self):

    self.assertEqual(sg.utilsg.simple_glob(['foo','bar','vroom'],'*oo*'  ), ['foo','vroom'] )

  def test_simple_glob_function_no_wildcard(self):

    self.assertEqual(sg.utilsg.simple_glob(['foo','bar','vroom'],'oo'  ), [] )

  def test_end_of_filepath_function(self):

    self.assertEqual(sg.utilsg.end_of_filepath('/test/foo/bar'), 'bar')
    self.assertEqual(sg.utilsg.end_of_filepath('/foo/bar/'), 'bar')
    self.assertEqual(sg.utilsg.end_of_filepath('foo/bar/'), 'bar')

class TestExper(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info_dict()
    P = sg.Project(D['my_project']);
    #P.load('O_temp')
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture


  def test_load_method_non_existent_var(self):
    P = self.fixture
    E = P['DPO']
    varname = 'this_doesnt_exist'

    # attempt to load non-existent field    
    P.load(varname)
    
    self.assertEqual(len(E.vars),0)


  def test_load_method_existent_var(self):
    P = self.fixture
    E = P['DPO']
    varname = 'O_temp'

    # attempt to load non-existent field    
    P.load(varname)
    
    self.assertEqual(len(E.vars),1)

  def test_load_method_multiple_existent_var(self):
    P = self.fixture
    E = P['DPO']
    varnames = ['A_sat', 'A_slat' ]

    # attempt to load non-existent field    
    P.load(varnames)
    
    self.assertEqual(len(E.vars),2)


  def test_get_function_of_Exper_not_loaded(self):
    # try to get a Field that has not been loaded yet from the Exper object => None returned.
    E = self.fixture['DPO']
   
    self.assertEqual(E.get('O_temp'),None)

  def test_get_of_Exper(self):
    # try to get a Field that has been loaded from the Exper object => Field object.
    E = self.fixture['DPO']
    self.fixture.load('O_temp')
    self.assertEqual(str(E.get('O_temp')), 'O_temp')

  def test_delvar_method_of_Exper(self):  

    # try to delete a Field that has been loaded from the Exper object 

    E = self.fixture['DPO']
    self.fixture.load(['O_temp','O_sal','A_sat','A_shum'])
    E.delvar('O_temp')
    self.assertEqual(E.get('O_temp') is None, True)

    E.delvar(['O_sal','A_sat'])
    self.assertEqual(E.get('O_sal') is None, True)
    self.assertEqual(E.get('A_sat') is None, True)

    del E['A_shum']
    self.assertEqual(E.get('A_shum') is None, True)

# tests around coord and grid aspects of fields

class TestCoordField(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info_dict()
    P = sg.Project(D['my_project']);P.load('O_temp')
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture


  def test_field_grid_len(self):

    self.assertEqual(len(self.fixture['DPO']['O_temp'].grid),3)

  def test_field_shape(self):

    self.assertEqual(self.fixture['DPO']['O_temp'].shape,self.fixture['DPO']['O_temp'].grid.shape())

  def test_coord(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'

    self.assertEqual( latitude*(longitude*latitude) , longitude*latitude )

  def test_coord_mult2(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'

    self.assertEqual( latitude_edges*(longitude*latitude) , longitude*latitude_edges )
    

    
  def test_coord_div(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'

    self.assertEqual( (longitude*latitude)/longitude , latitude**2 )

   

  def test_coord_dual(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'

    self.assertEqual( longitude.dual, longitude_edges )
   
  def test_coord_mul_field(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'

    self.assertEqual( (longitude*self.fixture['DPO']['O_temp']).shape, (19,100) )

  def test_coord_div_field(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'

    self.assertEqual( (self.fixture['DPO']['O_temp'] / longitude).shape, (19,100) )

  def test_coord_2D_div_field(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'

    self.assertEqual( (self.fixture['DPO']['O_temp'] / (longitude*latitude ) ).shape, (19,) )


  def test_ax_mul_field(self):
    
    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    self.assertEqual( (X*self.fixture['DPO']['O_temp']  ).shape, (19, 100) )


  def test_can_I_divide_field_by_ax_shape(self):
    
    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    self.assertEqual( (self.fixture['DPO']['O_temp'] / X ).shape, (19, 100) )

  def test_can_I_divide_field_by_ax2D_shape(self):
    
    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    self.assertEqual( (self.fixture['DPO']['O_temp'] / (X*Y ) ).shape, (19,) )


  def test_avg_temp_value(self):

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    self.assertAlmostEqual( self.fixture['DPO']['O_temp']/ (X*Y*Z) ,  3.9464440090035104 , places =2)


  def test_avg_temp_value_after_regrid(self):

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    # load velocity to get the velocity grid
    self.fixture.load('O_velX')

    TEMP_regrid = self.fixture['DPO']['O_temp'].regrid(self.fixture['DPO']['O_velX'].grid)

    self.assertAlmostEqual( TEMP_regrid/ (X*Y*Z) , 4.092108709111132  , places =2)



  def test_squeezed_dims_worked_on_loading(self):

    self.assertEqual( len(self.fixture['DPO']['O_temp'].squeezed_dims) , 1   )

  def test_if_unsqueezing_adds_dims(self):

    self.assertEqual( len( (sg.unsqueeze(self.fixture['DPO']['O_temp']) ).grid ) , 4   )

  def test_if_unsqueezing_removes_squeezed_dims(self):

    self.assertEqual( len( (sg.unsqueeze(self.fixture['DPO']['O_temp']) ).squeezed_dims ) , 0   )


  def test_Gr_squeeze_method(self):

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'


    TEMP = self.fixture['DPO']['O_temp']
    G = TEMP.grid

    self.assertEqual( len(G.squeeze()[0] ) , 3  )

    G = (TEMP[Y,50]).grid
    self.assertEqual( len(G.squeeze()[0] ) , 2  )

    G = (TEMP[Z,0,X,10]).grid
    self.assertEqual( len(G.squeeze()[0] ) , 1  )


  def test_squeeze_multiple_1dim(self):


    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    TEMP = self.fixture['DPO']['O_temp']


    K=TEMP[Y,50]
    self.assertEqual( (sg.squeeze(K)).shape , (19, 100)   )

    K=TEMP[Z,0,X,10]
    self.assertEqual( (sg.squeeze(K)).shape , (100,)   )



class TestFieldBasic(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info_dict()
    P = sg.Project(D['my_project']);
    P.load(['O_temp','A_sat'])
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture

  def test_slice_NH(self):

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'       

    SAT_sliced = SAT[Y,:50]

    self.assertEqual( SAT_sliced.shape ,  (50,100)  )

  def test_slice_one_lat(self):

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'       

    SAT_sliced = SAT[Y,50]

    self.assertEqual( SAT_sliced.shape ,  (1,100)  )

  def test_slice_everything(self):
    """ Slicing with : should yield the value attribute, an ndarray
    """

    SAT = self.fixture['DPO']['A_sat']

    SAT_sliced = SAT[:]

    self.assertEqual( isinstance(SAT_sliced, np.ndarray) ,  True  )




  def test_concatenate_arg_ax_None(self):
    """
    Test the sg.concatenate function with ax argument None.
    """
    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'       

    SAT1 = SAT[Y,:40]
    SAT2 = SAT[Y,40:55]
    SAT3 = SAT[Y,55:]

    SAT_combined = sg.concatenate((SAT1,SAT2,SAT3))

    self.assertEqual( SAT_combined.shape ,  (100,100)  )



  def test_concatenate_arg_ax_not_in_grid(self):
    """
    Test the sg.concatenate function with ax argument that points in a different axis direction from the grid. This should lead to a new Coord object that is added to the result grid and that we can examine.   
    """

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'       

    SAT1 = SAT[Y,:50]
    SAT2 = SAT[Y,50:]

    # Create test Coord to concatenate along.
    W = sg.Ax('W')
 
    SAT_combined = sg.concatenate([SAT1,SAT2 ], ax = W )

    self.assertEqual( SAT_combined.shape ,  (2,50,100)  )
    # concatenate has created a new Coord:
    self.assertEqual( np.array_equal(SAT_combined.grid[0].value, np.array([0.,1.]) ),  True  )


  def test_concatenate_arg_new_coord_given(self):
    """
    Test the sg.concatenate function with new_coord argument an indpendendent Coord.    
    """

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'       

    SAT1 = SAT[Y,:50]
    SAT2 = SAT[Y,50:]

    # Create test Coord to concatenate along.
    W = sg.Ax('W')
    w = sg.Coord('w' , axis = W, direction = 'W', value = np.array([0,1]))

    SAT_combined = sg.concatenate([SAT1,SAT2 ], new_coord = w )

    self.assertEqual( SAT_combined.shape ,  (2,50,100)  )




class TestVectorField(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info_dict()
    P = sg.Project(D['my_project']);
    P.load(['O_velX','O_velY','O_temp'])
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture

  def test_slice(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'       

    U = self.fixture['DPO']['O_velX']
    V = self.fixture['DPO']['O_velY']
    TEMP = self.fixture['DPO']['O_temp']

    # just to speed up multiplication (otherwise regridding takes place): 
    TEMP.grid = U.grid

    UV = U*V

    # Did multiplication yield a 2D vectorfield?
    self.assertEqual(len(UV),2)

    # scalar field with vector component should yield Field
    self.assertEqual(isinstance(TEMP*V, sg.fieldcls.Field),True)

    # check that vcumsum and vsum propagate to the Field members of VField:
    Ucs = U.vcumsum(coord = latitude_V)

    UVcs = UV.vcumsum(coord = latitude_V)

    R1 = Ucs.value
    R2 = UVcs[0].value
    
    R1[np.isnan(R1)] = 0
    R2[np.isnan(R2)] = 0

    self.assertEqual(np.array_equal(R1,R2) ,True )


    Ucs = U.vsum( )

    UVcs = UV.vsum( )

    R1 = Ucs
    R2 = UVcs[0]
 
    self.assertEqual(R1,R2)



class TestGrid(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info_dict()
    P = sg.Project(D['my_project']);
    P.load(['O_temp','A_sat'])
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture


  def test_division(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'   

    self.assertEqual((latitude*longitude)/X,latitude**2)

  def test_inflate(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    Igr = (depth*latitude*longitude).inflate()

    self.assertEqual(Igr[0].shape, (19, 100, 100))

  def test_grid_empty_grid_equal(self):

    self.assertEqual(sg.Gr() == sg.Gr(), True)
    

  def test_grid_sliced_method(self):

# Corresponds to CASE 1a in equal length grid case in fieldcls.py source code.
    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'   

    gr1 = depth*latitude*longitude

    gr1_sliced = gr1.sliced((X,slice(1,None,None)) )

    self.assertEqual(gr1_sliced.shape(), (19, 100, 99) )

    self.assertEqual(gr1_sliced[0] is depth, True )

    # Try single slab slice:
    gr1_sliced = gr1.sliced((X,10))
    self.assertEqual(gr1_sliced.shape(), (19,100,1) )


  def test_grid_permute_function_equal_len_and_coords(self):

# Corresponds to CASE 1a in equal length grid case in fieldcls.py source code.
    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*longitude
    gr2 = longitude*depth

    # define a np array consistent with gr1
    A = np.ones( gr1.shape()  )

    # gr1(gr2) should yield a function transposing ndarrays consistent with gr1 to ndarrays consistent with gr2

    self.assertEqual((gr1(gr2)(A)).shape, gr2.shape() )

  def test_grid_permute_function_equal_len_equiv_coords_only(self):

# Corresponds to CASE 1b in equal length grid case in fieldcls.py source code.

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

   # This time, we are going to a new grid that requires interpolation (on longitude).
 
    gr1 = depth*longitude
    gr2 = longitude_V*depth

    # define a np array consistent with gr1
    A = np.ones( gr1.shape()  )

    # gr1(gr2) should yield a function transposing ndarrays consistent with gr1 to ndarrays consistent with gr2, and interpolated onto it.

    self.assertEqual((gr1(gr2)(A)).shape, gr2.shape() )

  def test_grid_permute_function_equal_len_equiv_coords_only(self):

# Corresponds to CASE 1c in equal length grid case in fieldcls.py source code.

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

   # This time, we are going to a new grid that is incompatible, leading to a None result.
 
    gr1 = depth*longitude
    gr2 = latitude*depth

    self.assertEqual(gr1(gr2), None )

  def test_gr_interpret_slices_function(self):


    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   


    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'   


   # This time, we are going to a new grid that is incompatible, leading to a None result.
 
    gr1 = depth*latitude*longitude
 
    slNone = slice(None, None, None)

    self.assertEqual(sg.interpret_slices((longitude,10),gr1) == (slNone, slNone, slice(10,11,None) )  , True )


    self.assertEqual(sg.interpret_slices((X,10),gr1) == (slNone, slNone, slice(10,11,None) )  , True )

#    self.assertEqual(sg.interpret_slices(10 , slNone, 10   )

#    self.assertEqual( sg.interpret_slices((slNone, slNone),G) ,  (slNone, slNone)   )

  def test_gr_method_expand_size(self):
    """
    Test expand method of fieldcls.py


    SAT = P['DPO']['A_sat']
    SAT.shape is (100,100)
    W=SAT.grid.expand(SAT[:],depth**2)
    W.shape is (19,100,100)
    W contains 19 identical copies (slices) of SAT[:] 

    """

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   
    
    W=SAT.grid.expand(SAT[:],depth**2)

#    W has been expanded, and the other grid (depth**2) should be appended on the left side.         

    self.assertEqual(W.shape, (19,100,100)  )

  def test_gr_method_expand_broadcast(self):
    """
    Test expand method of fieldcls.py
    """

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   
    
    W=SAT.grid.expand(SAT[:],depth**2)

#    W contains 19 identical copies (slices) of SAT[:]          

    K=W[:,50,50]
    self.assertEqual((K == K[0]).all() , True  )

  def test_call_small_gr_on_big_gr(self):

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'   

    # need to do slice test earlier.
    SAT2 = SAT[Y,:50]

    gr1 = SAT2.grid
    gr2 = depth*SAT2.grid

    A = SAT2[:]
    B = gr1(gr2)(A)

    self.assertEqual(B.shape ,  (19, 50, 100) )

  def test_call_small_gr_on_big_gr_permute(self):

    """
    corresponds to case 2a of gr class call method in fieldcls.py
    """

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'   

    # need to do slice test earlier.
    SAT2 = SAT[Y,:50]

    gr1 = SAT2.grid
    # note that this does something different for a single coord left multiplicant:
    gr2 = (depth*longitude)*SAT2.grid   

    A = SAT2[:]
    B = gr1(gr2)(A)

    self.assertEqual(B.shape ,  (19, 100, 50) )

  def test_call_small_gr_on_big_gr_permute_interp(self):

    """
    corresponds to case 2b of gr class call method in fieldcls.py
    """


    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'   

    # need to do slice test earlier.
    SAT2 = SAT[Y,:50]

    gr1 = SAT2.grid
    # note that this does something different for a single coord left multiplicant:
    gr2 = (depth*longitude_V)*SAT2.grid   

    A = SAT2[:]
    B = gr1(gr2)(A)

    self.assertEqual(B.shape ,  (19, 100, 50) )
    
  def test_call_small_gr_on_big_gr_not_equiv(self):

    """
    corresponds to case 2c of gr class call method in fieldcls.py
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   


    self.assertEqual(depth(latitude*longitude) ,  None )


  def test_gr_method_reduce_dim1vs3_len_list(self):
    """
    Test reduce method of gr class
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth**2
    gr2 = depth*latitude*longitude

    A = np.ones(gr2.shape() )

    # should have the length of len(depth)
    self.assertEqual(len(gr1.to_slices(A,gr2)) ,  19 )
    
  def test_gr_method_reduce_dim1vs3_shape_element(self):
    """
    Test reduce method of gr class
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth**2
    gr2 = depth*latitude*longitude

    A = np.ones(gr2.shape() )

    # should have the shape of latitude*longitude
    self.assertEqual( gr1.to_slices(A,gr2)[0].shape  ,  (100,100) )

  def test_gr_method_reduce_dim2vs3_len_list(self):
    """
    Test reduce method of gr class
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude
    gr2 = depth*latitude*longitude

    A = np.ones(gr2.shape() )

    # should have the length of len(depth)*len(longitude) 
    self.assertEqual(len(gr1.to_slices(A,gr2)) ,  1900 )
    
  def test_gr_method_to_slices_dim2vs3_shape_element(self):
    """
    Test to_slices method of gr class
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude
    gr2 = depth*latitude*longitude

    A = np.ones(gr2.shape() )

    # should have the shape of longitude**2
    self.assertEqual( gr1.to_slices(A,gr2)[0].shape  ,  (100,) )

  def test_Gr_method_dual(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude

    gr1_dual = gr1.dual()

    self.assertEqual(np.array_equal(gr1_dual[0].value , depth_edges.value  ) , True )

  def test_gr_method_vsum(self):
    """
    Test vsum method of gr class
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude

    # should have the shape of longitude**2
    self.assertAlmostEqual( gr1.vsum(gr1.ones() )  , 121672626836.47124 , places =2 )


  def test_gr_method__find_args_coord(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    ctypes = {'x_coord':sg.XCoord,'y_coord':sg.YCoord,'z_coord':sg.fieldcls.Coord}
  
    self.assertEqual((latitude*longitude)._find_args_coord(coord_types = ctypes) , 
    [[], [latitude]] )

    
  def test_gr_method_der_type(self):

    """
    Test der method of gr class to see whether it returns a Field
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = longitude*latitude

    # should have the shape of longitude**2
    self.assertEqual( isinstance( gr1.der(longitude,gr1.ones() ) , sg.Field )  , True )

  def test_gr_method_der_X(self):

    """
    Test der method of gr class to see whether it returns a Field
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = longitude*latitude

    W = gr1.der(longitude,gr1.ones() )

    W.value[np.isnan(W.value)]=1

    # should have the shape of longitude**2
    self.assertEqual( W.value.sum()   , 0.0 )

  def test_gr_method_der_Y(self):

    """
    Test der method of gr class to see whether it returns a Field
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude

    W = gr1.der(latitude,gr1.ones() )

    W.value[np.isnan(W.value)]=1

    # should have the shape of longitude**2
    self.assertEqual( W.value.sum()   , 19.0 )

  def test_gr_method_vol(self):

    """
    Test volume method 
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude

    W = gr1.vol()
    # should have the shape of longitude**2
    self.assertAlmostEqual( W.value.sum()   , 121672626836.47124 , places = 2 )


class TesHigherFieldFunctionality(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info_dict()
    P = sg.Project(D['my_project']);P.load('F_heat')
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture


  def test_meridional_heat_transport(self):

    P = self.fixture

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'   


    # obtain oceanic heat flux as sg field object HF from project.
    HF = P['DPO']['F_heat']
    HF2 = P['DPC']['F_heat']

    PHT = Y|(HF*X)*1e-15
    PHT2 = Y|(HF2*X)*1e-15

    self.assertEqual(PHT.shape, (100,))
    self.assertEqual(PHT2.shape, (100,))






# --------- run the classes ------------

if __name__ == '__main__':
    unittest.main()

