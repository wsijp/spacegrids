#encoding:utf-8

""" unit tests.
Data required for these tests is the example project my_project mentioned in the tutorial.

"""

import unittest
import numpy as np
import spacegrids as sg


class TestInfo(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__

    self.fixture = sg.info(verbose = False)

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture

  def test_type(self):
    
    self.assertEqual(type(self.fixture),dict)


  def test_type2(self):
    D = self.fixture
    if len(D) > 0:
      self.assertEqual(type(D.keys()[0]),str)

  
class TestCoordField(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info(verbose = False)
    P = sg.project(D['my_project']);P.load('O_temp')
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture


  def test_field_grid_len(self):

    self.assertEqual(len(self.fixture['DPO']['O_temp'].gr),3)

  def test_field_shape(self):

    self.assertEqual(self.fixture['DPO']['O_temp'].shape,self.fixture['DPO']['O_temp'].gr.shape())

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


  def test_div_ax(self):
    
    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    self.assertEqual( (self.fixture['DPO']['O_temp'] / X ).shape, (19, 100) )

  def test_2D_div_ax(self):
    
    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    self.assertEqual( (self.fixture['DPO']['O_temp'] / (X*Y ) ).shape, (19,) )


  def test_avg_temp_value(self):

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'

    self.assertAlmostEqual( self.fixture['DPO']['O_temp']/ (X*Y*Z) ,  3.9464440090035104 )





  def test_squeezed_dims_worked_on_loading(self):

    self.assertEqual( len(self.fixture['DPO']['O_temp'].squeezed_dims) , 1   )

  def test_if_unsqueezing_adds_dims(self):

    self.assertEqual( len( (sg.unsqueeze(self.fixture['DPO']['O_temp']) ).gr ) , 4   )

  def test_if_unsqueezing_removes_squeezed_dims(self):

    self.assertEqual( len( (sg.unsqueeze(self.fixture['DPO']['O_temp']) ).squeezed_dims ) , 0   )


class TestFieldBasic(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info(verbose = False)
    P = sg.project(D['my_project']);
    P.load(['O_temp','A_sat'])
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture

  def test_slice(self):

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'       

    SAT_sliced = SAT[Y,:50]

    self.assertEqual( SAT_sliced.shape ,  (50,100)  )

  def test_cat(self):

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].axes:
      exec c.name + ' = c'       

    SAT1 = SAT[Y,:40]
    SAT2 = SAT[Y,40:55]
    SAT3 = SAT[Y,55:]

    SAT_combined = sg.concatenate((SAT1,SAT2,SAT3))

    self.assertEqual( SAT_combined.shape ,  (100,100)  )




class TestGrid(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info(verbose = False)
    P = sg.project(D['my_project']);
    P.load(['O_temp','A_sat'])
    self.fixture = P

  def tearDown(self):
    print 'Tearing down %s'%type(self).__name__
    del self.fixture



  def test_inflate(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    Igr = (depth*latitude*longitude).inflate()

    self.assertEqual(Igr[0].shape, (19, 100, 100))


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


  def test_gr_method_expand_size(self):
    """
    Test expand method of fieldcls.py


    SAT = P['DPO']['A_sat']
    SAT.shape is (100,100)
    W=SAT.gr.expand(SAT[:],depth**2)
    W.shape is (19,100,100)
    W contains 19 identical copies (slices) of SAT[:] 

    """

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   
    
    W=SAT.gr.expand(SAT[:],depth**2)

#    W has been expanded, and the other grid (depth**2) should be appended on the left side.         

    self.assertEqual(W.shape, (19,100,100)  )

  def test_gr_method_expand_broadcast(self):
    """
    Test expand method of fieldcls.py
    """

    SAT = self.fixture['DPO']['A_sat']

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   
    
    W=SAT.gr.expand(SAT[:],depth**2)

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

    gr1 = SAT2.gr
    gr2 = depth*SAT2.gr

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

    gr1 = SAT2.gr
    # note that this does something different for a single coord left multiplicant:
    gr2 = (depth*longitude)*SAT2.gr   

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

    gr1 = SAT2.gr
    # note that this does something different for a single coord left multiplicant:
    gr2 = (depth*longitude_V)*SAT2.gr   

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

  def test_gr_method_vsum(self):

    """
    Test vsum method of gr class
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude

    # should have the shape of longitude**2
    self.assertAlmostEqual( gr1.vsum(gr1.ones() )  , 121672626836.47124 )


  def test_gr_method_find_args_coord(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    ctypes = {'x_coord':sg.xcoord,'y_coord':sg.ycoord,'z_coord':sg.coord}
  
    self.assertEqual((latitude*longitude).find_args_coord(coord_types = ctypes) , 
    [[], [latitude]] )

    
  def test_gr_method_der_type(self):

    """
    Test der method of gr class to see whether it returns a field
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = longitude*latitude

    # should have the shape of longitude**2
    self.assertEqual( isinstance( gr1.der(longitude,gr1.ones() ) , sg.field )  , True )

  def test_gr_method_der_X(self):

    """
    Test der method of gr class to see whether it returns a field
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
    Test der method of gr class to see whether it returns a field
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
    Test der method of gr class to see whether it returns a field
    """

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*latitude

    W = gr1.vol()
    # should have the shape of longitude**2
    self.assertAlmostEqual( W.value.sum()   , 121672626836.47124 )



# --------- run the classes ------------

if __name__ == '__main__':
    unittest.main()

