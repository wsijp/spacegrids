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


class TestGrid(unittest.TestCase):

  def setUp(self):
    print 'Setting up %s'%type(self).__name__
    D = sg.info(verbose = False)
    P = sg.project(D['my_project']);P.load('O_temp')
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


# CHECK WHAT THE SHAPE OF A SHOULD BE!
    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*longitude
    gr2 = longitude*depth

    A = np.ones( gr1.shape()  )

    self.assertEqual((gr2(gr1)(A)).shape, (100, 19))

  def test_grid_permute_function_equal_len_diff_coords(self):

    for c in self.fixture['DPO'].cstack:
      exec c.name + ' = c'   

    gr1 = depth*longitude
    gr2 = longitude_V*depth

    A = np.ones( gr2.shape()  )

    self.assertEqual((gr2(gr1)(A)).shape, (100, 19))





# --------- run the classes ------------

if __name__ == '__main__':
    unittest.main()

