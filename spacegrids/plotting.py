#encoding:utf-8

""" Functions relating to plotting Fields.

Field objects contain more information than Numpy ndarrays. For instance, information about grids is stored alongside the actual data defined on that grid. This information can be used to automatically choose figure axes and labels etc. This is done in the plot, contour and contourf functions contained in this module. These functions take Field arguments, and call their plt counterparts by the same name. Arguments can be passed on directly to these counterparts.
"""

import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings

from _config import *
from fieldcls import *

warnings.formatwarning = warning_on_one_line

def auto_cont(m,M, num_cont, max_try = 5):
  """
  Automatically pick contours.

  Pick a quantity of around num_cont contours between integers m and M. 

  Int max_try indicates how many times to try and refine the boundaries of the interval.

  Returns:
    1D ndarray containing the desired contours: np.arange(m_new,M_new,step)    
  """

  if M < m:
    rng = auto_cont(m = -m, M = - M, num_cont = num_cont)
    return -rng

  raw_step = abs(M - m)/float(num_cont)

  # compute orders of magnitude:
  m_order = order_mag(m)
  M_order = order_mag(M)
  step_order = order_mag(raw_step)

  step = round_order(raw_step)
  m_new = round_order(m,step_order) 
  M_new = round_order(M, step_order) + step

  # sometimes grey areas remain in the plot beyond the m and M margins. This is to get rid of that by adding a few more contours on either side.
  tries = 0
  while (m-step < m_new) and (tries < max_try):
    m_new = m_new - step
    tries += 1

  tries = 0
  while (M+step > M_new) and (tries < max_try):
    M_new = M_new + step
    tries += 1

  return np.arange(m_new,M_new,step)      

def _prep_axes(fld, minus_z=True,xl=None,yl=None,xscale = 1.,yscale = 1.,ax_units=True , grid = None):
  """
  Prepare axes names etc for plotting.

  Supplies data and label information ready to use in plotting.

  Args:
    fld: (Field) to be plotted
    minus_z: (Boolean) display negative numbers on vertical scale (Z) if True
    xl: (None, list/ tuple of length 2) plot limits in the x-direction
    yl: (None, list/ tuple of length 2) plot limits in the y-direction
    xscale: (float) scale in x direction
    yscale: (float) scale in y direction
    ax_units: (Boolean) display units on axes if True
    grid: (None or Gr) grid to draw information from.

  Returns:
    body,mbody,M,m,X,Y,xlbl,ylbl,xscale,yscale
  """

  xlbl=''
  ylbl = '' 
 
  if grid is None:
    grid = fld.grid


  X, Y, x_name, y_name = grid[1], grid[0],grid[1].name, grid[0].name


  if minus_z:

    if hasattr(grid[0],'axis'):
      if grid[0].axis.name == 'Z':
        yscale *= -1
      elif grid[1].axis.name == 'Z':
        xscale *= -1

# determine full label display names (e.g. 'Longitude')
  if hasattr(grid[0],'axis'):
  
    ylbl = grid[0].axis.long_name
  else:
    ylbl = grid[0].name
 
  if hasattr(grid[1],'axis'):
    if hasattr(grid[1].axis, 'long_name'):
      if isinstance(grid[1].axis.long_name,str) or isinstance(grid[1].axis.long_name,unicode):
        xlbl = grid[1].axis.long_name
      else:
        warnings.warn('Label not added: not a string.')  
        xlbl = ''
    else:
      warnings.warn('Label not added: not available.')     
      xlbl = ''
  else:
    xlbl = grid[1].name

  if ax_units:
    if hasattr(grid[0],'units'):
      if isinstance(grid[0].units,str) or isinstance(grid[0].units,unicode):
        ylbl += ' (' +grid[0].units +')'
      else:
        warnings.warn('Label not added: not a string.')
        ylbl = ''
    else:
      warnings.warn('Label not added: not available.')
      ylbl = ''

    if hasattr(grid[1],'units'):
      if isinstance(grid[1].units,str) or isinstance(grid[1].units,unicode):
        xlbl += ' (' +grid[1].units +')'
      else:
        warnings.warn('Label not added: not a string.')
        xlbl = ''
    else:
      warnings.warn('Label not added: not available.')
      xlbl = ''


  if xl is None:
    xl = [0,len(X)]

  if yl is None:
    yl = [0,len(Y)]

  X = X[xl[0]:xl[1]]
  Y = Y[yl[0]:yl[1]]

  body = fld[:]
  mbody = ma.array(body)[yl[0]:yl[1],xl[0]:xl[1]]

  M = np.nanmax(body)
  m = np.nanmin(body)


  return body,mbody,M,m,X,Y,xlbl,ylbl,xscale,yscale



def quiver(vfld, showland=True, xlabel = True,ylabel = True, minus_z=True,xl=None,yl=None,xscale = 1.,yscale = 1.,ax_units=True,ax=None,greyshade = '0.65',kmt = None, **kwargs):
  """
  Function that calls Matplotlib quiver and passes on arguments, but takes a Field as argument (instead of an ndarray)
  
  Assuming that z vs x and z vs y quiver plots are rare: expects x-y plots.  

  Args: 
    vfld: (VField) to be plotted
    showland: (Boolean) show land (nan) is True (not working yet)
    xlabel: (Boolean) show x label if True
    ylabel: (Boolean) show y label if True  
    minus_z: (Boolean) display negative numbers on vertical scale (Z) if True
    xl: (None, list/ tuple of length 2) plot limits in the x-direction
    yl: (None, list/ tuple of length 2) plot limits in the y-direction
    xscale: (float) scale in x direction
    yscale: (float) scale in y direction
    ax_units: (Boolean) display units on axes if True
    ax: (None or figure axis handle) set to plt.gca() if argument is None
    greyshade: (str) shade of grey to use for land
    ***kwargs: normal kw arguments that can be passed on to quiver.
  """

  greyshade = float(greyshade)

  if len(vfld.direction()) == 2:
    if len(vfld[0].grid) == 2:

 # obtain prepared arrays and names from Field object. mbody is a masked array containing the Field data. mbody will be used in plotting. Some of the output will not be needed.
          
      if vfld.direction()[0].name == 'X':

        body_x,mbody_x,M,m,X,Y,xlbl,ylbl,xscale,yscale = _prep_axes(fld=vfld[0],  minus_z=minus_z, xl=xl,yl=yl,xscale = xscale,yscale = yscale,ax_units=ax_units)

        body_y,mbody_y,M,m,X,Y,xlbl,ylbl,xscale,yscale = _prep_axes(fld=vfld[1],  minus_z=minus_z, xl=xl,yl=yl,xscale = xscale,yscale = yscale,ax_units=ax_units)
        
      elif vfld.direction()[0].name == 'Y':

        body_y,mbody_y,M,m,X,Y,xlbl,ylbl,xscale,yscale = _prep_axes(fld=vfld[0],   minus_z=minus_z, xl=xl,yl=yl,xscale = xscale,yscale = yscale,ax_units=ax_units)

        body_x,mbody_x,M,m,X,Y,xlbl,ylbl,xscale,yscale = _prep_axes(fld=vfld[1],   minus_z=minus_z, xl=xl,yl=yl,xscale = xscale,yscale = yscale,ax_units=ax_units)


      if showland:
        msk = copy.deepcopy(body_x)
        msk[np.isnan(msk) == False] = 1.
        msk[np.isnan(msk) == True] = 0.
    
        BW_pair = ((0.,0.,0.),(1., 1.,1.))
       
        cdict = {'red':BW_pair,'green':BW_pair,'blue':BW_pair}

        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict)
        plt.pcolormesh(xscale*X,yscale*Y,msk,cmap = my_cmap)

      Q = plt.quiver(xscale*X,yscale*Y,mbody_x,mbody_y, **kwargs)    

 
      if xlabel:
        plt.xlabel(xlbl)
      if ylabel:
        plt.ylabel(ylbl) 

      if ax is not None:
        ax = plt.gca()


      return Q

    else:
      warnings.warn('Define quiver input on 2D grid.',RuntimeWarning)


def _scale_prep_deco(func):
  """Decorator used for contour and contourf
  """

  def plot_wrapper(fld, num_cont =15, xlabel = True,ylabel = True, minus_z=True,xl=None,yl=None,xscale = 1.,yscale = 1.,ax_units=True, num_xticks = 6, num_yticks = 6, greyshade = '0.65', showland = False,grid = None,*args, **kwargs):
    """  
    Returned function for _scale_prep_deco decorator.


    Wraps contour and contourf, who take slightly different args: arguments specific to only one function:
  
    contourf
    greyshade = '0.65'

    contour
    showland = False


    Args:
      fld: (Field) to be plotted
      num_cont: (int) number of contours to be plotted
      xlabel: (Boolean) show x label if True
      ylabel: (Boolean) show y label if True  
      minus_z: (Boolean) display negative numbers on vertical scale (Z) if True
      xl: (None, list/ tuple of length 2) plot limits in the x-direction
      yl: (None, list/ tuple of length 2) plot limits in the y-direction
      xscale: (float) scale in x direction
      yscale: (float) scale in y direction
      ax_units: (Boolean) display units on axes if True
      num_xticks: (None or int) number of ticks on x-axis, or disabled
      num_yticks: (None or int) number of ticks on y-axis, or disabled
      greyshade: (str) shade of grey to use for land
      showland: (Boolean) show land (nan) is True
      grid: (None or Gr) grid to draw information from.
      *args,**kwargs: normal kw arguments that can be passed on to contour(f).
    """

  # obtain prepared arrays and names from Field object. mbody is a masked array containing the Field data. mbody will be used in plotting.
  # M, m are the max and min of the data.


#    body,mbody,M,m,X,Y,xlbl,ylbl,xscale,yscale = _prep_axes({k:kwargs[k] for k in ['fld', 'num_cont', 'xlabel' ,'ylabel', 'minus_z', 'xl','yl','xscale','yscale','ax_units']})

    body,mbody,M,m,X,Y,xlbl,ylbl,xscale,yscale = _prep_axes(fld=fld,   minus_z=minus_z, xl=xl,yl=yl,xscale = xscale,yscale = yscale,ax_units=ax_units, grid = grid)

  # Use of X, Y confusing here, as they refer to coord objects, whereas X,Y,... usually refer to ax objects. Change in later version

#    print num_cont 
    if (num_cont > 0) and not('levels' in kwargs):
      levels =  auto_cont(m,M,num_cont)

    else:
      levels = kwargs['levels']
      del kwargs['levels']


    X_scaled = xscale*X
    Y_scaled = yscale*Y

    cset = func(X_scaled,Y_scaled,mbody, levels = levels, num_cont =num_cont, xlabel = xlabel,ylabel = ylabel, minus_z=minus_z,xl=xl,yl=yl,xscale = xscale ,yscale = yscale ,ax_units=ax_units, num_xticks = num_xticks, num_yticks = num_yticks, greyshade = greyshade, showland = showland, grid = grid,**kwargs)
    
    if xlabel:
      plt.xlabel(xlbl)
    if ylabel:
      plt.ylabel(ylbl) 

    if num_xticks:
      conts = [e for e in auto_cont(X_scaled[0],X_scaled[-1],num_xticks) if e > cset.ax.get_xlim()[0] and e < cset.ax.get_xlim()[1]   ]
   
      plt.xticks(conts)
  
    if num_yticks:
   
      conts = [e for e in auto_cont(Y_scaled[0],Y_scaled[-1],num_yticks) if e > cset.ax.get_ylim()[0] and e < cset.ax.get_ylim()[1]   ]
   
      plt.yticks(conts)

    return cset

  return plot_wrapper


@_scale_prep_deco
def contourf(X_scaled,Y_scaled,mbody, levels = None, num_cont =15, xlabel = True,ylabel = True, minus_z=True,xl=None,yl=None,xscale = 1.,yscale = 1.,ax_units=True, num_xticks = 6, num_yticks = 6, greyshade = '0.65', showland = False,grid = None,*args, **kwargs):
  """
  Contourf function with Field attributes .

  Before decoration: takes axes and numpy array argument.
  After decoration: a function that calls Matplotlib contourf and passes on arguments, but takes a Field as argument (instead of a numpy array)
  

    Args: (for decorated version)
      fld: (Field) to be plotted
      num_cont: (int) number of contours to be plotted
      xlabel: (Boolean) show x label if True
      ylabel: (Boolean) show y label if True  
      minus_z: (Boolean) display negative numbers on vertical scale (Z) if True
      xl: (None, list/ tuple of length 2) plot limits in the x-direction
      yl: (None, list/ tuple of length 2) plot limits in the y-direction
      xscale: (float) scale in x direction
      yscale: (float) scale in y direction
      ax_units: (Boolean) display units on axes if True
      num_xticks: (None or int) number of ticks on x-axis, or disabled
      num_yticks: (None or int) number of ticks on y-axis, or disabled
      greyshade: (str) shade of grey to use for land
      showland: (Boolean) show land (nan) is True
      grid: (None or Gr) grid to draw information from.
      *args,**kwargs: normal kw arguments that can be passed on to contour(f).
  """


  cmap = plt.cm.jet
  cmap.set_bad('w',1.)

  cset = plt.contourf(X_scaled,Y_scaled,mbody, levels = levels, **kwargs) 

# create the patch and place it in the back of countourf (zorder!)
  if greyshade != '':
    ax = plt.gca()
  
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xy = (xmin,ymin)
    pwidth = xmax - xmin
    pheight = ymax - ymin

    p = mpl.patches.Rectangle(xy, pwidth, pheight, fill=1,color=greyshade, zorder=-10)
    ax.add_patch(p)
  
  return cset


@_scale_prep_deco
def contour(X_scaled,Y_scaled,mbody, levels = None, num_cont =15, xlabel = True,ylabel = True, minus_z=True,xl=None,yl=None,xscale = 1.,yscale = 1.,ax_units=True, num_xticks = 6, num_yticks = 6, greyshade = '0.65', showland = False,grid = None,*args, **kwargs):
  """
  Contour function with Field attributes .

  Before decoration: takes axes and numpy array argument.
  After decoration: a function that calls Matplotlib contourf and passes on arguments, but takes a Field as argument (instead of a numpy array)
  

    Args: (for decorated version)
      fld: (Field) to be plotted
      num_cont: (int) number of contours to be plotted
      xlabel: (Boolean) show x label if True
      ylabel: (Boolean) show y label if True  
      minus_z: (Boolean) display negative numbers on vertical scale (Z) if True
      xl: (None, list/ tuple of length 2) plot limits in the x-direction
      yl: (None, list/ tuple of length 2) plot limits in the y-direction
      xscale: (float) scale in x direction
      yscale: (float) scale in y direction
      ax_units: (Boolean) display units on axes if True
      num_xticks: (None or int) number of ticks on x-axis, or disabled
      num_yticks: (None or int) number of ticks on y-axis, or disabled
      greyshade: (str) shade of grey to use for land
      showland: (Boolean) show land (nan) is True
      grid: (None or Gr) grid to draw information from.
      *args,**kwargs: normal kw arguments that can be passed on to contour(f).  
  """

  if showland:
    msk = copy.deepcopy(mbody)
    msk[np.isnan(msk) == False] = 1.
    msk[np.isnan(msk) == True] = 0.
    
    BW_pair = ((0.,0.,0.),(1., 1.,1.))
       
    cdict = {'red':BW_pair,'green':BW_pair,'blue':BW_pair}

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict)
    plt.pcolormesh(X_scaled,Y_scaled,msk,cmap = my_cmap)
 
  cmap = plt.cm.jet
  cmap.set_bad('w',1.)

  cset = plt.contour(X_scaled,Y_scaled,mbody, levels = levels, **kwargs) 

 
  return cset


def plot(fld0 = None,fld1=None, minus_z=True,xlbl='',ylbl='', grid = None,start_zero=False, **kwargs):
   """
   Function that calls Matplotlib plot, but takes fields as arguments (instead of numpy arrays).

   Adds additional plot elements based on iformation from Field.

   Args:
     fld0: (Field) to be plotted if Field fld1 is None, otherwise x-coord with fld1 as y-coord
     minus_z: (Boolean) z-axis points downard if true
     xlbl, ylbl: (str) override Field labels if not ''
     grid: (None or Gr) replaces Field grid if not None
     start_zero: (Boolean) if True, x-axis starts at 0  
     **kwargs: additional args to pass to plot function.
   """
 
   if fld1 is None:
     # In this case, no x-axis is given, using Field x-axis or explicitly specified grid


     if grid is None:
       crd = fld0.grid[0]
     else:
       crd = grid[0]

     if start_zero is True:
       # Start x-axis at 0
       crd = crd.start_zero()

     if minus_z: 
       # z-axis points downard
       if hasattr(crd,'axis'):
         if crd.axis.name == 'Z':
           xax = fld0
           yax = -crd
         else:
           xax = crd
           yax = fld0
       else:
         xax = crd
         yax = fld0
     else:
       xax = crd
       yax = fld0
            
   else:
     xax = fld0
     yax = fld1

   if not(xlbl == 'No'):
     if xlbl == '':
       if hasattr(xax,'axis'):
         if hasattr(xax.axis, 'long_name'):
           if isinstance(xax.axis.long_name, str) or isinstance(xax.axis.long_name, unicode):
             xlbl = xax.axis.long_name
           else:
             xlbl = ''
             warnings.warn('Label not added: not a string.')

       elif isinstance(xax,Field):
         xlbl = xax.name
       if hasattr(xax,'units'):
         if isinstance(xax.units,str) or isinstance(xax.units,unicode):
           xlbl += ' ('+ xax.units + ')'
         else:
           warnings.warn('Label units not added: not a string.')


   if not(ylbl == 'No'):
     if ylbl == '':
       if hasattr(yax,'axis'):
        
         ylbl = yax.axis.long_name
       elif isinstance(yax,Field):
         if isinstance(yax.name,str) or isinstance(yax.name,unicode):
           ylbl = yax.name
         else:
           warnings.warn('Label not added. Not a string.')
       if hasattr(yax,'units'):
         if isinstance(yax.units,str) or isinstance(yax.units,unicode):
           ylbl += ' ('+ yax.units + ')'
         else:
           warnings.warn('Label units not added. Not a string.')

   else:
     ylbl = ''


   p = plt.plot(xax[:],yax[:], **kwargs)
   plt.xlabel(xlbl)
   plt.ylabel(ylbl)
 
   return p






def treat_kmt(kmt):
  """Returns a finer (interpolated) version of ndarray argument kmt.

  Used in plotting bathymetry.
  """

  sh = kmt.shape

  grx = np.array([ [x for x in range(sh[1])] for y in range(sh[0]) ]).ravel()
  gry = np.array([ [y for x in range(sh[1])] for y in range(sh[0]) ]).ravel()

  grxi = np.array([ [x/5 for x in range(5*sh[1])] for y in range(5*sh[0]) ]).ravel()
  gryi = np.array([ [y/5 for x in range(5*sh[1])] for y in range(5*sh[0]) ]).ravel()


  kmti = griddata(np.array([grx,gry]).transpose(),kmt.ravel(),np.array([grxi,gryi]).transpose(),method = 'nearest')


  kmti = kmti.reshape((5*sh[0],5*sh[1]))

  return kmti


def add_kmt(kmt):
  """Plot 1 contour of ndarray argument kmt.

  Used in plotting bathymetry. Calls treat_kmt.
  """
  kmti = treat_kmt(kmt)
  sh = kmti.shape

  plt.contour(np.linspace(0,360,sh[1]),np.linspace(-90,90,sh[0]) , kmti,levels = [0.5],interpolation='nearest',colors='k',linewidths=1.5);





