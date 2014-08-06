""" 
"""

import os
import numpy as np


from fieldcls import *
from  expercls import *
from projectcls import *
from ops import *
from plotting import *
from _config import *
from _iosg import *
from utilsg import *

# masked arrays are done as follows:
#  S=[1.,2.,3.,np.nan,5.]
#  R=ma.masked_array(S,np.isnan(S))
#  ma.sum(R)

# Land is usually nan in fields. To obtain a field from F where land is nan as in F, but the rest is 1, divide F by itself: F/F.
# If F is on grid yt*xt, say, and dA=dy*dx is a field containing the surface areas of the grid cells (dx = xt.d(yt), dy = yt.d()), then dA*(F/F) yields the surface area of the ocean area (or the non-nan area).



