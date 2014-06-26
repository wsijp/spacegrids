import warnings
import numpy as np
import os

# choose from netcdf4, scientificio, scipyio
cdf_lib = 'netcdf4'
#cdf_lib = 'scipyio'

#use_scientificio = False

# Set the path here. This path will be used to find your experiments.

if os.getenv("LOC") == "standard":
# For now, the same behaviour as other locations
  home_path = os.environ['HOME']
else:
  home_path = os.environ['HOME']

# Name of text file designating project directories and containing project name
projnickfile = "projname"

# flag for field multiplication. Multiplication of directional fields leads to vector fields if true.
strict_vector = True

# name of subdirectory containing ocean masks. Path is relative to project directories and experiment directories.

# Names of axis direction names, to be used in figure labels and such:
# Corresponds to axis property of coord object.
ax_disp_names = {'X':'Longitude','Y':'Latitude','Z':'Depth','T':'Time'}

# global node stack for search function to find out if coord objects are equivalent:
nodes_done = []

mask_dir = '/masks/'
# these hard coded axes names (to indicate T-grid) are only used for masks
tgrid_names = {'UVic2.8_t':('yt','xt'),'UVic2.9_t':('latitude','longitude')}

# keywords by which coordinate grid type (e.g. tracer grid) can be identified/ guessed from description (long_name attribute) in netcdf file:
grid_type_names = {'ts_grid':['t grid','TS grid'],'uv_grid':['u grid','UV grid']}

# globs by which netcdf files are recognised.
cdf_file_extensions = ['*.nc','*.cdf']



# flag determining whether scripts are verbose. Default is not.

verbosity = 0

# Some useful constants:

R   = 6370.0e3;
g = 9.806
zero_kelvin = -273.15
omega  = np.pi/43082.0
angleStep = np.pi * 2 / 360
month_len=(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


#  ---------- names by which Netcdf files are searched --------------

# keywords by which coordinate spatial/ time direction (e.g. x-direction) can be identified/ guessed from description (long_name attribute) in netcdf file.
# strings starting with ! indicate keywords that should NOT appear in the descriptions. These denied keywords have precedence over the allowed keywords.


coord_dir_names = {'x_dir_names' : ['eastward','Eastward','zonal','Zonal','longitude','Longitude'], 'y_dir_names' : ['northward','Northward','meridional','Meridional','latitude','Latitude'], 'z_dir_names' : ['upward','Upward','vertical','Vertical','depth','Depth','pressure','!surface','level','Level'], 't_dir_names' : ['T','time','Time','year','Year']}

edge_names = ['edges','bounds']

fval_names = ['missing_value','fill_value','_FillValue']


# Formatting function used for warnings.

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s:%s \n' % (filename, lineno, category.__name__, message)









