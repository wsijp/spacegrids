# -*- coding: utf-8 -*-

import imp 

#import spacegrids as sg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#sg = imp.load_source('spacegrids', '/path/PACKAGES/spacegrids/__init__.py')

import spacegrids as sg

# <headingcell level=3>

# Overview of the available data

# <markdowncell>

# The info function returns a dictionary of all the available project names and their paths. Here, we have one project (named ```my_project```). Make sure there is at least one project listed when following these examples. 

# <codecell>

print 'module ', sg


D = sg.info()

# <codecell>

D

# <markdowncell>

# To start a project instance, named ```P```, it is easy to use:

# <codecell>

P = sg.project(D['my_project'])

# <markdowncell>

# Ignore the harmless warnings. They relate to the interpretation of the axes and dimension names found in the Netcdf files. Upon project initialization, meta data associated with the grids, dimensions and axis directions is collected from the data files and interpreted. The project is now initialized. Subdirectories of the ```test_project``` directory correspond to (numerical) experiments and also ```exper``` objects. If Netcdf files had been present in the ```test_project``` directory, they would have been interpreted as experiments. In this case, an ```exper``` object (experiments) corresponds not to subdirectory, but a single Netcdf file. This latter arrangement corresponds to the structure of a project saved with ```P.write()```, where all experiments are written to experiment files. If a path to a file (not project directory) is provided to ```sg.project```, a project is initialised assuming all project is contained in that single file.

# <headingcell level=1>

# The field class

# <headingcell level=3>

# Introduction to fields

# <markdowncell>

# To look at an example of fields (e.g. temperature or humidity), it is best to load one into the project.

# <codecell>

P.load('O_temp')

# <markdowncell>

# 
# 
# The ocean temperature field is now loaded into memory and associated with the project P. How do you know what variable name to pick? Inspection on the Bash command line of the Netcdf file with ```ncdump -h``` usually provides this information. Alternatively, experiments have an ls function: ```P['DPO'].ls()``` would yield a long list of available variables to be loaded as fields.

# <markdowncell>

# Notice that ```O_temp``` has not been loaded for Lev. This is because the Netcdf file inside the ```Lev``` directory follows different naming conventions. Inspection of the Netcdf file in the ```Lev``` directory reveals that temperature is referred to as ```otemp```. To load it, we use:

# <codecell>

P.load('otemp')

# <markdowncell>

# Now all three temperature fields are loaded. We access the temperature fields of the ```DPO``` experiment and the observations ```Lev``` and test resulting field instances:

# <codecell>

TEMP = P['DPO']['O_temp'] 
TEMP2 = P['Lev']['otemp']

TEMP # a field

# <markdowncell>

# The shape of a field corresponds to the shape of the Numpy array it contains and is consistent with the shape of its grid:

# <codecell>

TEMP.shape # show the shape of the numpy array containing the data. This is a 3D field

# <codecell>

TEMP.gr # show the grid on which the temp data is defined

# <markdowncell>

# The names of the grid components depth, latitude and longitude correspond to the names used in the Netcdf file. For instance, there are different naming conventions for the Lev data:

# <codecell>

TEMP2.gr

# <markdowncell>

# Fields can be saved to Netcdf again, along with all their metadata and coordinates via their write method, e.g. ```TEMP2.write()```. To specify the path for instance: ```TEMP2.write(path = '/home/me/')```. Experiment and coord objects also have save methods. In the case of an experiment, all loaded fields are saved to one file. An entire project can be saved with ```P.save()```, yielding a project file structure (directory containing a projname text file and Netcdf files corresponding to experiments). Generally, it is good to provide a ```name``` argument different from the project name, and safeguards exist to avoid overwriting the original project on disk, but no safeguards exist yet to check for multiple instances of the same project name (defined in the projname file): for now, this is up to the user to check. Finally, an arbitrary list of fields can be saved with ```sg.nugget```.

# <markdowncell>

# An overview of the experiments and loaded fields corresponding to the project are shown by:

# <codecell>

P

# <markdowncell>

# And for an experiment:

# <codecell>

P['DPO']

# <headingcell level=4>

# *** A note of warning: whenever the sg module is reloaded (e.g. with reload(sg)), all the above steps and the below steps need to be done again in sequence. Stale objects will lead to strange errors! ***

# <headingcell level=3>

# Slicing, plotting and interpolation

# <markdowncell>

# Objects representing general coordinate directions (X,Y axes etc) have been constructed during initialization of project P, and can be accessed as follows:

# <codecell>

P['DPO'].axes

# <markdowncell>

# Bring these axes into the namespace under their own names:

# <codecell>

for c in P['DPO'].axes:
  exec c.name + ' = c'

# <markdowncell>

# So that we can reference them:

# <codecell>

X

# <markdowncell>

# Field objects allow slicing by reference to axis name and using ndarray index values. Here, we use the axis object, followed by a slice object (e.g. ```TEMP[X,5]``` or ```TEMP[Z,:]``` or ```TEMP[Y,55:]``` or ```TEMP[X,slice(55,None,None))]```. To plot the surface temperature using sg version of contourf:

# <codecell>

sg.contourf(TEMP[Z,0])
clb = plt.colorbar()
clb.set_label(TEMP.units)
plt.show()

# <markdowncell>

# The sg module automatically selects sensible contour levels and x,y ticks. Versions older than 1.1.4 lack this functionality. The number of levels can be set using argument ```num_cont```. Alternatively, the Matplotlib levels arguments can be passed. Labels and orientation are extracted from the grid metadata. The field grid is used to extract plotting information, unless another grid is specified in the grid argument.

# <codecell>

sg.contourf(TEMP[X,50]) # a section at index 50.
clb = plt.colorbar()
clb.set_label(TEMP.units)
plt.show()

# <markdowncell>

# Slicing can also be done via indices using the normal ndarray notation:

# <codecell>

(TEMP[0,:,:]).shape == (TEMP[Z,0]).shape

# <markdowncell>

# A note of warning on stale objects. Slicing using ax objects is a common case where stale objects after a module reload lead to strange errors. If we were to reload the module in the above case, and bring the ```X```,```Y```,```Z```,```T``` ```ax``` objects into the namespace again, but use the now stale old TEMP field, the above slicing with ```Z``` will lead to a perplexing error. All objects will need to be refreshed in this case, which means redoing all the previous steps. In most common cases a module reload is not required, and this type of confusion arises mostly in module code development situations.  

# <markdowncell>

# Means are calculated easily by division.

# <codecell>

# Determine zonal mean mT of TEMP.

mTEMP = TEMP/X
mTEMP2 = TEMP2/X

sg.contourf(mTEMP,levels = range(-2,34,2)) # plot zonal mean
clb =plt.colorbar(ticks = range(0,34,4))
clb.set_label(TEMP.units)
plt.show()
# <codecell>

mTEMP.gr  # the dimensionality of the zonal mean is less

# <markdowncell>

# A calculated field such as ```mTEMP``` can be inserted into the project using the ```insert``` method. This is generally not needed, but can be useful when the project or experiment is later written to disk.

# <codecell>

P['DPO'].insert(mTEMP,name = 'mtemp')

# <markdowncell>

# Now, ```temp``` is part of the DPO experiment:

# <codecell>

P['DPO']

# <markdowncell>

# To see the effect of opening Drake Passage in the climate model: 

# <codecell>

sg.contourf(P['DPO']['O_temp'][Z,0]-P['DPC']['O_temp'][Z,0],cmap = mpl.cm.RdYlBu_r,levels = range(-8,9,1))
clb = plt.colorbar(ticks = range(-8,9,2))
clb.set_label(TEMP.units)
plt.show()
# <markdowncell>

# Data is generally defined on different grids, so interpolation is needed. In our example, the observations are defined on a finer grid than the model data:

# <codecell>

TEMP2.shape

# <markdowncell>

# To overlay a plot of the model geography on atmospheric fields, we can load the ocean bathymetry field G_kmt, interpolate it on a finer grid with ```sg.finer_field``` (to avoid angled contours and show true resolution) and use contour plot over the 0.5 contour to capture the transition from land (value 0) to ocean(values 1 to 19). Here, we overlay this geographic outline over sea level air temperature:

# <codecell>

P.load(['G_kmt','A_slat'])
KMT = P['DPO']['G_kmt']
SAT = P['DPO']['A_slat']
nc = 10
sg.contourf(SAT, num_cont = nc,cmap = mpl.cm.RdYlBu_r)
cs=sg.contour(SAT,linestyles='solid', num_cont = nc, colors = 'k')
plt.clabel(cs,fmt='%1.0f')
sg.contour(sg.finer_field(KMT), colors='k',levels=[0.5,]  )
plt.show()
# <markdowncell>

# To interpolate a field ```F``` in Spacegrids, simply call ```F``` on the desired grid ```gr``` as ```F(gr)```. For instance, to compute the difference between the zonal mean of the model and the observations:

# <codecell>

dmTEMP = mTEMP(mTEMP2.gr) - mTEMP2

# <codecell>

cont_int = np.arange(-3.,6.,1.)
cmap = mpl.cm.RdYlBu_r	# Colormap red yellow blue reversed
pl1 = sg.contourf(dmTEMP,  levels = cont_int, cmap = cmap);

clb = plt.colorbar(ticks = cont_int)
clb.set_label(TEMP.units)
plt.xlim([-80.,70.])
plt.xticks([-60,-30,0,30,60])
plt.ylabel('Depth (m)')
tle = plt.title(' Zonal avg. T difference model - observations')
plt.show()
# <markdowncell>

# The model is pretty good! It is a bit too warm in the upper 1000m and a bit too cool in the deep ocean.

# <codecell>

dmTEMP.gr # the resulting grid is that of the Lev data

# <markdowncell>

# To compute the meridional streamfunction, we load the velocity in the y direction and apply the calculations in simple notation:

# <codecell>

# Give the experiment objects convenient names E, E2:
E = P['DPO'];E2 = P['DPC']
varname = 'O_velY'

# load velocity by netcdf name.
P.load(varname)

# obtain velocity as sg field object V from project.
V = E[varname]
V2 = E2[varname]

# take the meridional stream function by taking zonal grid-integral X*V 
# and the vertical primitive of the result along Z using the pipe.
PSI = Z|(X*V)*1e-6
PSI2 = Z|(X*V2)*1e-6

# <markdowncell>

# Now that we have the required fields, we can plot them. 

# <codecell>

lvls = np.arange(-100,100,4)
xlims = [-80,70]
ylims = [-4500., -0.]

# --- start figure ---

lbl = ord('a')

pan = 0

height = 2
width = 1

rows = 2
cols = 1

F = PSI
ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

# use the sg wrapper for the contour function.
cs= sg.contour(F,linestyles='solid',showland=True,levels = lvls,colors='k');

plt.clabel(cs,fmt='%1.1f');

plt.xlim(xlims)
plt.ylim(ylims)
plt.tick_params(axis='x',labelbottom='off')
plt.xlabel('')
plt.ylabel('Depth (m)')
tle = plt.title('('+ chr(lbl) +') MOC Drake Passage open ')

lbl += 1
pan += 1

F = PSI2
ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

cs= sg.contour(F,linestyles='solid',showland=True,levels = lvls,colors='k');
plt.clabel(cs,fmt='%1.1f');
plt.xlim(xlims)
plt.ylim(ylims)

# Overiding the xlabel assigned by sg:
plt.xlabel('Latitude')
tle = plt.title('('+ chr(lbl) +') MOC Drake Passage closed ')
plt.show()
# <markdowncell>

# There is strong southern sinking in the DP closed case, and a switch to northern sinking in the DP open case. That's why there was warming in the NH.

# <markdowncell>

# Field concatenation is done with the concatenate function, similar to Numpy. If we split surface air temperature ```SAT``` latitudinally into three pieces:

# <codecell>

SAT1=SAT[Y,:40];SAT2=SAT[Y,40:50];SAT3=SAT[Y,50:]

# <markdowncell>

# And re-assemble along the default axis ```ax = None```, sg will guess to use the incomplete (differing) axes:

# <codecell>

SAT_cat = sg.concatenate([SAT1,SAT2,SAT3])

# <markdowncell>

# (This yields the same result as ```sg.concatenate([SAT1,SAT2,SAT3],ax=Y)```.) Drawing the zonal mean shows that ```SAT``` has been re-assembled correctly:

# <codecell>

(SAT_cat/X).draw()
plt.show()
# <markdowncell>

# If multiple files are found in an experiment directory, and they contain the same field, sg will assume they are pieces of the same field and concatenate them upon loading. An example would be files containing different time slices of the same field (e.g. temperature).

# <markdowncell>

# Tip: field objects have a squeeze and unsqueeze method: use help to learn more about them. Fields are squeezed upon loading and unsqueezed when saved.

# <headingcell level=3>

# Grids and coord objects

# <markdowncell>

# Grid objects model the grids on which the field data is defined. They consist of tuples containing ```coord``` objects. ```Coord``` objects consist of points along an axis and are read and named from the Netcdf files. To bring these objects into the namespace under their own name: 

# <codecell>

for c in E.cstack:
  exec c.name +' = c'

# <markdowncell>

# This allows us to reference the ```coord``` objects built from the Netcdf file by the names extracted from the metadata:

# <codecell>

depth[:] # depth coordinate points

# <codecell>

latitude[:10] # latitudinal coordinate points (truncated)

# <markdowncell>

# The result of the ```__getitem__``` method here is the ndarray containing the data of the coordinate points. It is important to realize that these ```coord``` names (depth, latitude etc) have been read from the Netcdf file and have not been hard coded into Spacegrids. 

# <markdowncell>

# A shorthand for grid construction from ```coord``` objects is via multiplication:

# <codecell>

latitude*longitude

# <codecell>

latitude*longitude*depth

# <markdowncell>

# The latter 3 dimensional grid is associated with ndarray data where the first index represents latitude, the second longitude and the third depth. Field objects contain a grid object that is always congruent with the ndarray holding the field data. The grid objects contained in fields are used to always ensure data alignment, for instance under multiplication.

# <markdowncell>

# Coord objects and grid objects are not the same. To construct a one dimensional grid object: 

# <codecell>

latitude**2

# <headingcell level=3>

# Axis alignment

# <markdowncell>

# We had:

# <codecell>

TEMP.gr

# <markdowncell>

# With data in the ndarray:

# <codecell>

(TEMP[:]).shape

# <markdowncell>

# where we recognise the length of the depth coord as 19: 

# <codecell>

len(depth[:])

# <markdowncell>

# Rearrangement of the coordinates:

# <codecell>

TEMP_shuffled = TEMP(latitude*longitude*depth)
TEMP_shuffled.gr

# <codecell>

(TEMP_shuffled[:]).shape

# <markdowncell>

# Data axis alignment is automatic. For instance, multiplying temperature with velocity:

# <codecell>

F = TEMP*V

# <codecell>

F.gr

# <markdowncell>

# Multiplication alignment works regardless of the order of the coordinates:

# <codecell>

F2 = TEMP_shuffled*V
F2.gr

# <markdowncell>

# Here, ```V``` is actually defined on a different grid: the velocity grid:

# <codecell>

V.gr

# <markdowncell>

# Why then is the result of ```TEMP*V``` defined on the temperature grid? This is because multiplication automatically triggers interpolation of the right multiplicant to the grid of the left multiplicant. Vice versa:

# <codecell>

(V*TEMP).gr

# <headingcell level=3>

# Axis objects

# <markdowncell>

# The objects ```P['DPO']```, ```P['DPC']``` and ```P['Lev']``` are experiment objects, and correspond to the directories containing the Netcdf files. Different data ("experiment") files may contain data defined on different grid and use different coordinate naming conventions. A list of the coord objects associated with an experiment object is accessed from the .cstack attribute:

# <codecell>

P['Lev'].cstack

# <codecell>

P['DPO'].cstack[:4] # truncated list of coord objects

# <markdowncell>

# The depth coordinates differ for the model and the observational data:

# <codecell>

P['Lev'].cstack[-1][:] # Access the depth levels for the Lev experiment.

# <markdowncell>

# We have the depth levels for ```P['DPO']``` in the namespace:

# <codecell>

depth[:] # Access depth levels for the model data 

# <markdowncell>

# However, both coord objects point in the same direction: ```Z```. The name ("Z") of this direction has also been read from the Netcdf files. This is encapsulated in the ```ax``` class (axis).

# <codecell>

Z.__class__.__name__

# <codecell>

depth.axis

# <markdowncell>

# Multiplication of a grid with an ```ax``` object (representing an axis) picks the ```coord``` component along that axis.

# <codecell>

X*(latitude*longitude)

# <codecell>

Y*(P['Lev'].cstack[0]*P['Lev'].cstack[1])

# <markdowncell>

# We also saw the ```ax``` object being used to calculate zonal means. In the following example, we load the surface heat flux, demonstrate some operations using ```ax``` objects and compute the poleward heat transport.

# <codecell>

# load heat flux by netcdf name.
P.load('F_heat')

# obtain oceanic heat flux as sg field object HF from project.
HF = P['DPO']['F_heat']
HF2 = P['DPC']['F_heat']

# <markdowncell>

# Again, ignore the harmless warning: the Lev data contains no heat flux data. The heat flux data is 2 dimensional, as can be seen from the grid:

# <codecell>

HF.gr

# <markdowncell>

# And looks like this:

# <codecell>

pl = sg.contourf(HF,levels = range(-200,200,20),cmap = mpl.cm.RdYlBu)
clb = plt.colorbar()
clb.set_label(HF.units)
plt.show()
# <markdowncell>

# Shift fields by reference to the ```ax``` objects:

# <codecell>

HF_shifted = sg.roll(HF,coord = X,shift = 50)

# <codecell>

pl = sg.contourf(HF_shifted,levels = range(-200,200,20),cmap = mpl.cm.RdYlBu)
clb = plt.colorbar()
clb.set_label(HF.units)
plt.show()
# <codecell>

HF_shifted.gr

# <markdowncell>

# Note that the grid is updated along with the data (the zero meridian still passes through Greenwich), resulting in a different grid (```longitude_rolled``` appears in place of ```longitude```). This behaviour can be disabled.

# <markdowncell>

# Compute the zonal integral (using multiplication notation ```HF*X```) followed by the primitive in the ```Y``` direction (using the pipe | notation):

# <codecell>

# Determine the total oceanic poleward heat transport via the y-primitive and X-intergral, and scale to Petawatt.

PHT = Y|(HF*X)*1e-15
PHT2 = Y|(HF2*X)*1e-15

# <codecell>

pl1, = sg.plot(PHT,color = 'b');
pl2, = sg.plot(PHT2,color='r');

plt.grid()
plt.xlim([-80.,80.])
plt.xticks([-60,-30,0,30,60])
plt.ylabel('PW')
tle = plt.title(' Ocean PHT (PW) ')
lgnd = plt.legend([pl1,pl2],['DP open','DP closed'],loc=2)
plt.show()
# <markdowncell>

# Interesting, the Drake Passage closed configuration has more southward poleward heat transport. We could have used the ```coord``` names to compute the PHT:

# <codecell>

PHT = latitude|(HF*longitude)*1e-15

# <markdowncell>

# The axis names X,Y,... are a shorthand for these coord names, and are easier to use. 

# <headingcell level=3>

# More grid operations

# <markdowncell>

# Operations such as derivatives and integration are methods of the ```coord``` class:

# <codecell>

latitude.der # the derivative on the sphere in the latitudinal direction

# <codecell>

X.der

# <codecell>

HF_shifted2 = longitude.coord_shift(HF,50)

# <codecell>

HF_shifted2.gr

# <markdowncell>

# Do we get the same result from both shifts? This is a shorthand to access the X ```coord``` of ```HF_shifted2```:

# <codecell>

# obtain coord in X direction
X*(HF_shifted2.gr) 

# <markdowncell>

# The two shifted ```coord``` objects are different objects

# <codecell>

# compare the ndarray content of the X components of the two shifted grids
X*(HF_shifted2.gr) == X*(HF_shifted.gr)

# <markdowncell>

# But contain the same values:

# <codecell>

(X*(HF_shifted2.gr))[:] == (X*(HF_shifted.gr))[:]

# <markdowncell>

# Other derived ```coord``` objects are obtained from slicing:

# <codecell>

HF_sliced = HF_shifted[Y,:36]

# <codecell>

pl = sg.contourf(HF_sliced,levels = range(-200,200,20),cmap = mpl.cm.RdYlBu)
clb = plt.colorbar()
clb.set_label(HF.units)
tle = plt.title('Heat flux in Southern Ocean')
plt.show()
# <markdowncell>

# Similar to the roll function, slicing gives a new grid:

# <codecell>

HF3 = HF[X,:50,Y,:36]

# <codecell>

HF3.gr

# <codecell>

len(Y*HF3.gr)

# <markdowncell>

# The sum method of a ```coord``` object only sums the ```ndarray``` content of a ```field```:

# <codecell>

sTEMP = longitude.sum(TEMP)
sTEMP.gr

# <markdowncell>

# Similar to Pandas, the NaN (```np.nan```) values are not counted in these operations. Similarly, the vsum method is the sum weighted with the volumes (lengths, areas) of the grid cells:

# <codecell>

help(longitude.vsum)

# <markdowncell>

# Similar for cumsum and vcumsum.

# <markdowncell>

# Grid objects also have volume weighted sum methods. These methods are preferred over their ```coord``` counterparts as coordinate cell dimensions may depend on other coordinates (e.g. the distance between longitudinal arcs diminishes towards the poles), and grid objects take care of this. 

# <codecell>

my_grid = latitude*longitude
my_grid.vsum(TEMP)[:]

# <codecell>

my_grid.vsum(TEMP).gr

# <markdowncell>

# The volume of the domain of non NaN values inside TEMP is found from: 

# <codecell>

TEMP.gr.vsum(TEMP.ones())

# <markdowncell>

# This can also be achieved with the multiplication notation:

# <codecell>

TEMP.ones()*(X*Y*Z)

# <markdowncell>

# The ```ones``` method of fields gives a new field where all non- NaN values (of ```TEMP``` in this case) are set to 1 and the NaN values remain unchanged. Here, the one values constitute the domain of ocean water, as NaN values indicate land (rock).

# <codecell>

# slice at surface and show part of array to illustratie only 2 values
TEMP[Z,0].ones()[27:37,:7]

# <markdowncell>

# This allows a calculation of the mean ocean temperature in the model:

# <codecell>

TEMP.gr.vsum(TEMP)/TEMP.gr.vsum(TEMP.ones())

# <codecell>

# more succinctly:
TEMP.gr.mean(TEMP)

# <markdowncell>

# This can also be achieved with the division notation:

# <codecell>

TEMP/(X*Y*Z)

# <markdowncell>

# The ```ax``` class has a derivative method:

# <codecell>

# Calculate derivative in X direction:
dTEMPdX = X.der(TEMP)
dTEMPdX.gr

# <codecell>

# compute the vertical temperature profile and plot it.
vpTEMP = TEMP/(Y*X)
pl, = sg.plot(vpTEMP)
lb = plt.xlabel('degrees C');tle = plt.title('Vertical profile of average ocean temp')
plt.show()
# <codecell>

# take the vertical derivative of the vertical temperature profile and plot it.
pl, = sg.plot(Z.der(vpTEMP))
lb = plt.xlabel('degs C per m');tle = plt.title('Vertical temperature gradient')
plt.show()
# <markdowncell>

# Another notation, involving ```ax``` objects, uses the ^ (carat) symbol:

# <codecell>

# test whether nparray content is the same. Excluding the first nan value.
(Z^vpTEMP)[1:] == (Z.der(vpTEMP))[1:]

# <markdowncell>

# The volume weighted cumulative sum method ```vcumsum``` of the ```coord``` class, the primitive integral along that axis, is easily accessed with the pipe | notation and the ```ax``` class (here the instance ```Z```). In the latter case the name of the particular coord need not be known.

# <codecell>

# test whether the depth coord method vcumsum yields the same result as Z|
depth.vcumsum(vpTEMP)[:] == (Z|vpTEMP)[:] 

# <headingcell level=1>

# Operators, divergences and vector fields.

# <markdowncell>

# Spacegrids allows representations of vector fields (e.g. velocity fields). We're loading the zonal velocity to examine a 2 dimensional vector fields constructed from ```U``` and ```V```.

# <codecell>

# load the ocean velocity in the X direction.
P.load('O_velX')

# <markdowncell>

# Earlier, multiplication of two fields yielded a new field. In contrast, multiplying the two fields ```U``` and ```V``` does not result in a new field here. Instead, we get a container class holding the two fields. This container is the ```vfield``` class (vector field).

# <codecell>

# obtain horizontal velocity fields
U = P['DPC']['O_velX']
V = P['DPO']['O_velY']
# create a 2 dimensional vector field defined on 3 dimensional space
my_vector_field = U*V
my_vector_field

# <markdowncell>

# Why do these fields behave differently under multiplication? Fields have a direction attribute:

# <codecell>

U.direction

# <codecell>

V.direction

# <codecell>

TEMP.direction

# <markdowncell>

# The velocity fields ```U``` and ```V``` have directions ```X``` and ```Y```, obviously indicating the direction they are pointing in: they are "directional" fields. In contrast, temperature ```TEMP``` is marked as a scaler, this is a "scalar" field. In general, fields of like direction multiply to yield a new field containing the multiplied nparray data. For instance:

# <codecell>

# This gives a field
U*U

# <markdowncell>

# This gives a vector field of the velocity components squared:

# <codecell>

(U*V)*(U*V)

# <markdowncell>

# Vector fields have a direction method that yields a container of the directions of its component fields.

# <codecell>

my_vector_field.direction()

# <markdowncell>

# The sg module has its own wrapper to the quiver plot method of Matplotlib that takes vector fields:

# <codecell>

# plot surface velocity vector at every other grid cell and up to latitude index of 32.
step = 2
show_j = 32
sg.quiver(U[X,slice(0,100,step),Y,slice(0,show_j,step),Z,0]*V[X,slice(0,100,step),Y,slice(0,show_j,step),Z,0],showland=False)
plt.show()
# <markdowncell>

# The (field) ```Operator``` class and its derived classes (e.g. sg.Der, sg.Integ, sg.Field_Mul etc.) allows the construction of operators that act on fields and vector fields. E.g. div and curl operators, etc. Ideas for useful operators are shown below. Operators, say ```K1``` and ```K2```, can be multiplied before they are applied. E.g. if ```K = K1*K2```, no evaluation has occurred yet, and evaluation yields: ```K*F = K1(K2(F))```. So ```K2``` acts on ```F``` first, then ```K1``` acts on the resulting field. the curl of a 2D vectorfield U*V would then be ```curl*(U*V)```. If ```U*V``` is defined on 3D space, the ```curl``` operator acts on each z-level. For instance, the X-derivative object can be instantiated as follows:

# <codecell>

# Instatiate an X-derivative operator object
ddX = sg.Der(X)
ddX

# <codecell>

ddX*TEMP

# <markdowncell>

# This is equivalent to calling the ```ddX``` field operator object on ```TEMP```:

# <codecell>

ddX(TEMP)

# <markdowncell>

# Applying field operators in succession:

# <codecell>

ddX*ddX

# <markdowncell>

# Where:

# <codecell>

# Create second order X-derivative
d2dX2 = ddX*ddX
# Apply it to the temperature field
d2dX2*TEMP

# <markdowncell>

# Is the same as:

# <codecell>

# Succesive calls of the ddX object
ddX(ddX(TEMP))

# <markdowncell>

# We could have obtained the Z derivative of the vertical temperature profile above with an operator class:

# <codecell>

# Create Z-derivative object
ddZ = sg.Der(Z)
# Apply the derivative operator and plot the result
pl, = sg.plot(ddZ*vpTEMP)
lb = plt.xlabel('degs C per meter')
plt.show()
# <markdowncell>

# Another useful field operator is ```Pick```. It picks a component of a vector field:

# <codecell>

sg.Pick(X)*my_vector_field

# <markdowncell>

# Here are definitions of some useful field operators:

# <codecell>

ddX = sg.Der(X) # derivative in X
ddY = sg.Der(Y) # etc
ddZ = sg.Der(Z)

pX = sg.Pick(X) # pick the component with direction X from vfield object (projection)
pY = sg.Pick(Y) # etc.
pZ = sg.Pick(Z)

mX = sg.Mean(X) # take zonal mean of any field argument or right-multiplicant.
mY = sg.Mean(Y) # etc
mZ = sg.Mean(Z)

sX = sg.Integ(X) # take integral
sY = sg.Integ(Y)
sZ = sg.Integ(Z)

intX = sg.Prim(X) # take primitive function of any field argument
intY = sg.Prim(Y)
intZ = sg.Prim(Z)

set2X = sg.Set_Direction(X) # change the direction attribute of field argument to X
set2Y = sg.Set_Direction(Y)
set2Z = sg.Set_Direction(Z)
set2scalar = sg.Set_Direction(sg.ID())

Avg = mZ*mX*mY  # a 3D average

curl = ddX*set2scalar*pY - ddY*set2scalar*pX      #  curl operator expressed in basic operators

div = ddX*set2scalar*pX + ddY*set2scalar*pY       # divergence 

div3 = ddX*set2scalar*pX + ddY*set2scalar*pY + ddZ*set2scalar*pZ

grad = set2X*ddX + set2Y*ddY    # gradient. To be used on single field objects.

# <codecell>

# Take surface slices
sU = U[Z,0]
sV = V[Z,0]
surface_vectors = sU*sV
(div*surface_vectors).direction

# <codecell>

pl = sg.contourf((div*surface_vectors)[Y,10:-10])

# <codecell>

plt.show()

P.write()

print 'FINISHED.'

