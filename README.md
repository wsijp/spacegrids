### Overview
[![Build Status](https://travis-ci.org/willo12/spacegrids.svg)](https://travis-ci.org/willo12/spacegrids)

Spacegrid is a single python module containing classes to work with and organise data defined on grids (e.g. climate model and observational). It is a bit like "Numpy with grids". The goal is to provide convenient coordinate, grid and data (field) classes that link together, and have short and intuitive data manipulation expressions. In addition, a small framework for project data management is provided. The module provides informative error messages and warnings, and is documented here and extensively inside the python code. 

Installs from pypi with ```pip install spacegrids``` (see https://pypi.python.org/pypi/spacegrids/)

For now, to satisfy dependencies, install on Ubuntu/ Debian:

```apt-get install python-{tk,numpy,matplotlib,scipy} ```

The module consists only of the single file sg.py that can be downloaded from this Github page. The main classes are:

* **coord**	  - grid points along a coordinate axis (discretisation of 1D space).
*  **ax**	  - axis: represents general coord direction: a meta coord. E.g. X or Y.
*  **gr**	  - grid: derived from tuple class, containing coord objects, representing lattices.	
*  **ax_gr**	  - tuple of axes: as gr, but containing ax elements: a meta grid. E.g. (X, Y).
*  **field**	  - includes numpy array of data, has grid in gr attribute.
*  **vfield**   - vector field. tuple-derived container of field elements.
*  **exper**	  - experiment. Represents an experiment directory.
*  **project**  - project. Represents a project directory. 


A dataset is stored in an object, called a "field", that references the grid on which the data is defined, also an object. Data vs grid consistency is ensured at the time of the creation of the field object, and this integrity is maintained throughout  with each grid-related operation on the data, avoiding errors relating to the interpretation of array indices. The aim is to achieve grid operations in a logical intuitive way and with minimal notational effort, with succinct expressions such as F/X denoting zonal averaging. Functionality falls in two main categories: 

1. Object manipulation. This is the main category, and can involve the construction of grid objects from coord objects and or other grid objects, or grid related operations on a field F such as integration and differentiation or interpolation of a field F via the expression F(new_grid) etc.

2. File IO related. Collections of datesets (read from Netcdf) is organized in units named "experiments", and is interpreted and converted into SG objects (fields, grids etc) based on information read from the Netcdf file (not the user). These objects are then arranged in memory within a project structure that reflects a subset of the data directory structure on disk and is easy to reference.

SG also contains wrappers for some of the standard Matplotlib plotting functions that take fields as argument. This can be advantageous as fields contain more information than a regular numpy array, for instance axis labels (that are displayed from within the SG wrappers) and whether the field is defined on a lat-lon grid or whether z should point downward in the plot. Finally, operators such as grad and curl acting on vector fields can be defined in a natural way.

### Getting started and data organization on disk

Getting started is very simple: download the sg.py file from this Github or http://web.science.unsw.edu.au/~wsijp/ and allow Python to find it, and organize some Netcdf data inside directories along the lines described below. 

The module has been developed with Python v. 2.7.3 on Linux with the intention of platform independence (e.g. trying to use os.path.join etc.). The following python modules are required for the SG module: numpy, matplotlib, Scientific.IO (for Netcdf), scipy

A useful page for Netcdf is:
http://gfesuite.noaa.gov/developer/netCDFPythonInterface.html

To manage your data using the SG module, Netcdf data must be organized inside directories you need to create in the following way (the grid related functionality can be used without doing this). Create a directory named PROJECTS inside your home directory ($HOME, the expected default) or in any location of your choosing (see below). PROJECTS will contain all the data for all the projects. Then create a subdirectory for each of your projects. Inside each project, say subdirectory test_project, create a small text file named projname containing only a name, say "glacial", to identify your project. In this case, type: echo glacial > projname within the test_project dir. All subdirectories inside test_project will now be regarded as potential experiment directories and will be searched for netcdf files when a SG project object P is created (see below).  If a netcdf file is found in any of these subdirectories, an exper (experiment) object of the same name as the directory is created and added to the dictionary of experiments belonging to the project object P. Alternatively, if a Netcdf file is found directly inside test_project, it will be interpreted as the output of a single experiment, and an experiment object will be created for it. Experiment directories can have multiple netcdf files, and can be located inside subdirectories of project directories. 

Inside PROJECTS, create a project named test with an experiment named flx_BL via:

```bash
mkdir -p test_project/flx_BL 
echo test > test_project/projname
```

Place Netcdf data files inside the directory flx_BL. If there are time slices, it is best to have them all in one file for each variable. Multiple variables may be spread over multiple files. Different experiments inside the same project may contain output from different models. Different variables within the same experiment may be defined on different grids, first and foremost the tracer and velocity grid, but also different resolutions.

If you have a UVic version 2.8 output file available, all the examples below will work without modification. Otherwise, specific coord and Netcdf variable names will need to be adjusted. 

### Starting a project

Start by obtaining an inventory D (dictionary of project names vs paths) of all the projects with:

```python
D = sg.info()			
```

By default, this expects a directory `PROJECTS` in `$HOME`, but this can be changed via the rootdir arg to `sg.info()`, allowing you to locate `PROJECTS` inside a directory different from `$HOME`.

Now that we know what projects are available, start a project `P` corresponding to a project named 'test'. This name comes from the `projname` file placed inside the project directory (see above):

```python
P = sg.project(D['test']) 	# create project using path. Collects coords etc.
```

We could have used the full path to the project directory, but this is more convenient. Project creation also takes an optional argument "expnames" that can be used to only load from a subset of experiments within the project directory on disk. The "expnames" argument allows globbing via wildcards and other symbols, allowing functionality that is familiar from working with the filesystem. For instance, if the project directory "test" contains additional experiments "SD_BL", "glacial01" and "glacial02", where we are only interested in "flx_BL" and "SD_BL", `P = sg.project(D['test'], expnames = '*_BL')` will only load those experiments. 

The project `P` contains an experiment named 'flx_BL'. `E = P['flx_BL']` provides access to this exper object. To load data from disk into the project attributes (located in memory), use e.g.:

```python
P.load(['sat','temperature','v'])	# v is meridional ocean velocity
```

The data for flx_BL is now accessible as field objects via e.g. `F = P['flx_BL']['sat'] or TEMP = P['flx_BL']['temperature']`.

Different models usually have different naming conventions for variables. If the same project contains experiments done with different models, only those experiments where the variable names match will have their variables loaded. Say model 1 (used in, say, exper E1) refers to certain data by the name of 'temperature', while model 2 (used in exper E2) calls it 'temp', then `P.load(['temperature','temp'])` will load the variable for both experiments (and trigger some benign warnings).  Subsequent care should be taken here in referring to the variable name inside the project P according to the model naming conventions, e.g. use `E2['temp']` and `E1['temperature']`.

Most classes have a copy method that allows convenient duplication of objects and templating (via passed arguments). E.g. `G = F.copy(value = new_value)`.

### Some examples

The following examples are based on two simulations: DPO and DPC, representing Drake Passage open and closed. For these examples to work, UVic v2.9 netcdf diagnostic output data needs to be organized inside subdirectories named DPO and DPC, along with a text file named projname containing the single line "my_project", inside a subdirectory of ~/PROJECTS.

Compressed data files in the correct file hierarchy and the scripts can be downloaded from http://http://web.science.unsw.edu.au/~wsijp/code/examples.tar.gz


Example 1. A figure of a vertical temperature profile based on model output.

```python

# =============================
# This python script creates a figure of vertical temperature profiles of the test model output. The output is not intended to be an accurate approximation of reality, but is for test purposes.
# =============================

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

# This is the SG module:
import spacegrids as sg

# run sg.info() to obtain a dictionary of project names vs project paths. The names are read from text file called "projname" in each (sub)directory of ~/PROJECTS/. The presence of this small text file is used to indicate a project directory. 

D = sg.info()

# start a project object P using the path to that project as argument. This path is easily referenced via the project name and the dictionary D obtained above.
# SG will look through all the subdirectories of that path that contain netcdf files to create experiment objects. (If you put a directory masks in your project path, it will load the masks inside that directory as well.)
# note that D['something'] in the following is a path.

P = sg.project(D['my_project'],expnames = 'DP*')

# Give the experiment objects convenient names E, E2:

E = P['DPO'];E2 = P['DPC']

# bring axes into namespace. e.g. X,Y,Z. The experiment objects will contain a list of axes in their axes attribute.
for c in E.axes:
  exec c.name + ' = c'

varname = 'O_temp'

# We can load fields from the netcdf file. This is where we do the actual IO and load these fields into the memory:

P.load([varname])

#Take vertical profiles of global horizontal averages (convert to Celcius):

mT = E[varname]/(X*Y) 
mT2 = E2[varname]/(X*Y) 

dmT = mT2 - mT

mST = E[varname][Z,0]/X
mST2 = E2[varname][Z,0]/X

# Informative report on deep ocean temperature difference on stdout:
kk=16
print 'T difference 1 at depth '+str((Z*mT.gr)[kk]) +' is ' +str(dmT[kk])

# Parameter preparation for figure layout:

lbl = ord('a')

pan = 0

height = 2
width = 2

rows = 2
cols = 2

# start creating subplots and plot:

ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

plt.yticks(np.arange(-7000,0,1000),np.arange(-7,0,1))

plt.title('('+ chr(lbl) +') Ocean temperature profile ')

plt.grid()

p1, = sg.plot(mT, color = 'k')
p2, = sg.plot(mT2, color = 'r')

plt.xlabel('$^\circ C$')
plt.ylabel('Depth (m)')

plt.legend([p1,p2],['DP open','DP closed'],loc=4)

lbl += 1
pan += 1

ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

plt.title('('+ chr(lbl) +') Zonal avg. SST  ')

plt.grid()

p1, = sg.plot(mST, color = 'k')
p2, = sg.plot(mST2, color = 'r')

plt.xticks([-90,-60,-30,0,30,60,90])

plt.ylabel('$^\circ C$')
plt.xlabel('Latitude')


# could use this to prepare for a next panel in future:
lbl += 1
pan += 1



plt.show()

```


Example 2. Plot the Meridional Overturning Circulation (MOC).

```python
# =============================
# This python script creates a figure of the meridional overturning circulation of the test model output. The output is not intended to be an accurate approximation of reality, but is for test purposes.
# =============================

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

# This is the SG module:
import sg


# run sg.info() to obtain a dictionary of project names vs project paths. The names are read from text file called "projname" in each (sub)directory of ~/PROJECTS/. The presence of this small text file is used to indicate a project directory. 

D = sg.info()

# start a project object P using the path to that project as argument. This path is easily referenced via the project name and the dictionary D obtained above.
# SG will look through all the subdirectories of that path that contain netcdf files to create experiment objects. (If you put a directory masks in your project path, it will load the masks inside that directory as well.)
# note that D['something'] in the following is a path.

P = sg.project(D['my_project'],expnames = 'DP*')


# Give the experiment objects convenient names E, E2:
E = P['DPO'];E2 = P['DPC']

# bring axes into namespace. e.g. X,Y,Z. The experiment objects will contain a list of axes in their axes attribute.
for c in E.axes:
  exec c.name + ' = c'

varname = 'O_velY'

# load velocity by netcdf name.
P.load(varname)

# obtain velocity as sg field object V from project.
V = E[varname]
V2 = E2[varname]


# take the meridional stream function by taking zonal grid-integral X*V and the vertical primitive of the result along Z.
PSI = Z|(X*V)*1e-6
PSI2 = Z|(X*V2)*1e-6

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
cs= sg.contour(F,linestyles='solid',levels = lvls,colors='k');

plt.clabel(cs,fmt='%1.1f');

plt.xlim(xlims)
plt.ylim(ylims)
plt.tick_params(axis='x',labelbottom='off')
plt.xlabel('')
plt.yticks([-0.,-1000,-2000,-3000,-4000])
plt.ylabel('Depth (m)')
plt.title('('+ chr(lbl) +') MOC Drake Passage open ')

lbl += 1
pan += 1

F = PSI2
ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

cs= sg.contour(F,linestyles='solid',levels = lvls,colors='k');

plt.clabel(cs,fmt='%1.1f');

plt.xlim(xlims)
plt.ylim(ylims)
plt.xticks([-60,-30,0,30,60])
plt.yticks([-0.,-1000,-2000,-3000,-4000])
# Overiding the xlabel assigned by sg:
plt.xlabel('Latitude')
plt.title('('+ chr(lbl) +') MOC Drake Passage closed ')

plt.show()
```


Example 3. Plotting the meridional heat transport based on ocean surface heat flux output.


```python
# =============================
# This python script creates a figure of the oceanic poleward heat transport (PHT) of the test model output, for two cases (Drake Passage open and closed). The output is not intended to be an accurate approximation of reality, but is for test purposes.
# =============================

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

# This is the SG module:
import sg


# run sg.info() to obtain a dictionary of project names vs project paths. The names are read from text file called "projname" in each (sub)directory of ~/PROJECTS/. The presence of this small text file is used to indicate a project directory. 

D = sg.info()

# start a project object P using the path to that project as argument. This path is easily referenced via the project name and the dictionary D obtained above.
# SG will look through all the subdirectories of that path that contain netcdf files to create experiment objects. (If you put a directory masks in your project path, it will load the masks inside that directory as well.)
# note that D['something'] in the following is a path.

P = sg.project(D['my_project'],expnames = 'DP*')


# bring axes into namespace. e.g. X,Y,Z. The experiment objects will contain a list of axes in their axes attribute.
for c in P['DPO'].axes:
  exec c.name + ' = c'


# load heat flux by netcdf name.
P.load('F_heat')

# obtain oceanic heat flux as sg field object HF from project.
HF = P['DPO']['F_heat']
HF2 = P['DPC']['F_heat']

# Determine the total oceanic poleward heat transport via the y-primitive and X-intergral, and scale to Petawatt.

PHT = Y|(HF*X)*1e-15
PHT2 = Y|(HF2*X)*1e-15


# --- start figure ---

lbl = ord('a')

pan = 0

height = 2
width = 1

rows = 1
cols = 1

ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

# use the sg wrapper for the plot function.
pl1 = sg.plot(PHT);
pl2 = sg.plot(PHT2);

plt.grid()
plt.xlim([-80.,80.])
plt.xticks([-60,-30,0,30,60])
plt.ylabel('PW')
plt.title('('+ chr(lbl) +') Ocean PHT (PW). ')

lbl += 1
pan += 1


plt.show()
```


### Working with data: fields and their grids

The field object `F` contains the actual data values (e.g. temperature or humidity). This data is accessed directly via `F[:]`, yielding a numpy array. F also has a grid attribute `gr`, e.g. `F.gr = (yt,xt,)`, that points to a spatial lattice `gr` object (e.g. a 100x110 lat lon array of positions where humidity is defined) on which the `F` data is defined. `gr` objects are derived from the tuples class. The example elements `yt` and `xt` here are `coord` objects (other names could be 'longitude', 'latitude' etc) obtained from the Netcdf file. In this case, `xt` contains the X-direction tracer grid positions e.g. `0.,3.75,...,356.25` (similar for `yt` in the Y direction) and the example `coord` name `xt` etc here derives from the UVic v2.8 naming convention. These names are specific to the Netcdf files used and do not need to be supplied by the user. The module code interprets the coordinate data in the Netcdf file, and provides a `name` attribute to the coord object corresponding to the name of the source data as recorded in the Netcdf file. The coord objects yield the grid coordinates along a particular axis `(X,Y,Z,T)` via `xt[:]` etc. 

Grids are tuples of coord objects, e.g. `(yt,xt)`, with extra methods and satisfying certain multiplication rules. They represent latices (arranged discrete collections) of N dimensional spatial coordinates and so represent the familiar 100x200 grid or a 70x100x20 grid etc on which model data might b defined.  For simplicity of use, grids can be constructed by multiplying coord objects. For instance, `yt*xt` yields the 2D lat-lon grid `(yt,xt,)`. 1D grids contain a single `coord` element, e.g. `(xt,)`. A coord element, say `xt`, does not constitute a grid unless contained within a `gr` object (i.e. it is a tuple member). To obtain this 1D `gr` object, use `xt**2`. This yields `(xt,)`.

`coord` objects have an `axis` attribute, which points to an `ax` object that indicates the direction along which the `coord` object is defined, e.g. `X` or `Y` for the X-axis or the Y-axis. Other than specifying a coordinate direction, these `ax` objects can be thought of a meta-`coord` objects: grouping all coord objects pointing in the same direction together. For `xt` and `xu`, where `xu` is the velocity grid along the `X` direction, the `ax` attribute is `X`. `yt.axis` is `Y`. In contrast to `coord` objects, `ax` objects have no values. They retain the multiplication rules of coord objects. In principle there is no limit to the amount of orthogonal axes that can be defined, and N-dimensional spaces can be described. `coord` and `ax` objects from experiment `flx_BL` can be pulled into the namespace by their names via:

```python
for c in P['flx_BL'].cstack:
    exec c.name +' = c'

for l in P['flx_BL'].axes:
    exec l.name + ' = l '
```

where all the `coord` objects for `flx_BL` are in `P['flx_BL'].cstack`, and `ax` objects in `P['flx_BL'].axes`.

Typically, `ax` objects have name attributes `X`, `Y`, `Z` and `T`, and the above `exec` statements bring these objects into the namespace by these names (`T` signifies time here, and is not to be confused with temperature). You can then type something like `X*Y`, which yields the two dimensional grid `(X,Y,)`.



### Operations with fields, grids and axes

Integrals of `F` can be obtained via `X*F`, `Y*F` etc, resulting in `field` objects with grid `(yt,)` resp. `(xt,)`. Primitive functions are obtained via `X|F` etc (yielding `gr` `yt*xt`), mean via `F/X` etc. (`gr` `(yt,)` ) `X^F` yields the derivative in the X-direction.

Examples: 

 * For ocean temperature field `TEMP`, the zonal mean temperature is `TEMP/X`. 
   Resulting field is defined on the 2 dimensional grid `zt*yt`. 

 * A global vertical profile of temperature can be obtained via `TEMP/(X*Y)`. Resulting field is defined on the 1D grid `zt**2`. 

 * For `V = P['flx_BL']['v']`, defined on the 3D `zt*yu*xu` grid, the meridional overturning stream function is `PSI = Z|(X*V)`, defined on the `zt*yu` grid. 

Regrid (interpolate) `F` to a new grid `gr2` by calling `F` on `gr2` via `F(gr2)`. If `F` and `G` correspond the the same physical quantity but loaded from different models (perhaps under different names), `F - G(F.gr)` yields a difference plot for these differently gridded data sets. (This is not recommended for 3D data where all 3 axes are different as this can take a long time to interpolate.) To specify the interpolation method, the argument 'methods' can be specified, e.g. F(gr2,method = 'nearest'). The default is 'linear'.

Examples: 

 * `V(zt*yt*xt)` yields `V` regridded to the tracer grid.

 * Multiplicants of different dimensions. For `F` defined on `yt*xt` and `gr2 = zt*yt*xt`, `F(gr2)` is defined on `gr2`, with `F[:]` constant along the (new) `zt`-axis: the original has been copied `len(zt)` times. 

 * If `gr2 = xt*yt`, `F(gr2)[:]` will be the transposed of `F[:]`. 

 * If `gr2` is a different grid `yt2*xt2`, then `F(yt2*xt2)` yields `F` interpolated on that new grid. 

Field multiplication `F*G` is element-wise on the data when `F.gr == G.gr`. If the grids match up to a permutation, the data arrays associated with the field objects are automatically tranposed to match the grid of the left multiplicant. So the order of the coordinates in the multiplicants need not match. 

Example:  

 *  If `F.gr` is `zt*yt*xt` and `G.gr` is `yt*xt*zt`, then `G` is re-arranged before multiplication and `(F*G).gr` is `zt*yt*xt`. If `G.gr` is `yt*xt`, then every z-level of `F` is multiplied with the 2D data array `G[:]` and the product grid is `zt*yt*xt`.

If a field `G` has a different grid `gr2` but spanning the same subspace as `F.gr`, then `G*F = G*F(gr2)`, while `F*G` yields `F*G(F.gr)`: multiplication triggers interpolation here.

slicing can be done in the numpy way, or using coord objects. `F[xt,1:]` excludes the first slice in the `X`-direction, yielding a new field with a new grid `xt_sliced*yt`, where `xt_sliced` has been reduced accordingly. The specific `coord` object can be replaced by its `ax` object, so that the same slicing can be achieved by `F[X,1:]`. This substitution of specific `coord` objects with their general `ax` object can often be done in operations. 

`F.shape` yields `F[:].shape` and must be identical to `F.gr.shape()`. Field creation via `G = sg.field(value = A, grid = my_gr)`  will not proceed unless `A.shape = my_gr.shape()` for numpy array `A`. Field `gr` attributes are therefore always consistent with the data array attribute, where the order of the grid `coord` elements corresponds to the order of the array indices. 



### More about the coord and gr class

The `coord` class has subsclasses `ycoord` and `xcoord`. `coord` objects have an equivalence relation. It is set via `xt|xu`, where `xu` is another `coord`. `xt^xu` tests for equivalance (True here). Usually, `^` is used to indicate that `xt` and `xu` are defined along the same direction (`X`). Here, `xu` could be the velocity grid and `xt` the tracer grid that is offset from it. `xu*xt` yields `xu`, whereas `xt*xu` yields `xt`: the left multiplicant is always selected. `xt*yt` yields `(xt,yt,)`, as `xt^yt` is `False`. Also, `X*xt = xt`, `X*xu = xu`, `X*(xt*yt) = xt` etc. Much of the `X*F` type functionality is based on these rules: the user does not need to know the specific coord objects, just the `ax` objects. 
`gr` objects can be multiplied: e.g. `(zt*yu*xu)*(yt*xt) = zt*yu*xu`, where equivalent `coord` objects interact by retaining the left multiplicant. 

In the most general case, field multiplication works as follows:


```python
F * G = F(F.gr * G.gr) * G(F.gr * G.gr)
```

As a practical example, if `M` were a field object containing a horizontal mask containing values of 1 and 0 or 1 and nan defined on `yt*xt` then `M*TEMP` masks out the required values in the same locations for all z-levels (each z-slice is masked the same way).

gr objects take other `gr` objects as argument. This yields a (lambda) function that takes a numpy array argument and re-arranges and or interpolates it from the self grid to the other grid. E.g. `xt*yt(yt*xt) `yields a transpose operation on an array. Another example is the function `xt*yt(xu*yu)`. This is an interpolation acting on numpy arrays. Furthermore, the value of `(F*G).value = F.gr(F.gr*G.gr)(F.value)*G.gr(F.gr*G.gr)(G.value)`. Recall that `F.value` and `G.value` are numpy arrays.



### Vector fields

A vector field object is a collection of directional fields (derived from tuples and with additional attributes and methods). Fields have a direction attribute that contains (points to) an `ax` objects.  In directional fields, or vector field components, this is an `ax` object representing the `X`,`Y` or `Z` direction. E.g. `V.direction = Y` and `U.direction = X`. This works through in their multiplication rules: `U*V = (U,V)`, yielding a 2D `vfield` object, while `U*U` yields a `field`. The direction attribute of non-directional fields is named "scalar", the `ax` object defined as `sg.ID()`, representing no direction. Enforcement of grid consistency between `vfield` components has not been implemented yet: it is possible to define a `vfield` using (directional) fields each defined on different grids (e.g. an unhelpful situation where `U` is defined on a 100x200 grid and `V` on 200x150). Their creation invokes only a warning, and care should be taken.



### Operators

The Operator class allows the definition of operators that generally act on and yield (v)fields. E.g. div and curl operators, etc. Ideas for useful operators are shown below. Operators, say `K1` and `K2`, can be multiplied before they are applied. E.g. if `K = K1*K2`, no evaluation has occurred yet, and evaluation yields: `K*F = K1(K2(F))`. So `K2` acts on `F` first, then `K1` acts on the resulting field. the curl of a 2D vectorfield `U*V` would then be `curl*(U*V)`. If `U*V` is defined on 3D space, the `curl` operator acts on each z-level.



### Plotting

The field class has a draw method to allow a quick visual inspection of the field. The SG module also contains wrappers for the standard Matplotlib plotting functions that take field arguments and can make informed guesses about the direction of axes in the figure and label them.



### Useful operators

```python
ddX = sg.Der(X)	# derivative in X
ddY = sg.Der(Y) # etc
ddZ = sg.Der(Z)

pX = sg.Pick(X) # pick the component with direction X from vfield object
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

Avg = mZ*mX*mY  # a 3D average

curl = ddX*pY - ddY*pX 		#  curl operator expressed in basic operators

div = ddX*pX + ddY*pY 		# divergence 

div3 = ddX*pX + ddY*pY + ddZ*pZ

grad = set2X*ddX + set2Y*ddY	# gradient. To be used on single field objects.

```
