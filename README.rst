Spacegrids
==========

Spacegrids is an open source library providing a Numpy array with grids, labelled axes and associated grid-related mathematical methods such as regridding and integration. Spacegrids provides an object data model of Netcdf data that ensures consistency between a Numpy data array and its grid under common operations (and so avoiding common pitfalls related to axis interpretation), and much more. It is a write less do more library for everyday use.

The Field, Gr (grid) and Coord objects make everyday use easy:

    >>> import spacegrids as sg		
    >>> D = sg.info(nonick = True)  
    >>> P = sgPproject(D['my_project'] , nonick = True)  
    >>> P.load(['temperature','u'])  
    >>> # obtain the axes under their names T,X,Y,Z in namespace:
    >>> for c in P['some_experiment'].axes:
    >>>   exec c.name + ' = c'	
    >>> TEMP = P['some_experiment']['temperature'] 
    >>> U = P['some_experiment']['u'] # zonal velocity
    >>> TEMP_sliced = TEMP[Y,:50] # slice in Y-direction
    >>> m_TEMP = TEMP_sliced/(X*Y) # take zonal mean
    >>> TEMP_regridded = TEMP.regrid(U.gr)  # U on different grid
 

Features
--------

- A numpy array with grid allowing automatic alignment and dimension broadcasting
- Easy to use and intuitive regridding functionality
- A data object model corresponding closely to Netcdf
- Easier IO via abstraction of IO with multiple Netcdf files
- Makes working with output of many experiments easy via aggregation methods
- The Field class eliminates errors arising from picking the wrong array index
- Quicker plotting due to automatic labels, axes etc.
- Distance-related methods such as spatial differentiation and integration on sphere
- Extensive unit tests and documentation

There is lots of documentation, both in the source code and elsewhere. Other documentation can be found at: 

- `a practical tutorial <http://nbviewer.ipython.org/github/willo12/spacegrids/blob/master/Spacegrids.ipynb>`_ 
- `a more advanced tutorial <http://nbviewer.ipython.org/github/willo12/spacegrids/blob/master/advanced.ipynb>`_ 
- `an overview of all classes, methods and functions <http://web.maths.unsw.edu.au/~wsijp/html/index.html>`_ 


Installation
------------

Install spacegrids simply by running (on command line):

    pip install spacegrids

On Mac, pip can be installed via "sudo easy_install pip". On Ubuntu/ Debian, install dependencies via package manager if pip install fails:

    apt-get install python-{tk,numpy,matplotlib,scipy}


Contribute
----------

- Issue Tracker: github.com/willo12/spacegrids/issues
- Source Code: github.com/willo12/spacegrids

Support
-------

If you are having issues, please let us know.

License
-------

The project is licensed under the BSD license.
