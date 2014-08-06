Spacegrids
==========

Spacegrids will solve your problem of customizing your Python scripts for every new data analysis project and avoid common pitfalls related to grid definition by providing an object data model of Netcdf data that preserves grid definition under common operations and much more. It is a write less do more library for everyday use to enhance productivity.

The Field, Gr (grid) and Coord objects make everyday use easy:

    import spacegrids as sg		
    D = sg.info(nonick = True)  
    P = sgPproject(D['my_project'] , nonick = True)  
    P.load(['temperature','u'])  
    # obtain the axes under their names T, X, Y, Z in namespace:
    for c in P['some_experiment'].axes:
      exec c.name + ' = c'	
    TEMP = P['some_experiment']['temperature'] 
    U = P['some_experiment']['u'] # zonal velocity
    TEMP_sliced = TEMP[Y,:50] # slice in Y-direction
    m_TEMP = TEMP_sliced/(X*Y) # take zonal mean
    TEMP_regridded = TEMP.regrid(U.gr)  # U is on a different grid to TEMP
 

Features
--------

- A numpy array with grid allowing automatic alignment and dimension broadcasting
- This leads to easy to use and intuitive regridding functionality
- A data object model corresponding closely to Netcdf
- Easier IO via abstraction of IO with multiple Netcdf files
- Makes working with output of many experiments easy
- The Field class eliminates errors arising from picking the wrong array index
- Quicker plotting due to automatic labels, axes etc.
- Takes grid geometry and definition into account
- Distance-related methods such as spatial differentiation and integration on sphere


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
