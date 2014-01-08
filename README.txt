README. 

The SG module contains classes to make organizing and analysing Netcdf data more convenient. 

Works with: UVic (all versions), CSIRO Mk3L, FAMOUS, CCM3.6 (Atmosphere), Levitus Data

=============================================================

The SG module consists of one small file: sg.py
This can be obtained from http://web.science.unsw.edu.au/~wsijp/code/sg.py
The manual is located at: http://web.science.unsw.edu.au/~wsijp/code/Manual.txt

The following python modules are required:

numpy
matplotlib
Scientific.IO (for Netcdf)
scipy

A useful page for Netcdf:
http://gfesuite.noaa.gov/developer/netCDFPythonInterface.html


DATA ORGANISATION ON DISK

Before using the SG module, Netcdf data must be organised inside directories you need to create in the following way. Create a directory named PROJECTS inside your home directory ($HOME) or on a storage disk (see below). This will contain all the data for all the projects. Then create a subdirectory for each of your projects. Inside each project, create a small text file named projname containing only a name, say "glacial", to identify your project. In this case, type: echo glacial > projname within the project dir. All subdirectories inside the project directory will now be regarded as potential experiment directories and will be searched for netcdf files when a SG project object P is created (see below). If a netcdf file is found in any of these subdirectories, an exper (experiment) object of the same name as the directory is created and added to the dictionary of experiments belonging to the project object P. Experiment directories can have multiple netcdf files, and can be located inside subdirectories of project directories. 




