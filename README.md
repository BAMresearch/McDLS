## McDLS

Software for the retrieval of model parameter distributions from
dynamic light scattering data.

Output files of a Monte Carlo run are stored in a directory named after
the input file followed with a timestamp to avoid overwriting existing
results.

### Literature

*Coming Soon*

### Current status ###

The package should run on a Python 3 installation, for example an Anaconda environment.
Standalone packages are available for Windows. 

### Requirements ###

To run McDLS from the source code repository (i.e. using a Python interpreter), the following items are required:

- [Python 3](https://www.python.org/downloads/), with the following packages:

- [Numpy](http://www.scipy.org/scipylib/download.html) 

- [Scipy](http://www.scipy.org/scipylib/download.html) 

- [matplotlib](http://matplotlib.org/downloads.html) 

- [PySide2](https://pypi.org/project/PySide2/) 

### Notes on DLS data handling

- The program expects `*.ASC` data files, containing `ALV-7004 CGS-8F Data`
in the first line of the file. It shows a warning otherwise.

- It is recommended to load multiple files from measurements of the same
 sample at once. They are averaged automatically in order to infer
 uncertainties for all values. In a second step, the averaged data set is
 split into one data set for every angle. For example, loading six files of
 measurements at eight angles will result in eight averaged data sets being
 listed in the program.

- Files from more than one sample can be loaded, all files sharing
 the same sample name are averaged along matching angles.

- The measurement numbers ("indices") are taken from the file name of
 each data file. They consist of a group index and a measurement index which
 are shown in the column "Measurements".

- The program processes and fits only those data sets which are selected.
 If none are selected, it processes all data sets. The selection of all data
 sets can be toggled by <Ctrl-A> key presses or via context menu at the
 data set list (right mouse button).

- The menu "Data Settings" shows configuration options for the first of
 the selected data sets. Changes are transferred to all data sets sharing its
 sample name.

- The option "Calc. series statistics" in the optimisation menu combines the
 distribution moments of all data sets in one output file
 (`[...]_seriesStats.dat`) for each active parameter. Finally, it shows a
 simple plot of the mean and its standard deviation across
 the scattering angles.
 
### What to do in case of unsuccessful fits

On success, the resulting size distribution
and data fit are stored to files with uncertainties.
If convergence is not reached no output is generated, and only the
log is stored in a file. The following can be tried:

 1) Adjust the convergence criterion to a larger value. 
    The convergence criterion can be adjusted by the user, to support
    data whose uncertainties are too large or too small.
 2) Verify that the parameter range of the model is appropriate. Very
    wide ranges may prevent the correct solution from appearing within
    the limited number of iterations attempted by the program.
 
