### Foreword ###
Welcome to McSAS: a tool for analysis of SAS patterns. 
This tool can extract form-free size distributions from small-angle scattering data using the Monte-Carlo method described in:

Brian R. Pauw, Jan Skov Pedersen, Samuel Tardif, Masaki Takata, and Bo B. Iversen. *“Improvements and Considerations for Size Distribution Retrieval from Small-angle Scattering Data by Monte Carlo Methods.”* Journal of Applied Crystallography 46, no. 2 (February 14, 2013    ). [doi:10.1107/S0021889813001295](http://dx.doi.org/10.1107/S0021889813001295).

The GUI and latest improvements are described in:
Ingo Breßler, Brian R. Pauw, Andreas Thünemann, *"McSAS: A package for extracting quantitative form-free distributions"*. Submitted to J. Appl. Cryst., manuscript available on arXiv:1412.1900(http://arxiv.org/abs/1412.1900).

### Features ###

Several form factors have been included in the package, including:

- Spheres

- Cylinders (spherically isotropic)

- Ellipsoids (spherically isotropic)

- Core-shell spheres and ellipsoids

- Gaussian chain

- Kholodenko worm

- Densely packed spheres (LMA-PY structure factor). 

### Current status ###
The package is currently in beta state, no major issues remain with this version and it should run on a Python 2.7 installation (also works on Enthough Canopy Python).
We are currently working on making standalone packages available for Windows, Linux and Mac OS X (a Windows binary is already available in the downloads section). 
The most up-to-date branch is the "restructuring"-branch, and runs on Windows, Linux and Mac OS X. 
A quick start guide and example data is included in the "doc"-directory that comes with the distribution. 

### Requirements ###
To run McSAS from the source code repository (i.e. using a Python interpreter), the following items are required:

- [Python (latest 2.7 branch)](https://www.python.org/downloads/), with the following packages:

- [Numpy](http://www.scipy.org/scipylib/download.html) 

- [Scipy](http://www.scipy.org/scipylib/download.html) 

- [matplotlib](http://matplotlib.org/downloads.html) 

- [PySide](http://qt-project.org/wiki/PySide) 

### Installation on systems with a working Python distribution ###
For those unfamiliar with the Git versioning system, it is recommended to start by installing Altassian SourceTree (and perhaps reading [Bitbucket 101](https://confluence.atlassian.com/display/BITBUCKET/Bitbucket+101) ). 
This is a GUI around the Git versioning system that simplifies the usage and allows you to get started quickly. 

Using the "clone" button on the top left side of this page, you can download a copy of the latest version. 
Make sure when downloading to select the "restructuring"-branch. 
Following this, McSAS can be started on Unix(-like) systems by opening a terminal window, changing directory to the location of McSAS, and typing 

```
#!shell

$ ./main.py
```

On Windows systems, double-clicking the "main.py" file should open python and start McSAS.

### Standalone packages ###
Standalone packages are available in the [downloads](https://bitbucket.org/pkwasniew/mcsas/downloads) section of this page. 
These are available for Mac OS X (tested on 10.6, 10.8 and 10.10), Windows and Linux. 
These do not require any additional software to be installed on the host computer. 

### Screenshots: ###
![McSAS20150111.png](https://bitbucket.org/repo/jkGXGq/images/2699194750-McSAS20150111.png)
![McSAS20150111Result.png](https://bitbucket.org/repo/jkGXGq/images/4000224154-McSAS20150111Result.png)