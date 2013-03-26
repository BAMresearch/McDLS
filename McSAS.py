# -*- coding: utf-8 -*-
# McSAS.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

r"""
Overview
========
A class and supplementary functions for Monte-Carlo fitting of SAXS patterns.
It is released under a `Creative Commons CC-BY-SA license
<http://creativecommons.org/licenses/by-sa/3.0/>`_.
Please cite as::

    Brian R. Pauw et al., J. Appl. Cryst. 46, (2013), pp. 365--371
        doi: http://dx.doi.org/10.1107/S0021889813001295

Classes and Functions Defined in This File
------------------------------------------

 - :py:class:`McSAS() <McSAS.McSAS>`:
   A class containing all the functions required to perform a
   Monte Carlo analysis on small-angle scattering data.
 - :py:func:`binning_array`:
   Can be used to do n-by-n pixel binning of 2D detector
   images. The returned uncertainty is the larger of either the binned
   uncertainty or the sample standard deviation in the bin.
 - :py:func:`binning_1D`:
   bins the data and propagates errors, or calculates errors
   if not initially provided
 - :py:func:`binning_weighted_1D`:
   Weighted binning, where the intensities of a
   pixel are divided between the two neighbouring bins depending on the
   distances to the centres. If error provided is empty, the standard
   deviation of the intensities in the bins are computed.
 - :py:func:`pickle_read`:
   Reads in pickled data from a file (by filename)
 - :py:func:`pickle_write`:
   write a block or dictionary to a file (by filename)
    
Made possible with help from (amongst others)
---------------------------------------------

 - | Samuel Tardif
   | Derivations (mostly observability) and checking of mathematics
 - | Jan Skov Pedersen
   | checking of mathematics
 - | Pawel Kwasniewski <kwasniew@esrf.fr>
   | Code cleanup and documentation
 - | Ingo Bressler <ingo.bressler@bam.de>
   | Code cleanup, modification and documentation

A Note on Units
---------------

Internally, all length units are in meters, all angle units in degrees
clockwise from top. *Intensity* is in
:math:`\left[ 1 \over {m \cdot sr} \right]`,
*q* in :math:`\left[ 1 \over m \right]`.
The electron density contrast squared,
*drhosqr* is in :math:`\left[ m^{-4} \right]`.
Other units may be used, but if absolute units are supplied and absolute
volume fractions required, meters are necessitated.

Example Usage
-------------

*For detailed usage, please see the* :doc:`quickstart`

Fitting a single dataset using all automatic and default parameters
(may go wrong on poorly conditioned input, needs sensibly-spaced datapoints
and good uncertainty estimates).
The dataset is considered to consist of three variables *Q*, *I* and *IERR*::

 McSAS(Q = Q, I = I, IERR = IERR, Plot = True)

Optional parameters can be supplied in parameter-value pairs to finetune
optimisation behaviour::

 A = McSAS(Q = Q, I = I, IERR = numpy.maximum(0.01 * I, E),
           Ncontrib = 200, Convcrit = 1,
           Bounds = array([0.5e-9, 35e-9]),
           Maxiter = 1e5, Histscale = 'log',
           drhosqr = 1e30, Nreps = 100, Plot = True)

Module Documentation
====================
"""

import scipy # For many important functions
from scipy import optimize # For the leastsq optimization function
import numpy # For arrays
from numpy import (inf, array, isfinite, reshape, prod, shape, pi, diff, zeros,
                  arange, size, sin, cos, sum, sqrt, linspace, logspace, log10,
                  isnan, ndim)
import os # Miscellaneous operating system interfaces
import time # Timekeeping and timing of objects
import sys # For printing of slightly more advanced messages to stdout
import pickle #for pickle_read and pickle_write

class McSAS(object):
    r"""
    Main class containing all functions required to do Monte Carlo fitting.

    **Required input parameters:**

        - *Q*: 1D or 2D array of q-values
        - *I*: corresponding intensity values of the same shape
        - *IERR*: corresponding intensity uncertainty values of the same shape

    **Optional input parameters:**

        - *PSI*: 2D array
            Detector angle values, only required for 2D pattern fitting.
        - *Bounds*: list
            Two-element vector or list indicating upper and lower size
            bounds of the particle radii used in the fitting procedure. If
            not provided, these will be estimated as:
            :math:`R_{max} = {pi \over q_{min}}` and
            :math:`R_{min} = {pi \over q_{max}}`. Units in meter.
        - *Nsph*: int, default: 200
            Number of spheres used for the MC simulation
        - *Maxiter*: int, default: 1e5
            Maximum number of iterations for the :py:func:`MCFit` function
        - *Rpfactor*: float, default: :math:`1.5 \over 3`
            Parameter used to compensate the :math:`volume^2` scaling of each
            sphere contribution to the simulated I(q).
        - *Nreps*: int, default: 100
            Number of repetitions of the MC fit for determination of final
            histogram uncertainty.
        - *qlims*: list, default: [0, inf]
            Limits on the fitting range in q.
            Units in :math:`m^{-1}`
        - *Histbins*: int, default: 50
            Number of bins used for the histogramming procedure.
        - *Histscale*: string, default: 'log'
            Can be set to 'log' for histogramming on a logarithmic size scale, 
            recommended for q- and/or size-ranges spanning more than a decade.
        - *Histweight*: string, default: 'volume'
            Can be set to 'number' to force plotting of number-weighted
            distributions
        - *drhosqr*: float, default: 1
            Scattering contrast - when known it will be used to calculate the
            absolute volume fraction of each contribution.
            Units in :math:`m^{-4}`
        - *Convcrit*: float, default: 1
            Convergence criterion for the least-squares fit. The fit converges
            once the :math:`normalized \chi^2 < Convcrit`. If convergence is
            reached with `Convcrit == 1`, the model describes
            the data (on average) to within the uncertainty, and thus all
            information has been extracted from the scattering pattern.
        - *StartFromMin*: bool, default: False
            If set to False, the starting configuration is a set of spheres
            with radii uniformly sampled between the given or estimated
            bounds. If set to True, the starting configuration is a set of
            spheres with radii set to the lower given or estimated Bound
            (if not zero). Practically, this makes little difference and this
            feature might be depreciated.
        - *Maxntry*: int, default: 5
            If a single MC optimization fails to reach convergence within
            *Maxiter*, it may just be due to bad luck. The procedure will try
            to redo that MC optimization for a maximum of *Maxntry* tries
            before concluding that it is not bad luck but bad input.
        - *Plot*: Bool, default: False
            If set to True, will generate a plot showing the data and fit, as
            well as the resulting size histogram.
        - *Memsave*: Bool, default: False
            For 2D pattern fitting, or for fitting patterns with a very large
            number of datapoints or contributions, it may make sense to turn
            this option on in order for intensity generating functions not to
            take up much memory. The cost for this is perhaps a 20-ish percent
            reduction in speed.
        - *BOUNDS*: string
            The McSAS function to use for calculating random number generator 
            bounds based on input (f.ex. q and I).
            default: :py:func:`SphBounds`
        - *FF*: string
            The McSAS function to use for calculating the form factors.
            default: :py:func:`FF_sph_1D`
        - *RAND*: string
            the McSAS function to use for generating random numbers
            default: :py:func:`random_uniform_sph`
        - *SMEAR*: string
            the McSAS function to use for smearing of intensity
            default: :py:func:`_passthrough`
        - *VOL*: string
            the McSAS function to use for calculating the base object volume
            default: :py:func:`vol_sph`

    **Returns:**

    A McSAS object with the following results stored in the *result* member
    attribute. These can be extracted using
    McSAS.getresult('Keyword',VariableNumber=0)
    where the *VariableNumber* indicates which shape parameter information is
    requested for
    (some information is only stored in *VariableNumber = 0* (default)).

    **Keyword** may be one of the following:

        *Imean*: 1D array (*VariableNumber = 0*)
            The fitted intensity, given as the mean of all Nreps results.
        *q*: 1D array (*VariableNumber = 0*)
            Corresponding q values
            (may be different than the input q if *qlims* was used).
        *Istd*: array (*VariableNumber = 0*)
            Standard deviation of the fitted I(q), calculated as the standard 
            deviation of all Nreps results.
        *Rrep*: size array (Nsph x Nreps) (*VariableNumber = 0*)
            Collection of Nsph sphere radii fitted to best represent the
            provided I(q) data. Contains the results of each of *Nreps*
            iterations. This can be used for rebinning without having to
            re-optimize.
        *Screps*: size array (2 x Nreps) (*VariableNumber = 0*)
            Scaling and background values for each repetition.
            Used to display background level in data and fit plot.
        *VariableNumber*: int
            Shape parameter index.
            E.g. an ellipsoid has 3: width, height and orientation.
        *Hx*: array
            Histogram bin left edge position (x-axis in histogram).
        *Hmid*: array
            Center positions for the size histogram bins
            (x-axis in histogram, used for errorbars).
        *Hwidth*: array
            Histogram bin width
            (x-axis in histogram, defines bar plot bar widths).
        *Hmean*: array
            Volume-weighted particle size distribution values for
            all Nreps results (y-axis bar height).
        *Hnmean*: array
            Number-weighted analogue of the above *Hmean*.
        *Hy*: size array (Histbins x Nreps)
            Volume-weighted particle size distribution bin values for
            each MC fit repetition (the mean of which is *Hmean*, and the
            sample standard deviation of which is *Hstd*).
        *Hny*: size array (Histbins x Nreps)
            Number-weighted particle size distribution bin values for
            each MC fit repetition.
        *Hstd*: array
            Standard deviations of the corresponding volume-weighted size
            distribution bins, calculated from *Nreps* repetitions of the
            :py:meth:`McSAS.MCfit_sph` function.
        *Hnstd*: array
            Standard deviation for the number-weigthed distribution.
        *Vf*: size array (Ncontrib x Nreps)
            Volume fractions for each of Ncontrib contributions in each of
            *Nreps* iterations.
        *Nf*: size array (Ncontrib x Nreps)
            Number fraction for each contribution.
        *Vft*: size array (Nreps)
            Total scatterer volume fraction for each of the Nreps iterations.
        *Nft*: size array (Nreps)
            Total number fraction.
        *vfmin*: size array (Nsph x Nreps)
            Minimum required volume fraction for each contribution to become
            statistically significant.
        *nfmin*: size array (Nsph x Nreps)
            Number-weighted analogue to *vfmin*.
        *vfminbins*: size array (Hmid)
            Array with the minimum required volume fraction per bin to become
            statistically significant. Used to display minimum required level
            in histogram.
        *nfminbins*: size array (Hmid)
            Number-weighted analogue to *vfminbins*.
        *Screps*: size array (2 x Nreps)
            Scaling and background values for each repetition. Used to display
            background level in data and fit plot.
        *Vf*: size array (Nsph x Nreps)
            Volume fractions for each of *Nsph* spheres in each of *Nreps*
            iterations.
        *Vft*: size array (Nreps)
            Total scatterer volume fraction for each of the *Nreps*
            iterations.
        *vfmin*: size array (Nsph x Nreps)
            Minimum required volube fraction for each contribution to become
            statistically significant.
        *vfminbins*: size array (Hmid)
            Array with the minimum required volume fraction per bin to become
            statistically significant. Used to display minimum required level
            in histogram.
    """

    dataset=None #where Q, PSI, I and IERR is stored, original dataset
    fitdata=None #may be populated with a subset of the aforementioned dataset, limited to q-limits or psi limits and to positive I values alone
    parameters=None #where the fitting and binning settings are stored
    result=None #where all the analysis results are stored, I do not think this needs separation after all into results of analysis and results of interpretation. However, this is a list of dicts, one per variable (as the method, especially for 2D analysis, can deal with more than one random values. analysis results are stored along with the histogrammed results of the first variable with index [0]:
    functions=None #where the used functions are defined, this is where shape changes, smearing, and various forms of integrations are placed.

    def __init__(self,**kwargs):
        """
        The constructor, takes keyword-value input parameters. They can be
        one of the aforementioned parameter keyword-value pairs.
        This does the following:

            1. Initialises the variables to the right type
            2. Parses the input
            3. Stores the supplied data twice, one original and one for fitting 
                (the latter of which may be clipped or otherwise adjusted)
            4. Applies Q- and optional PSI- limits to the data
            5. Reshapes fitdata to 1-by-n dimensions
            6. Sets the function references
            7. Calculates the shape parameter bounds if not supplied
            8. Peforms simple checks on validity of input
            9. Runs the Analyse() function which applies the MC fit multiple times
            10. Runs the Histogram() procedure which processes the MC result 
            11. Optionally recalculates the resulting intensity in the same shape 
                as the original (for 2D datasets)
            12. Optionally displays the results graphically.

        """
        #initialize
        self.dataset=dict() #where Q, PSI, I and IERR is stored, original dataset
        self.fitdata=dict() #may be populated with a subset of the aforementioned dataset, limited to q-limits or psi limits and to positive I values alone
        self.parameters=dict() #where the fitting and binning settings are stored
        self.result=list() #where all the analysis results are stored, I do not think this needs separation after all into results of analysis and results of interpretation. However, this is a list of dicts, one per variable (as the method, especially for 2D analysis, can deal with more than one random values. analysis results are stored along with the histogrammed results of the first variable with index [0]:
        self.result.append(dict())
        self.functions=dict() #where the used functions are defined, this is where shape changes, smearing, and various forms of integrations are placed.

        #populate self with defaults
        self.set_defaults()
        #set supplied kwargs
        self.setpar(**kwargs) #passing on kwargs
        #set data values
        self.setdata(**kwargs)
        #apply q and psi limits and populate self.fitdata
        self.clip_dataset()
        #reshape fitdata to the correct 1-by-n dimensions
        self.reshape_fitdata()
        #apply input settings for fitting, setting the required function definitions
        self.setfunction(**kwargs)
        #calculate parameter bounds and store
        self.getfunction('BOUNDS')()
        #check and fix parameters where necessary
        self.check_parameters() #this is only a very simple function now in need of expansion

        #Analyse
        self.Analyse()

        ##Histogram
        self.Histogram()

        if ndim(kwargs['Q'])>1:
            #2D mode, regenerate intensity
            self.TwoDGenI()

        ##Plot
        if self.getpar('Plot'):
            #Odata=self.getdata(dataset='original')
            #Result=self.getresult()
            self.Plot()

    ############################################################################
    ##################### Pre-optimisation functions ###########################
    ############################################################################
    def _Iopt(self,I,Ic,E,Sc,ver=2,OutputI=False,Background=True):
        """
        Optimizes the scaling and background factor to match Ic closest to I.
        Returns an array with scaling factors. Input initial guess *Sc* has 
        to be a two-element array with the scaling and background.

        ** Required input arguments: **
            - *I* : An array of "measured" intensities
            - *Ic* : An array of intensities which should be scaled to match *I*
            - *E* : An array of uncertainties to match *I*
            - *Sc* : A 2-element array of initial guesses for scaling factor 
                and background

        ** Optional input arguments: **
            - *ver*: can be set to 1 for old version, more robust but slow, 
              default 2 for new version, 10x faster than version 1
            - *OutputI*: returns the scaled intensity as third output argument,
              default: False
            - *Background*: enables a flat background contribution, 
              default: True

        Output: 
            - *Sc*: a two-element array with intensity scaling factor and background
            - *cval*: the reduced chi-squared value
        """
        def csqr(Sc,I,Ic,E):
            """Least-squares error for use with scipy.optimize.leastsq"""
            cs=(I-Sc[0]*Ic-Sc[1])/E
            return cs
        
        def csqr_noBG(Sc,I,Ic,E):
            """Least-squares error for use with scipy.optimize.leastsq, without background """
            cs=(I-Sc[0]*Ic)/E
            return cs

        def csqr_v1(I,Ic,E):
            """Least-squares for data with known error, size of parameter-space not taken into account."""
            cs=sum(((I-Ic)/E)**2)/(size(I))
            return cs

        if ver==2:
            """uses scipy.optimize.leastsqr"""
            if Background:
                Sc,success=scipy.optimize.leastsq(csqr,Sc,args=(I.flatten(),Ic.flatten(),E.flatten()),full_output=0)
                cval=csqr_v1(I,Sc[0]*Ic+Sc[1],E)
            else:
                Sc,success=scipy.optimize.leastsq(csqr_noBG,Sc,args=(I.flatten(),Ic.flatten(),E.flatten()),full_output=0)
                Sc[1]=0.
                cval=csqr_v1(I,Sc[0]*Ic,E)
        else:
            """using scipy.optimize.fmin"""
            # Background can be set to False to just find the scaling factor.
            if Background:
                Sc=scipy.optimize.fmin(lambda Sc: csqr_v1(I,Sc[0]*Ic+Sc[1],E),Sc,full_output=0, disp=0)
                cval=csqr_v1(I,Sc[0]*Ic+Sc[1],E)
            else:
                Sc=scipy.optimize.fmin(lambda Sc: csqr_v1(I,Sc[0]*Ic,E),Sc,full_output=0, disp=0)
                Sc[1]=0.
                cval=csqr_v1(I,Sc[0]*Ic,E)
        if OutputI:
            return Sc,cval,Sc[0]*Ic+Sc[1]
        else:
            return Sc,cval

    def set_defaults(self):
        """
        Populates the default parameter settings
        """
        #fieldnames
        #fnames=list(['Bounds','Ncontrib','Maxiter','Rpfactor','Nreps','qlims','psilims','Histbins','Histscale','drhosqr','Convcrit','StartFromMin','Maxntry'])
        self.parameters={
                'Bounds':[],
                'Ncontrib':200,
                'Maxiter':1e5,
                'Rpfactor':0.5,
                'Nreps':100,
                'qlims':[],
                'psilims':[],
                'Priors':[], #of shape Rrep, to be used as initial guess for Analyse function, Analyse will pass on a Prior to MCFit.
                'Prior':[], #of shape Rset, to be used as initial guess for MCFit function
                'Histbins':50,
                'Histscale':'log',
                'Histweight':'volume', #can be set to "volume" or "number"
                'drhosqr':1,
                'Convcrit':1.,
                'StartFromMin':False,
                'Maxntry':5,
                'MaskNegI':False,
                'Memsave':False,
                'Plot':False}

        self.functions={
                'BOUNDS':self.SphBounds, #this function has to give a vector the size of the number of variables *2 (lower and upper)
                'RAND':self.random_uniform_sph,
                'FF':self.FF_sph_1D, #1D spheres
                'VOL':self.vol_sph,
                'SMEAR':self._passthrough} #none

    def setfunction(self,**kwargs):
        """Defines functions. In particular the following are specified:

        - The parameter bounds estimation function *BOUNDS*. Should be able
          to take input argument Bounds to update, should set the parameter
          bounds in ``self.parameter['Bounds']``

        - The random number generation function *RAND* This must take its
          parameters from self, and have an optional input argument specifying
          the number of sets to return (for MC initialization). It should
          return a set of Nsets-by-nvalues to be used directly in *FF*. This
          may be depreciated soon as is can be generated from within.

        - The Form-factor function *FF*. If called, this should get the
          required information from self and a supplied Nsets-by-nvalues
          shape-specifying parameter array. It should return an Nsets-by-q
          array. Orientational averaging should be part of the form-factor
          function (as it is most efficiently calculated there), so several
          form factor functions can exist for non-spherical objects.

        - The shape volume calculation function *VOL*, which must be able to
          deal with input argument *Rpfactor*, ranging from 0 to 1. Should
          accept an Nsets-by-nvalues array returning an Nsets number of
          (Rpfactor-compensated)-volumes. 

        - The smearing function *SMEAR*. Should take information from self
          and an input Icalc, to output an Ismear of the same length.

        This function will actually use the supplied function name as function
        reference.
        """
        for kw in kwargs:
            if kw in self.functions.keys():
                if callable(kwargs[kw]):
                    self.functions[kw]=kwargs[kw]
                else:
                    #make it into a function handle/pointer
                    self.functions[kw]=getattr(self,kwargs[kw])

    def getfunction(self,fname=None):
        '''
        returns the function handle or all handles (if no function name supplied).
        function name can be one of the following
        '''
        fname=fname.upper()
        if not fname in self.functions.keys():
            print "Unknown function identifier {}".format(fname)
            return None
        if fname==None:
            return self.functions
        else:
            return self.functions[fname]


    def random_uniform_ell(self,Nell=1):
        """Random number generator for generating uniformly-sampled
        size- and orientation parameters for ellipsoids.
        """
        #get parameters from self
        Bounds=self.getpar('Bounds') 
        #generate Nsph random numbers
        Rset=zeros((Nell,3))
        Rset[:,0]=numpy.random.uniform(numpy.min(Bounds[0]),numpy.max(Bounds[1]),Nell)
        Rset[:,1]=numpy.random.uniform(numpy.min(Bounds[2]),numpy.max(Bounds[3]),Nell)
        Rset[:,2]=numpy.random.uniform(numpy.min(Bounds[4]),numpy.max(Bounds[5]),Nell)

        #output Nsph-by-3 array
        return Rset

    def random_logR_ell(self,Nell=1):
        """Random number generator which behaves like its uniform counterpart,
        but with a higher likelihood of sampling smaller sizes.
        May speed up some fitting procedures.
        """
        #get parameters from self
        Bounds=self.getpar('Bounds') 
        #generate Nsph random numbers
        Rset=zeros((Nell,3))
        Rset[:,0]=10**(numpy.random.uniform(log10(numpy.min(Bounds[0])),log10(numpy.max(Bounds[1])),Nell))
        Rset[:,1]=10**(numpy.random.uniform(log10(numpy.min(Bounds[2])),log10(numpy.max(Bounds[3])),Nell))
        Rset[:,2]=numpy.random.uniform(numpy.min(Bounds[4]),numpy.max(Bounds[5]),Nell)

        #output Nsph-by-3 array
        return Rset

    def random_logR_sph(self,Nsph=1):
        """Random number generator with logarithmic probabilistic sampling."""
        #get parameters from self
        Bounds=self.getpar('Bounds') 
        #generate Nsph random numbers

        Rset=10**(numpy.random.uniform(log10(numpy.min(Bounds)),log10(numpy.max(Bounds)),Nsph))
        Rset=reshape(Rset,(prod(shape(Rset)),1))

        #output Nsph-by-1 array
        return Rset

    def random_uniform_sph(self,Nsph=1):
        """Random number generator with uniform distribution for
        the sphere form factor."""
        #get parameters from self
        Bounds=self.getpar('Bounds') 
        #generate Nsph random numbers

        Rset=numpy.random.uniform(numpy.min(Bounds),numpy.max(Bounds),Nsph)
        Rset=reshape(Rset,(prod(shape(Rset)),1))

        #output Nsph-by-1 array
        return Rset

    def vol_ell(self,Rset,Rpfactor=[]):
        """Calculates the volume of an ellipsoid, taking Rpfactor from input
        or preset parameters.
        """
        if Rpfactor==[]:
            Rpfactor=self.getpar('Rpfactor')
            
        return ((4.0/3*pi)*Rset[:,0]**(2*Rpfactor)*Rset[:,1]**(Rpfactor))[:,newaxis]

    def vol_sph(self,Rset,Rpfactor=[]):
        """Calculates the volume of a sphere, taking Rpfactor from input
        or preset parameters.
        """
        if Rpfactor==[]:
            Rpfactor=self.getpar('Rpfactor')
            
        return (4.0/3*pi)*Rset**(3*Rpfactor)

    def FF_sph_1D(self,Rset):
        """Calculate the Rayleigh function for a sphere.
        """
        q=self.getdata('Q')
        if size(Rset,0)>1: # multimensional matrices required, input Rsph has to be Nsph-by-1. q has to be 1-by-N
            qR=(q+0*Rset)*(Rset+0*q)
        else:
            qR=(q)*(Rset)

        Fsph=3*(sin(qR)-qR*cos(qR))/(qR**3)
        return Fsph

    def FF_ell_2D(self,Rset=[],Q=[],PSI=[]):
        """Calculates the form factor for oriented ellipsoids,
        normalized to 1 at Q=0.

        **Note**: All 2D functions should be able to potentially take
        externally supplied Q and PSI vectors.
        """
        #Rset is n-by-3. R1=Rset[:,0],R2=Rset[:,1],R3=Rset[:,2]
        #R1<R2, prolate ellipsoid (cigar-shaped), R1>R2, oblate ellipsoid (disk-shaped), rotation is offset from perfect orientation (psi-rot)
        d_to_r=1./360*2*pi #degrees to radians, forget the dot and get yourself into a non-floating point mess, even though pi is floating point...
        if Q==[]:
            q=self.getdata('Q')#1-by-N
            psi=self.getdata('PSI')#1-by-N
        else:
            #externally supplied data
            q=Q
            psi=PSI
        R1,R2,rot=Rset[:,0],Rset[:,1],Rset[:,2]
        NR=prod(shape(R1))
        if NR==1:
            #option 1:
            sda=sin((psi-rot)*d_to_r)
            cda=cos((psi-rot)*d_to_r)
            r=sqrt(R1**2*sda**2+R2**2*cda**2)
            qr=q*r
            Fell=3*(sin(qr)-qr*cos(qr))/(qr**3)
            ##quicker? no, 20% slower:
            #Fell=3*(
            #        sin(q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #            +R2**2*cos((psi-rot)*d_to_r)**2))
            #        -q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #            +R2**2*cos((psi-rot)*d_to_r)**2)
            #        *cos(q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #            +R2**2*cos((psi-rot)*d_to_r)**2)))/((q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #                +R2**2*cos((psi-rot)*d_to_r)**2))**3)
        else: #calculate a series
            Fell=zeros([NR,prod(shape(q))])
            for Ri in range(size(R1)):
                sda=sin((psi-rot[Ri])*d_to_r)
                cda=cos((psi-rot[Ri])*d_to_r)
                r=sqrt(R1[Ri]**2*sda**2+R2[Ri]**2*cda**2)
                qr=q*r
                Fell[Ri,:]=3*(sin(qr)-qr*cos(qr))/(qr**3)

        return Fell #this will be n-by-len(q) array

    def _passthrough(self,In):
        """A passthrough mechanism returning the input unchanged"""
        return In

    def reshape_fitdata(self):
        """This ensures that q, I, PSI and E are in 1-by-n shape"""
        for key in self.fitdata.keys():
            self.fitdata[key]=reshape(self.fitdata[key],(1,prod(shape(self.fitdata[key]))))

    def clip_dataset(self):
        """If q and/or psi limits are supplied in self.parameters,
        clips the dataset to within the supplied limits. Copies data to
        self.fitdata if no limits are set.
        """
        qlims=self.getpar('qlims')
        psilims=self.getpar('psilims')
        dataset=self.getdata(dataset='original')
        
        validbools=isfinite(dataset['Q'])
        # Optional masking of negative intensity
        if self.getpar('MaskNegI'):
            validbools=validbools*(I >= 0)
        if (qlims==[])and(psilims==[]):
            #q limits not set, simply copy dataset to fitdata
            validbools=validbools
        if (not(qlims==[])): #and qlims is implicitly set
            validbools = validbools*(dataset['Q']>numpy.min(qlims))&(dataset['Q']<=numpy.max(qlims)) #excluding the lower q limit may prevent q=0 from appearing
        if (not(psilims==[])): #we assume here that we have a dataset ['PSI']
            validbools = validbools*(dataset['PSI']>numpy.min(psilims))&(dataset['PSI']<=numpy.max(psilims)) #excluding the lower q limit may prevent q=0 from appearing

        for key in dataset.keys():
            dsk=dataset[key][validbools]
            #hey, this works!:
            self.setdata(**{key:dsk,'dataset':'fit'})
            #old version was direct addressing, which is to be discouraged to encourage flexibility in data storage
            #self.fitdata[key]=dataset[key][validbools]

        #self.fitdata=fitdata
        
    def getpar(self,parname=[]):
        """Gets the value of a parameter, so simple it is probably
        superfluous.
        """
        if parname==[]:
            return self.parameters
        else:
            return self.parameters[parname]

    def getdata(self,parname=[],dataset='fit'):
        """Gets the values of a dataset, retrieves from fitdata (clipped)
        by default. If the original data is wanted,
        use ``dataset = 'original'`` as *\*\*kwarg*.
        """
        if (parname==[]):
            if (dataset=='fit'):
                return self.fitdata
            else: 
                return self.dataset
        else:
            if (parname in self.fitdata.keys())and(dataset=='fit'):
                return self.fitdata[parname]
            else:
                return self.dataset[parname]

    def getresult(self,parname=[],VariableNumber=0):
        """Returns the specified entry from common result container."""
        if parname==[]:
            return self.result[VariableNumber]
        else:
            return self.result[VariableNumber][parname]

    def setresult(self,**kwargs):
        """Sets the supplied keyword-value pairs to the result. These can be
        arbitrary. Varnum is the sequence number of the variable for which
        data is stored. Default is set to 0, which is where the output of the
        MC routine is put before histogramming. The Histogramming procedure
        may populate more variables if so needed.
        """
        if 'VariableNumber' in kwargs.keys():
            varnum=kwargs['VariableNumber']
        else:
            varnum=0

        while len(self.result)<(varnum+1):
            #make sure there is a dictionary in the location we want to save the result
            self.result.append(dict())
        
        rdict=self.result[varnum]

        for kw in kwargs:
            rdict[kw]=kwargs[kw]

    def setpar(self,**kwargs):
        """Sets the supplied parameters given in keyword-value pairs for known
        setting keywords (unknown key-value pairs are skipped).
        If a supplied parameter is one of the function names, it is stored in
        the self.functions dict.
        """
        for kw in kwargs:
            if kw in self.parameters.keys():
                self.parameters[kw]=kwargs[kw]
            else:
                pass #no need to store unknown keywords.

    def setdata(self,**kwargs):
        """Sets the supplied data in the proper location. Optional argument
        *dataset* can be set to ``fit`` or ``original`` to define which
        dataset is set. Default is ``original``.
        """
        datasetlist=list(['Q','I','PSI','IERR']) #list of valid things
        if ('dataset' in kwargs):
            dataset=kwargs['dataset'].lower()
        else:
            dataset='original'
        if not ( dataset in ('fit','original')):
            dataset='original'

        if dataset=='original':
            for kw in kwargs:
                if kw in datasetlist:
                    self.dataset[kw]=kwargs[kw]
                else:
                    pass #do not store non-dataset values.
        else: #we want to store to fitdata: a clipped dataset
            for kw in kwargs:
                if kw in datasetlist:
                    self.fitdata[kw]=kwargs[kw]
                else:
                    pass #do not store non-dataset values.

    def EllBounds_2D(self):
        """This function will take the q and psi input bounds and outputs
        properly formatted two-element size bounds for ellipsoids. Ellipsoids
        are defined by their equatorial radius, meridional radius and axis
        misalignment (default -45 to 45 degrees in PSI).
        """
        Bounds=self.getpar('Bounds')
        q=self.getdata('Q')
        qBounds = array([pi/numpy.max(q),pi/numpy.min((abs(numpy.min(q)),abs(numpy.min(diff(q)))))]) # reasonable, but not necessarily correct, parameters
        if len(Bounds)==0:
            print 'Bounds not provided, so set related to minimum q or minimum q step and maximum q. Lower and upper bounds are {0} and {1}'.format(qBounds[0],qBounds[1])
            Bounds=numpy.array([qBounds[0],qBounds[1],qBounds[0],qBounds[1],-45,45])
        elif len(Bounds)==6:
            pass
            #print 'Bounds provided, set to {} and {}'.format(Bounds[0],Bounds[1])
        else:
            print 'Wrong number of Bounds provided, defaulting to {} and {} for radii, -45, 45 for misalignment'.format(qBounds[0],qBounds[1])
            Bounds=numpy.array([qBounds[0],qBounds[1],qBounds[0],qBounds[1],-45,45])
        Bounds=numpy.array([numpy.min(Bounds[0:2]),numpy.max(Bounds[0:2]),numpy.min(Bounds[2:4]),numpy.max(Bounds[2:4]),numpy.min(Bounds[4:6]),numpy.max(Bounds[4:6])])

        self.setpar(Bounds=Bounds)

    def SphBounds(self):
        """This function will take the q and input bounds and outputs properly
        formatted two-element size bounds.
        """
        Bounds=self.getpar('Bounds')
        q=self.getdata('Q')
        qBounds = array([pi/numpy.max(q),pi/numpy.min((abs(numpy.min(q)),abs(numpy.min(diff(q)))))]) # reasonable, but not necessarily correct, parameters
        if len(Bounds)==0:
            print 'Bounds not provided, so set related to minimum q or minimum q step and maximum q. Lower and upper bounds are {0} and {1}'.format(qBounds[0],qBounds[1])
            Bounds=qBounds
        elif len(Bounds)==1:
            print 'Only one bound provided, assuming it denotes the maximum. Lower and upper bounds are set to {0} and {1}'.format(qBounds[0],Bounds[1])
            Bounds=numpy.array([qBounds[0],Bounds])
        elif len(Bounds)==2:
            pass
            #print 'Bounds provided, set to {} and {}'.format(Bounds[0],Bounds[1])
        else:
            print 'Wrong number of Bounds provided, defaulting to {} and {}'.format(qBounds[0],qBounds[1])
            Bounds=qbounds
        Bounds=numpy.array([numpy.min(Bounds),numpy.max(Bounds)])

        self.setpar(Bounds=Bounds)

    def check_parameters(self):
        """Checks for the parameters, for example to make sure
        histbins is defined for all, or to check if all parameters fall
        within their limits.
        For now, all I need is a check that Histbins is a 1D vector
        with n values, where n is the number of parameters specifying
        a shape.
        """
        #testR=self.functions['RAND']()
        testR=self.getfunction('RAND')()
        NRval=prod(shape(testR))
        Histbins=self.getpar('Histbins')
        if not(isinstance(Histbins,list)): #meaning it will be a string
            HB=list()
            for ri in range(NRval):
                HB.append(int(Histbins))
            self.setpar(Histbins=HB)
        elif len(Histbins)<NRval:
            #histbins is a list but not of the right length
            while len(Histscale)<NRval:
                Histbins.append(Histbins[0])
            self.setpar(Histbins=Histbins)
        #now check histscale
        Histscale=self.getpar('Histscale')
        if not(isinstance(Histscale,list)): #meaning it will be a string
            HS=list()
            for ri in range(NRval):
                HS.append(Histscale) #repeat until we have enough
            #replace histscale
            self.setpar(Histscale=HS)
        elif len(Histscale)<NRval:
            #histscale is a list but not of the right length
            while len(Histscale)<NRval:
                Histscale.append(Histscale[0])
            self.setpar(Histscale=Histscale)

    ##########################################################################
    ####################### optimisation functions ###########################
    ##########################################################################

    def Analyse(self):
        """This function runs the Monte Carlo optimisation a multitude
        (*Nreps*) of times. If convergence is not achieved, it will try again
        for a maximum of *Maxntry* attempts.
        """
        #get data
        q=self.getdata('Q')
        I=self.getdata('I')
        E=self.getdata('IERR')
        #get settings
        Priors=self.getpar('Priors')
        Prior=self.getpar('Prior')
        Ncontrib=self.getpar('Ncontrib')
        Nreps=self.getpar('Nreps')
        Convcrit=self.getpar('Convcrit')
        Maxntry=self.getpar('Maxntry')
        #find out how many values a shape is defined by:
        #testR=self.functions['RAND']()
        testR=self.getfunction('RAND')()
        NRval=prod(shape((testR)))

        Rrep = zeros([Ncontrib,NRval,Nreps]) 
        #DEBUG:
        #print 'Rrep: {}'.format(shape(Rrep))
        Niters = zeros([Nreps])
        Irep = zeros([1,prod(shape(I)),Nreps])
        bignow = time.time() #for time estimation and reporting

        #This is the loop that repeats the MC optimization Nreps times, after which we can calculate an uncertainty on the results.
        priorsflag=False
        for nr in arange(0,Nreps):
            if ((Prior==[])and(Priors!=[]))or(priorsflag==True):
                priorsflag=True #this flag needs to be set as prior will be set after the first pass
                self.setpar(Prior=Priors[:,:,nr%size(Priors,2)])
            nt = 0 #keep track of how many failed attempts there have been 
            # do that MC thing! 
            ConVal=inf
            while ConVal>Convcrit:
                #retry in the case we were unlucky in reaching convergence within Maxiter.
                nt+=1
                Rrep[:,:,nr],Irep[:,:,nr],ConVal,Details = self.MCFit(OutputI=True,OutputDetails=True)
                if nt>Maxntry:
                    #this is not a coincidence. We have now tried Maxntry+2 times
                    print "could not reach optimization criterion within {0} attempts, exiting...".format(Maxntry+2)
                    return
            Niters[nr] = Details['Niter'] #keep track of how many iterations were needed to reach convergence

            biglap = time.time() #time management
            # in minutes:
            tottime = (biglap-bignow)/60. #total elapsed time
            avetime = (tottime/(nr+1)) #average time per MC optimization
            remtime = (avetime*Nreps-tottime) #estimated remaining time
            print "\t*finished optimization number {0} of {1} \r\n\t*total elapsed time: {2} minutes \r\n\t*average time per optimization {3} minutes \r\n\t*total time remaining {4} minutes".format(nr+1,Nreps,tottime,avetime,remtime)
        
        #at this point, all MC optimizations have been completed and we can process all Nreps results.
        Imean = numpy.mean(Irep,axis=2) #mean fitted intensity
        Istd = numpy.std(Irep,axis=2) #standard deviation on the fitted intensity, usually not plotted for clarity
        # store in output dict
        self.setresult(**{
            'Rrep':Rrep,
            'Imean':Imean,
            'Istd':Istd,
            'Qfit':q,# can be obtained from self.fitdata 
            'Niter':numpy.mean(Niters)}) #average number of iterations for all repetitions

    def MCFit(self,OutputI=False,OutputDetails=False,OutputIterations=False):
        """
        Object-oriented, shape-flexible core of the Monte Carlo procedure.
        Takes optional arguments:

        *OutputI*:
            Returns the fitted intensity besides the result

        *OutputDetails*:
            Details of the fitting procedure, number of iterations and so on

        *OutputIterations*:
            Returns the result on every successful iteration step, useful for
            visualising the entire Monte Carlo optimisation procedure for
            presentations.
        """
        #load dataset
        q=self.getdata('Q')
        I=self.getdata('I')
        E=self.getdata('IERR')
        #load parameters
        Ncontrib=self.getpar('Ncontrib')
        Bounds=self.getpar('Bounds')
        Convcrit=self.getpar('Convcrit')
        Rpfactor=self.getpar('Rpfactor')
        Maxiter=self.getpar('Maxiter')
        MaskNegI=self.getpar('MaskNegI')
        StartFromMin=self.getpar('StartFromMin')
        Memsave=self.getpar('Memsave')
        Prior=self.getpar('Prior')


        #find out how many values a shape is defined by:
        Randfunc=self.getfunction('RAND')
        FFfunc=self.getfunction('FF')
        VOLfunc=self.getfunction('VOL')
        SMEARfunc=self.getfunction('SMEAR')
        testR=Randfunc()
        NRval=prod(shape(testR))

        Rset=numpy.zeros((Ncontrib,NRval))


        # Intialise variables
        FFset = []
        Vset = []
        Niter = 0
        Conval = inf
        Details = dict()
        Ri = 0 #index of sphere to change. We'll sequentially change spheres, which is perfectly random since they are in random order.
        
        #generate initial set of spheres
        if size(Prior)==0:
            if StartFromMin:
                for Rvi in range(NRval): #minimum bound for each value
                    if numpy.min(Bounds[Rvi:Rvi+2])==0:
                        mb=pi/numpy.max(q)
                    else:
                        mb=numpy.min(Bounds[Rvi:Rvi+2])
                    Rset[:,Rvi]=numpy.ones(Ncontrib)[:]*mb/2.
            else:
                Rset=Randfunc(Ncontrib)
        elif (size(Prior,0)!=0)&(size(Ncontrib)==0):
            Ncontrib=size(Prior,0)
            Rset=Prior
        elif size(Prior,0)==Ncontrib:
            Rset=Prior
        elif size(Prior,0)<Ncontrib:
            print "size of prior is smaller than Ncontrib. duplicating random prior values"
            #while size(Prior)<Nsph:
            Addi=numpy.random.randint(size(Prior,0),size=Ncontrib-size(Prior,0))
            Rset=concatenate((Prior,Prior[Addi,:]))
            print "size now:", size(Rset)
        elif size(Prior,0)>Ncontrib:
            print "Size of prior is larger than Ncontrib. removing random prior values"
            Remi=numpy.random.randint(size(Prior,0),size=Ncontrib) #remaining choices
            Rset=Prior[Remi,:]
            print "size now:", size(Rset,0)
        
        if Memsave==False:
            #calculate their form factors
            FFset=FFfunc(Rset)
            Vset=VOLfunc(Rset,Rpfactor)
            #Vset=(4.0/3*pi)*Rset**(3*Rpfactor)
            #calculate the intensities
            Iset=FFset**2*(Vset+0*FFset)**2 #a set of intensities
            Vst=sum(Vset**2) # total volume squared
            It=sum(Iset,0) # the total intensity - eq. (1)
            It=reshape(It,(1,prod(shape(It)))) # reshaped to match I and q
        else:
            #calculate intensity in memory saving mode:
            #calculate volume for entire set, this does not take much space   
            Vset=VOLfunc(Rset,Rpfactor)

            FFset=FFfunc(Rset[0,:][newaxis,:])
            It=FFset**2*(Vset[0]+0*FFset)**2 #a set of intensities
            for ri in arange(1,Ncontrib):
                #calculate their form factors
                FFset=FFfunc(Rset[ri,:][newaxis,:])
                It=It+FFset**2*(Vset[ri]+0*FFset)**2 #a set of intensities
            Vst=sum(Vset**2) # total volume squared
            It=reshape(It,(1,prod(shape(It)))) # reshaped to match I and q

        # Optimize the intensities and calculate convergence criterium
        #SMEAR function goes here
        It=SMEARfunc(It)
        Sci = numpy.max(I)/numpy.max(It) #initial guess for the scaling factor.
        Bgi = numpy.min(I)
        Sc,Conval=self._Iopt(I,It/Vst,E,numpy.array([Sci,Bgi]),ver=1) 
        Sc,Conval=self._Iopt(I,It/Vst,E,Sc) # reoptimize with V2, there might be a slight discrepancy in the residual definitions of V1 and V2 which would prevent optimization.
        #print "Initial conval V1",Conval1
        print "Initial Chi-squared value",Conval

        if OutputIterations:
            # Output each iteration, starting with number 0. Iterations will be stored 
            # in Details['Itersph'], Details['IterIfit'], Details['IterConval'], 
            # Details['IterSc'] and Details['IterPriorUnaccepted'] listing the 
            # unaccepted number of moves before the recorded accepted move.
            Details['Itersph']=Rset[:,newaxis] #new iterations will (have to) be appended to this, cannot be zero-padded due to size constraints
            Details['IterIfit']=(It/Vst*Sc[0]+Sc[1])[:,newaxis] #ibid.
            Details['IterConVal']=Conval[newaxis]
            Details['IterSc']=Sc[:,newaxis]
            Details['IterPriorUnaccepted']=numpy.array(0)[newaxis]

        #start the MC procedure
        Now=time.time()
        Nmoves=0 #tracking the number of moves
        Nnotaccepted=0
        while (Conval>Convcrit) &(Niter<Maxiter):
            Rt=Randfunc()
            Ft=FFfunc(Rt)
            Vtt=VOLfunc(Rt,Rpfactor)
            Itt=(Ft**2*Vtt**2)
            # Calculate new total intensity
            if Memsave==False:
                Itest=(It-Iset[Ri,:]+Itt) # we do subtractions and additions, which give us another factor 2 improvement in speed over summation and is much more scalable
            else:
                Fo=FFfunc(Rset[Ri,:][newaxis,:])
                Io=(Fo**2*Vset[Ri]**2)
                Itest=(It-Io+Itt)

            #SMEAR function goes here
            Itest=SMEARfunc(Itest)
            Vstest = (sqrt(Vst)-Vset[Ri])**2+Vtt**2
            # optimize intensity and calculate convergence criterium
            Sct,Convalt = self._Iopt(I,Itest/Vstest,E,Sc) # using version two here for a >10 times speed improvement
            # test if the radius change is an improvement:
            if Convalt<Conval: # it's better
                if Memsave:
                    Rset[Ri,:],It,Vset[Ri],Vst,Sc,Conval=(Rt,Itest,Vtt,Vstest,Sct,Convalt)
                else:
                    Rset[Ri,:],Iset[Ri,:],It,Vset[Ri],Vst,Sc,Conval=(Rt,Itt,Itest,Vtt,Vstest,Sct,Convalt)
                print "Improvement in iteration number %i, Chi-squared value %f of %f\r" %(Niter,Conval,Convcrit),
                Nmoves+=1
                if OutputIterations:
                    # output each iteration, starting with number 0. 
                    # Iterations will be stored in Details['Itersph'], Details['IterIfit'], 
                    # Details['IterConval'], Details['IterSc'] and 
                    # Details['IterPriorUnaccepted'] listing the unaccepted 
                    # number of moves before the recorded accepted move.
                    Details['Itersph']=concatenate((Details['Itersph'],Rset[:,:,newaxis]),axis=1) #new iterations will (have to) be appended to this, cannot be zero-padded due to size constraints
                    Details['IterIfit']=concatenate((Details['IterIfit'],(Itest/Vstest*Sct[0]+Sct[1])[:,newaxis]),axis=1) #ibid.
                    Details['IterConVal']=concatenate((Details['IterConVal'],numpy.array(Convalt)[newaxis]))
                    Details['IterSc']=concatenate((Details['IterSc'],Sct[:,newaxis]),axis=1)
                    Details['IterPriorUnaccepted']=concatenate((Details['IterPriorUnaccepted'],numpy.array(Nnotaccepted)[newaxis]))
                Nnotaccepted=-1
            # else nothing to do
            Ri+=1 # move to next sphere in list
            Ri=Ri%(Ncontrib) # loop if last sphere
            Nnotaccepted+=1 # number of non-accepted moves, resets to zero after accepted move.
            Niter+=1 # add one to the iteration number           
        if Niter>=Maxiter:
            print "exited due to max. number of iterations (%i) reached" %(Niter)
        else:
            print "Normal exit"
        print "Number of iterations per second",Niter/(time.time()-Now+0.001) #the +0.001 seems necessary to prevent a divide by zero error on some Windows systems.   
        print "Number of valid moves",Nmoves
        print "final Chi-squared value %f" %(Conval)
        Details['Niter']=Niter
        Details['Nmoves']=Nmoves
        Details['elapsed']=(time.time()-Now+0.001)

        #Ifinal=sum(Iset,0)/sum(Vset**2)
        Ifinal=It/sum(Vset**2)
        Ifinal=reshape(Ifinal,(1,prod(shape(Ifinal))))
        #SMEAR function goes here
        Ifinal=SMEARfunc(Ifinal)
        Sc,Conval=self._Iopt(I,Ifinal,E,Sc)    
        if OutputI:
            if OutputDetails:
                #DEBUG:
                #print 'Rset: {}, I: {}, Conval: {}'.format(shape(Rset),shape((Ifinal*Sc[0]+Sc[1])),shape(Conval))

                return Rset,(Ifinal*Sc[0]+Sc[1]),Conval,Details
            else:
                return Rset,(Ifinal*Sc[0]+Sc[1]),Conval
        else:
            if OutputDetails:
                return Rset,Conval,Details # ifinal cannot be output with variable length intensity outputs (in case of masked negative intensities or q limits)
            else:
                return Rset,Conval # ifinal cannot be output with variable length intensity outputs (in case of masked negative intensities or q limits)

    ############################################################################
    #################### Post-optimisation functions ###########################
    ############################################################################

    def Histogram(self):
        """
        Takes the *Rrep* result from the :py:meth:`McSAS.Analyse` function
        and calculates the corresponding volume- and number fractions for each
        contribution as well as the minimum observability limits. It will
        subsequently bin the result across the range for histogramming purposes.

        While the volume-weighted distribution will be in absolute units
        (providing volume fractions of material within a given size range),
        the number distributions have been normalized to 1.
        
        Output a list of dictionaries with one dictionary per shape parameter:

            *VariableNumber*: int
                Shape parameter index. e.g. an ellipsoid has 3:
                width, height and orientation
            *Hx*: array
                Histogram bin left edge position (x-axis in histogram)
            *Hmid*: array
                Center positions for the size histogram bins
                (x-axis in histogram, used for errorbars)
            *Hwidth*: array
                Histogram bin width (x-axis in histogram,
                defines bar plot bar widths)
            *Hmean*: array
                Volume-weighted particle size distribution values for
                all *Nreps* results (y-axis bar height)
            *Hnmean*: array
                Number-weighted analogue of the above *Hmean*
            *Hy*: size (Histbins x Nreps) array
                Volume-weighted particle size distribution bin values for each
                MC fit repetition (the mean of which is *Hmean*, and the sample
                standard deviation of which is *Hstd*)
            *Hny*: size (Histbins x Nreps) array
                Number-weighted particle size distribution bin values
                for each MC fit repetition
            *Hstd*: array
                Standard deviations of the corresponding volume-weighted size
                distribution bins, calculated from *Nreps* repetitions of the
                MCfit_sph() function
            *Hnstd*: array
                Standard deviation for the number-weigthed distribution
            *Vf*: size (Ncontrib x Nreps) array
                Volume fractions for each of Ncontrib contributions 
                in each of Nreps iterations
            *Nf*: size (Ncontrib x Nreps) array
                Number fraction for each contribution
            *Vft*: size (Nreps) array
                Total scatterer volume fraction for each of the *Nreps*
                iterations
            *Nft*: size (Nreps) array
                Total number fraction 
            *vfmin*: size (Nsph x Nreps) array
                minimum required volume fraction for each contribution to
                become statistically significant.
            *nfmin*: size (Nsph x Nreps) array
                number-weighted analogue to *vfmin*
            *vfminbins*: size (Hmid) array 
                array with the minimum required volume fraction per bin to
                become statistically significant. Used to display minimum
                required level in histogram.
            *nfminbins*: size (Hmid) array
                number-weighted analogue to *vfminbins*
            *Screps*: size (2 x Nreps) array
                Scaling and background values for each repetition. Used to
                display background level in data and fit plot.
        """
        #get settings
        #set the bin edges for our radius bins either based on a linear division or on a logarithmic division of radii.
        Ncontrib=self.getpar('Ncontrib')
        Nreps=self.getpar('Nreps')
        Rpfactor=self.getpar('Rpfactor')
        Memsave=self.getpar('Memsave')
        drhosqr=self.getpar('drhosqr')
        Rrep=self.getresult('Rrep')
        Histbins=self.getpar('Histbins')
        Histscale=self.getpar('Histscale')
        Bounds=self.getpar('Bounds')

        #ov = zeros(shape(Rrep)) #observability
        Vf = zeros((Ncontrib,Nreps)) #volume fraction for each contribution
        Nf = zeros((Ncontrib,Nreps)) #number fraction for each contribution
        qm = zeros((Ncontrib,Nreps)) #volume fraction for each contribution
        vfmin = zeros((Ncontrib,Nreps)) #volume fraction for each contribution
        nfmin = zeros((Ncontrib,Nreps)) #number fraction for each contribution
        Vft = zeros([Nreps]) #total volume fractions
        Nft = zeros([Nreps]) #total number 
        Screps = zeros([2,Nreps]) #Intensity scaling factors for matching to the experimental scattering pattern (Amplitude A and flat background term b, defined in the paper)

        #functions!
        Randfunc=self.getfunction('RAND')
        FFfunc=self.getfunction('FF')
        VOLfunc=self.getfunction('VOL')
        SMEARfunc=self.getfunction('SMEAR')

        #data!
        q=self.getdata('Q')
        I=self.getdata('I')
        E=self.getdata('IERR')

        #loop over each repetition
        for ri in range(Nreps):
            Rset = Rrep[:,:,ri] #the single set of R for this calculation




            Vset = VOLfunc(Rset,Rpfactor) #compensated volume for each sphere in the set
            if Memsave==False:
                FFset = FFfunc(Rset) #Form factors, all normalized to 1 at q=0.
                # Calculate the intensities
                Iset = FFset**2*(Vset+0*FFset)**2 # Intensity for each contribution as used in the MC calculation
                It = sum(Iset,0) # the total intensity of the scattering pattern
            else:
                FFset=FFfunc(Rset[0,:][newaxis,:])
                It=FFset**2*(Vset[0]+0*FFset)**2 #a set of intensities
                for Rr in arange(1,Ncontrib):
                    #calculate their form factors
                    FFset=FFfunc(Rset[Rr,:][newaxis,:])
                    It=It+FFset**2*(Vset[Rr]+0*FFset)**2 #a set of intensities

            Vst = sum(Vset**2) # total compensated volume squared 
            It=reshape(It,(1,prod(shape(It))))
            It=SMEARfunc(It)
            
            # Now for each sphere, calculate its volume fraction (p_c compensated):
            Vsa=VOLfunc(Rset,Rpfactor) #compensated volume for each sphere in the set Vsa = 4./3*pi*Rset**(3*Rpfactor)
            # And the real particle volume:
            Vpa=VOLfunc(Rset,Rpfactor=1.) #compensated volume for each sphere in the set Vsa = 4./3*pi*Rset**(3*Rpfactor)

            Sci = numpy.max(I)/numpy.max(It) #initial guess for the scaling factor.
            Bgi = numpy.min(I)
            Sc,Cv = self._Iopt(I,It,E,[Sci,Bgi]) #optimize scaling and background for this repetition
            Screps[:,ri]=Sc #scaling and background for this repetition.
            Vf[:,ri] = (Sc[0]*Vsa**2/(Vpa*drhosqr)).flatten() # a set of volume fractions
            Vft[ri] = sum(Vf[:,ri]) # total volume 
            Nf[:,ri] = Vf[:,ri]/(Vpa.flatten()) #
            Nft[ri] = sum(Nf[:,ri]) # total number
            for isi in range(Ncontrib): #For each sphere
                #ov[isi,ri] = (Iset[isi,:]/(It)).max() #calculate the observability (the maximum contribution for that sphere to the total scattering pattern) NOTE: no need to compensate for p_c here, we work with volume fraction later which is compensated by default. additionally, we actually do not use this value.
                if Memsave:
                    FFset=FFfunc(Rset[isi,:][newaxis,:])
                    Ir=FFset**2*(Vset[isi]+0*FFset)**2
                    qmi = numpy.argmax(Ir.flatten()/It.flatten()) #determine where this maximum observability is of contribution isi (index)
                    qm[isi,ri] = q[0,qmi] #point where the contribution of isi is maximum
                    vfmin[isi,ri] = numpy.min(E*Vf[isi,ri]/(Sc[0]*Ir))
                    nfmin[isi,ri] = vfmin[isi,ri]/Vpa[isi]
                else:
                    qmi = numpy.argmax(Iset[isi,:].flatten()/It.flatten()) #determine where this maximum observability is of contribution isi (index)
                    qm[isi,ri] = q[0,qmi] #point where the contribution of isi is maximum
                    vfmin[isi,ri] = numpy.min(E*Vf[isi,ri]/(Sc[0]*Iset[isi,:]))
                    nfmin[isi,ri] = vfmin[isi,ri]/Vpa[isi]
                #close approximation:
                #vfmin[isi,ri] = (E[qmi]*Vf[isi,ri]/(Sc[0]*Iset[isi,qmi]))
                #or more precice but slower:
            Nf[:,ri]=Nf[:,ri]/Nft[ri]
            nfmin[:,ri]=nfmin[:,ri]/Nft[ri]

        #now we histogram over each variable
        #for each variable parameter we define, we need to histogram separately. 
        for vari in range(prod(shape(Histbins))):
            # Now bin whilst keeping track of which contribution ends up in which bin:
            #set bin edge locations
            if Histscale[vari] == 'lin':
                Hx = linspace(Bounds[0+2*vari],Bounds[1+2*vari],Histbins[vari]+1) #Hx contains the Histbins+1 bin edges, or class limits.
            else:
                Hx = 10**(linspace(log10(Bounds[0+2*vari]),log10(Bounds[1+2*vari]),Histbins[vari]+1))
            Hy = zeros([Histbins[vari],Nreps]) #total volume fraction contribution in a bin
            Hny = zeros([Histbins[vari],Nreps]) #total number fraction contribution in a bin
            vfminbin = zeros([Histbins[vari],Nreps]) #minimum required number of contributions /in a bin/ to make a measurable impact
            nfminbin = zeros([Histbins[vari],Nreps]) #minimum required number of contributions /in a bin/ to make a measurable impact
            Hmid = zeros(Histbins[vari])
            vfminbins = zeros(Histbins[vari])
            nfminbins = zeros(Histbins[vari])

            for ri in range(Nreps):
                
                Rset = Rrep[:,vari,ri] #the single set of R for this calculation

                for bini in range(Histbins[vari]):
                    findi = ((Rset>=Hx[bini])*(Rset<Hx[bini+1])) #indexing which contributions fall into the radius bin
                    #print 'findi: {} Vf: {}'.format(shape(findi),shape(Vf))
                    #findi = findi[:,0]
                    Hy[bini,ri] = sum(Vf[findi,ri]) #y contains the volume fraction for that radius bin
                    Hny[bini,ri] = sum(Nf[findi,ri]) #y contains the volume fraction for that radius bin
                    if sum(findi)==0:
                        vfminbin[bini,ri] = 0
                        nfminbin[bini,ri] = 0
                    else:
                        vfminbin[bini,ri] = numpy.max(vfmin[findi,ri])
                        vfminbin[bini,ri] = numpy.mean(vfmin[findi,ri])
                        nfminbin[bini,ri] = numpy.max(nfmin[findi,ri])
                        nfminbin[bini,ri] = numpy.mean(nfmin[findi,ri])
                    if isnan(Hy[bini,ri]):
                        Hy[bini,ri] = 0.
                        Hny[bini,ri] = 0.
            for bini in range(Histbins[vari]):
                Hmid[bini] = numpy.mean(Hx[bini:bini+2])
                vb = vfminbin[bini,:]
                vfminbins[bini] = numpy.max(vb[vb<inf])
                nb = nfminbin[bini,:]
                nfminbins[bini] = numpy.max(nb[vb<inf])
            Hmean = numpy.mean(Hy,axis=1)
            Hnmean = numpy.mean(Hny,axis=1)
            Hstd = numpy.std(Hy,axis=1)
            Hnstd = numpy.std(Hny,axis=1)
            self.setresult(**{
                'VariableNumber':vari, #this line will place the results in the dict at self.results[vari]
                'Hx':Hx,
                'Hmid':Hmid,
                'Hwidth':diff(Hx),
                'Hy':Hy,
                'Hny':Hny,
                'Hmean':Hmean,
                'Hstd':Hstd,
                'Hnmean':Hnmean,
                'Hnstd':Hnstd,
                'vfminbins':vfminbins,
                'vfmin':vfmin,
                'Vf':Vf,
                'Vft':Vft,
                'nfminbins':nfminbins,
                'nfmin':nfmin,
                'Nf':Nf,
                'Nft':Nft,
                'Screps':Screps})

    def TwoDGenI(self):
        """
        this function is optionally run after the histogram procedure for anisotropic images, and will calculate the MC fit intensity in imageform
        """
        Result=self.getresult()
        #load original dataset
        q=self.getdata('Q',dataset='original')
        I=self.getdata('I',dataset='original')
        E=self.getdata('IERR',dataset='original')
        PSI=self.getdata('PSI',dataset='original')
        #we need to recalculate the result in two dimensions
        kansas=shape(q) #we will return to this shape
        q=reshape(q,[1,-1]) #flatten
        I=reshape(I,[1,-1]) #flatten
        E=reshape(E,[1,-1]) #flatten
        PSI=reshape(PSI,[1,-1]) #flatten

        Randfunc=self.getfunction('RAND')
        FFfunc=self.getfunction('FF')
        VOLfunc=self.getfunction('VOL')
        SMEARfunc=self.getfunction('SMEAR')
        print 'Recalculating final result, this may take some time'
        #for each result
        Iave=zeros(shape(q))
        Nreps=self.getpar('Nreps')
        Rpfactor=self.getpar('Rpfactor')
        qlims=self.getpar('qlims')
        psilims=self.getpar('psilims')
        Memsave=self.getpar('Memsave')
        Ncontrib=self.getpar('Ncontrib')
        Screps=self.getresult('Screps')
        for nr in range(Nreps):
            print 'regenerating set {} of {}'.format(nr,Nreps)
            Rset=Result['Rrep'][:,:,nr]
            #calculate their form factors
            Vset=VOLfunc(Rset,Rpfactor)
            #Vset=(4.0/3*pi)*Rset**(3*Rpfactor)
            #calculate the intensities
            if Memsave==False:
                FFset = FFfunc(Rset,Q=q,PSI=PSI) #Form factors, all normalized to 1 at q=0.
                # Calculate the intensities
                Iset = FFset**2*(Vset+0*FFset)**2 # Intensity for each contribution as used in the MC calculation
                It = sum(Iset,0) # the total intensity of the scattering pattern
            else:
                FFset=FFfunc(Rset[0,:][newaxis,:],Q=q,PSI=PSI)
                It=FFset**2*(Vset[0]+0*FFset)**2 #a set of intensities
                for ri in arange(1,Ncontrib):
                    #calculate their form factors
                    FFset=FFfunc(Rset[ri,:][newaxis,:],Q=q,PSI=PSI)
                    It=It+FFset**2*(Vset[ri]+0*FFset)**2 #a set of intensities
            Vst=sum(Vset**2) # total volume squared
            It=reshape(It,(1,-1)) # reshaped to match I and q
            # Optimize the intensities and calculate convergence criterium
            #SMEAR function goes here
            It=SMEARfunc(It)
            Iave=Iave+It*Screps[0,nr]+Screps[1,nr] #add to average
        #print "Initial conval V1",Conval1
        Iave=Iave/Nreps
        #mask (lifted from clip_dataset)
        validbools=isfinite(q)
        # Optional masking of negative intensity
        if self.getpar('MaskNegI'):
            validbools=validbools*(I >= 0)
        if (qlims==[])and(psilims==[]):
            #q limits not set, simply copy dataset to fitdata
            validbools=validbools
        if (not(qlims==[])): #and qlims is implicitly set
            validbools = validbools*(q>numpy.min(qlims))&(q<=numpy.max(qlims)) #excluding the lower q limit may prevent q=0 from appearing
        if (not(psilims==[])): #we assume here that we have a dataset ['PSI']
            validbools = validbools*(PSI>numpy.min(psilims))&(PSI<=numpy.max(psilims)) #excluding the lower q limit may prevent q=0 from appearing
        Iave=Iave*validbools
        #shape back to imageform
        I2D=reshape(Iave,kansas)
        self.setresult(I2D=I2D)

    def CSVwrite(self,filename,*args,**kwargs):
        """
        This function writes a semicolon-separated csv file to [filename]
        containing an arbitrary number of output variables *\*args*. in case of
        variable length columns, empty fields will contain ''.

        Optional last argument is a keyword-value argument:
        VarableNumber=[integer], indicating which shape parameter it is
        intended to draw upon. VariableNumber can also be a list or array
        of integers, of the same length as the number of output variables
        *\*args* in which case each output variable is matched with a shape
        parameter index. Default is zero.

        Input arguments should be names of fields in *self.result*.
        For example::

            A.McCSV('hist.csv', 'Hx', 'Hwidth', 'Hmean', 'Hstd',
                    VariableNumber = 0)

        I.e. just stick on as many columns as you'd like. They will be
        flattened by default. A header with the result keyword names will be
        added.
        
        Existing files with the same filename will be overwritten by default.
        """
        vna=zeros(len(args),dtype=int)
        if 'VariableNumber' in kwargs:
            vni=kwargs['VariableNumber']
            if isinstance(vni,(list,ndarray)):
                if len(vni)!=len(args):
                    print('Error in CSVwrite, supplied list of variablenumbers does not have the length of 1 or the same length as the list of output variables.')
                    return
                for vi in range(len(args)):
                    vna[vi]=vni[vi]
            else:
                #single integer value
                vna=vna+vni
                
        #uses sprintf rather than csv for flexibility
        ncol=len(args)
        #make format string used for every line, don't need this
        #linestr=''
        #for coli in range(ncol):
        #    linestr=linestr+'{'+'};'
        #linestr=linestr[0:-1]+'\n' #strip the last semicolon, add a newline

        inlist=list()
        for argi in range(len(args)):
            inlist.append(self.getresult(args[argi],VariableNumber=vna[argi]).flatten())
        #find out the longest row
        nrow=0
        for argi in range(len(args)):
            nrow=numpy.max((nrow,len(inlist[argi])))

        #now we can open the file:
        fh=open(filename,'w')
        emptyfields=0
        #write header:
        linestr=''
        for coli in range(ncol):
            linestr=linestr+'{};'.format(args[coli])
        linestr=linestr[0:-1]+'\n'    
        fh.write(linestr)
        for rowi in range(nrow):
            linestr=''
            for coli in range(ncol):
                #print 'rowi {} coli {} len(args[coli]) {}'.format(rowi,coli,len(args[coli]))
                if len(inlist[coli])<=rowi: #we ran out of numbers for this arg
                    linestr=linestr+';' #add empty field
                    emptyfields+=1
                else:
                    linestr=linestr+'{};'.format(inlist[coli][rowi])
            linestr=linestr[0:-1]+'\n'

            fh.write(linestr)

        fh.close()
        print '{} lines written with {} columns per line, and {} empty fields'.format(rowi,ncol,emptyfields)

    def Plot(self,AxisMargin=0.3):
        """
        This function plots the output of the Monte-Carlo procedure in two
        windows, with the left window the measured signal versus the fitted
        intensity (on double-log scale), and the righthand window the size
        distribution.
        """
        import matplotlib.font_manager as fm
        def setaxis(ah):
            import matplotlib.font_manager as fm
            plotfont = fm.FontProperties(
                        #this only works for macs, doesn't it?
                        #family = 'Courier New Bold', fname = '/Library/Fonts/Courier New Bold.ttf')
                        family = 'Arial')
            textfont = fm.FontProperties( #Baskerville.ttc does not work when saving to eps
                        #family = 'Times New Roman', fname = '/Library/Fonts/Times New Roman.ttf')
                        family = 'Times')
            "sets the axes parameters. axtyp can be one of 'q' or 'R'"
            #setaxis font and ticks
            ah.set_yticklabels(ah.get_yticks(), fontproperties = plotfont,size='large')
            ah.set_xticklabels(ah.get_xticks(), fontproperties = plotfont,size='large')
            ah.set_xlabel(ah.get_xlabel(), fontproperties=textfont,size='x-large')
            ah.set_ylabel(ah.get_ylabel(), fontproperties=textfont,size='x-large')
            #q_ax.set_yticklabels(q_ax.get_yticks(), fontproperties = plotfont)
            #q_ax.set_xticklabels(q_ax.get_xticks(), fontproperties = plotfont)
            #R_ax.spines['bottom'].set_color('black')
            ah.spines['bottom'].set_lw(2)
            ah.spines['top'].set_lw(2)
            ah.spines['left'].set_lw(2)
            ah.spines['right'].set_lw(2)
            ah.tick_params(axis='both',colors='black',width=2,which='major',direction='in',length=6)
            ah.tick_params(axis='x',colors='black',width=2,which='minor',direction='in',length=3)
            ah.tick_params(axis='y',colors='black',width=2,which='minor',direction='in',length=3)
            #q_ax.spines['bottom'].set_lw(2)
            #q_ax.spines['top'].set_lw(2)
            #q_ax.spines['left'].set_lw(2)
            #q_ax.spines['right'].set_lw(2)
            #q_ax.tick_params(axis='both',colors='black',width=2,which='major',direction='in',length=6)
            #q_ax.tick_params(axis='x',colors='black',width=2,which='minor',direction='in',length=3)
            #q_ax.tick_params(axis='y',colors='black',width=2,which='minor',direction='in',length=3)
            locs,labels = xticks()
            xticks(locs, map(lambda x: "%g" % x, locs))
            locs,labels = yticks()
            yticks(locs, map(lambda x: "%g" % x, locs))
            return ah

        #load parameters
        Histscale=self.getpar('Histscale')
        Histweight=self.getpar('Histweight')
        #load result
        Result=self.getresult()
        #check how many result plots we need to generate: maximum three.
        nhists=len(Histscale)

        #set plot font
        plotfont = fm.FontProperties(
                    size='large',
                    family = 'Arial')
        textfont = fm.FontProperties( #Baskerville.ttc does not work when saving to eps
                    size='large',
                    family = 'Times')
        #initialize figure and axes
        fig=figure(figsize=(7*(nhists+1),7),dpi=80,facecolor='w',edgecolor='k')
        #load original dataset
        q=self.getdata('Q',dataset='original')
        I=self.getdata('I',dataset='original')
        E=self.getdata('IERR',dataset='original')
        TwoDMode=False
        if ndim(q)>1:
            #2D data
            TwoDMode=True
            PSI=self.getdata('PSI',dataset='original')
            #we need to recalculate the result in two dimensions
            #done by TwoDGenI function
            I2D=self.getresult('I2D')
            Ishow=I.copy()
            #quadrant 1 and 4 are simulated data, 2 and 3 are measured data
            Ishow[(PSI>0)*(PSI<=90)]=I2D[(PSI>0)*(PSI<=90)]
            Ishow[(PSI>180)*(PSI<=270)]=I2D[(PSI>180)*(PSI<=270)]
            #xalimits=(-numpy.min(q[:,0]),numpy.max(q[:,-1]))
            #yalimits=(-numpy.min(q[0,:]),numpy.max(q[-1,:]))
            xmidi=int(round(size(q,1)/2))
            ymidi=int(round(size(q,0)/2))
            QX=numpy.array([-q[ymidi,0],q[ymidi,-1]])
            QY=numpy.array([-q[0,xmidi],q[-1,xmidi]])
            extent=(QX[0],QX[1],QY[0],QY[1])

            q_ax=fig.add_subplot(1,(nhists+1),1,axisbg=(.95,.95,.95),xlim=QX,ylim=QY,xlabel='q_x, 1/m',ylabel='q_y, 1_m')
            imshow(log10(Ishow),extent=extent,origin='lower')
            q_ax=setaxis(q_ax)
            colorbar()
        else:
            q_ax=fig.add_subplot(1,(nhists+1),1,axisbg=(.95,.95,.95),xlim=(numpy.min(q)*(1-AxisMargin),numpy.max(q)*(1+AxisMargin)),ylim=(numpy.min(I)*(1-AxisMargin),numpy.max(I)*(1+AxisMargin)),xscale='log',yscale='log',xlabel='q, 1/m',ylabel='I, 1/(m sr)')
            q_ax=setaxis(q_ax)
            errorbar(q,I,E,zorder=2,fmt='k.',ecolor='k',elinewidth=2,capsize=4,ms=5,label='Measured intensity',lw=2,solid_capstyle='round',solid_joinstyle='miter')
            grid(lw=2,color='black',alpha=.5,dashes=[1,6],dash_capstyle='round',zorder=-1)
            #xscale('log')
            #yscale('log')
            aq=sort(Result['Qfit'][0,:])
            aI=Result['Imean'][0,argsort(Result['Qfit'][0,:])]
            plot(aq,aI,'r-',lw=3,label='MC Fit intensity',zorder=4)
            plot(aq,numpy.mean(Result['Screps'][1,:])+0*aq,'g-',linewidth=3,label='MC Background level:\n\t ({0:03.3g})'.format(numpy.mean(Result['Screps'][1,:])),zorder=3)
            leg=legend(loc=1,fancybox=True,prop=textfont)
        title('Measured vs. Fitted intensity',fontproperties=textfont,size='x-large')

        R_ax=list()
        for histi in range(nhists):
            #get data:
            Hx=self.getresult(parname='Hx',VariableNumber=histi)
            Hmid=self.getresult(parname='Hmid',VariableNumber=histi)
            Hwidth=self.getresult(parname='Hwidth',VariableNumber=histi)
            if Histweight=='volume':
                Hmean=self.getresult(parname='Hmean',VariableNumber=histi)
                vfminbins=self.getresult(parname='vfminbins',VariableNumber=histi)
                Hstd=self.getresult(parname='Hstd',VariableNumber=histi)
            elif Histweight=='number':
                Hmean=self.getresult(parname='Hnmean',VariableNumber=histi)
                vfminbins=self.getresult(parname='nfminbins',VariableNumber=histi)
                Hstd=self.getresult(parname='Hnstd',VariableNumber=histi)
            else: 
                print 'Incorrect value for Histweight: should be either "volume" or "number"'

            #prep axes
            if Histscale[histi]=='log': #quick fix with the [0] reference. Needs fixing, this plotting function should be rewritten to support multiple variables.
                R_ax.append(fig.add_subplot(1,(nhists+1),histi+2,axisbg=(.95,.95,.95),xlim=(numpy.min(Hx)*(1-AxisMargin),numpy.max(Hx)*(1+AxisMargin)),ylim=(0,numpy.max(Hmean)*(1+AxisMargin)),xlabel='Radius, m',ylabel='[Rel.] Volume Fraction',xscale='log'))
            else:
                R_ax.append(fig.add_subplot(1,(nhists+1),histi+2,axisbg=(.95,.95,.95),xlim=(numpy.min(Hx)-(1-AxisMargin)*numpy.min(Hx),numpy.max(Hx)*(1+AxisMargin)),ylim=(0,numpy.max(Hmean)*(1+AxisMargin)),xlabel='Radius, m',ylabel='[Rel.] Volume Fraction'))

            R_ax[histi]=setaxis(R_ax[histi])
            #fill axes
            bar(Hx[0:-1],Hmean,width=Hwidth,color='orange',edgecolor='black',linewidth=1,zorder=2,label='MC size histogram')
            plot(Hmid,vfminbins,'ro',ms=5,markeredgecolor='r',label='Minimum visibility limit',zorder=3)
            errorbar(Hmid,Hmean,Hstd,zorder=4,fmt='k.',ecolor='k',elinewidth=2,capsize=4,ms=0,lw=2,solid_capstyle='round',solid_joinstyle='miter')
            legend(loc=1,fancybox=True,prop=textfont)
            title('Radius size histogram',fontproperties=textfont,size='x-large')
            #reapply limits in x
            xlim((numpy.min(Hx)*(1-AxisMargin),numpy.max(Hx)*(1+AxisMargin)))


        fig.subplots_adjust(left=0.1,bottom=0.11,right=0.96,top=0.95,wspace=0.23,hspace=0.13)
        
    def Rangeinfo(self,ParameterRange=[0,inf],Parameter=0):
        """Calculates the total volume or number fraction of the MC result
        within a given range, and returns the total numer or volume fraction
        and its standard deviation over all nreps as well as the first four
        distribution moments: mean, variance, skewness and kurtosis
        (Pearson's definition).
        Will use the *Histweight* parameter for determining whether to return
        the volume or number-weighted values.

        Input arguments are:

            *ParameterRange*
              The radius range in which the moments are to be calculated
            *Parameter*
              Which shape parameter the moments are to be calculated for
              (e.g. 0 = width, 1 = length, 2 = orientation)

        Returns a 4-by-2 array, with the values and their sample standard
        deviations over all *Nreps*.
        """
        Rrep=self.getresult('Rrep')
        Ncontrib=size(Rrep,0)
        NRval=size(Rrep,1)
        Nreps=size(Rrep,2)
        Rpfactor=self.getpar('Rpfactor')
        Memsave=self.getpar('Memsave')
        drhosqr=self.getpar('drhosqr')
        Histbins=self.getpar('Histbins')
        Histscale=self.getpar('Histscale')
        Histweight=self.getpar('Histweight')
        Bounds=self.getpar('Bounds')
        
        #ov = zeros(shape(Rrep)) #observability
        Vf = self.getresult('Vf') #volume fraction for each contribution
        Nf = self.getresult('Nf') #number fraction for each contribution
        Vft = self.getresult('Vft') #total volume fractions
        Nft = self.getresult('Nft') #total number 
        Screps = self.getresult('Screps') #Intensity scaling factors for matching to the experimental scattering pattern (Amplitude A and flat background term b, defined in the paper)

        Val=zeros(Nreps) #total value
        Mu=zeros(Nreps) #moments..
        Var=zeros(Nreps) #moments..
        Skw=zeros(Nreps) #moments..
        Krt=zeros(Nreps) #moments..

        #loop over each repetition
        for ri in range(Nreps):
            Rset = Rrep[:,Parameter,ri] #the single set of R for this calculation
            validi=(Rset>numpy.min(ParameterRange))*(Rset<numpy.max(ParameterRange))
            Rset=Rset[validi]
            Vset = Vf[validi,ri] #compensated volume for each sphere in the set
            Nset = Nf[validi,ri] #compensated volume for each sphere in the set

            if Histweight=='volume':
                Val[ri] = sum(Vset)
                Mu[ri] = sum(Rset*Vset)/sum(Vset)
                Var[ri] = sum( (Rset-Mu[ri])**2*Vset )/sum(Vset)
                sigma=sqrt(abs(Var[ri]))
                Skw[ri] = sum( (Rset-Mu[ri])**3*Vset )/(sum(Vset)*sigma**3)
                Krt[ri] = sum( (Rset-Mu[ri])**4*Vset )/(sum(Vset)*sigma**4)
            elif Histweight=='number':
                Val[ri]=sum(Nset)
                Mu[ri] = sum(Rset*Nset)/sum(Nset)
                Var[ri] = sum( (Rset-Mu[ri])**2*Nset )/sum(Nset)
                sigma=sqrt(abs(Var[ri]))
                Skw[ri] = sum( (Rset-Mu[ri])**3*Nset )/(sum(Nset)*sigma**3)
                Krt[ri] = sum( (Rset-Mu[ri])**4*Nset )/(sum(Nset)*sigma**4)
            else:
                print('Error in moment calculation, unrecognised Histweight value')
                return None

        return numpy.array([[numpy.mean(Val),numpy.std(Val,ddof=1)],
                [numpy.mean(Mu),numpy.std(Mu,ddof=1)],
                [numpy.mean(Var),numpy.std(Var,ddof=1)],
                [numpy.mean(Skw),numpy.std(Skw,ddof=1)],
                [numpy.mean(Krt),numpy.std(Krt,ddof=1)]])
        

########################## END McSAS OBJECT ############################

#some quick pickle functions to make my life easier

def pickle_read(filename):
    """*\*args* can be 1-4, indicates number of output variables.
    If it is even possible to extract more from pickle."""
    fh=open(filename)
    O=pickle.load(fh)
    fh.close()
    return O

def pickle_write(filename,DBlock):
    """Writes DBlock to a file."""
    fh=open(filename,'w')
    pickle.dump(DBlock,fh)
    fh.close()
    return

def binning_array(Q,PSI,I,IERR,S=2):
    """This function applies a simple S-by-S binning routine on images.
    It calculates new error based on old error superseded by standard
    deviation in a bin."""
    def isodd(x):
        #checks if a value is even or odd
        return bool(x&1)
    ddi={'Q':Q,'PSI':PSI,'I':I,'IERR':IERR}
    
    sq=shape(Q)
    if isodd(sq[0]):
        #trim edge
        for it in ddi.keys():
            ddi[it]=ddi[it][1:,:]
    if isodd(sq[1]):
        #trim edge
        for it in ddi.keys():
            ddi[it]=ddi[it][:,1:]
    #now we can do n-by-n binning
    sq=shape(Q)
    qo=zeros((sq[0]/S,sq[1]/S))
    ddo={'Q':qo.copy(),'PSI':qo.copy(),'I':qo.copy(),'IERR':qo.copy()}
    for it in ['Q','PSI','I']:
        for ri in range(sq[0]/S):
            for ci in range(sq[1]/S):
                ddo[it][ri,ci]=numpy.mean(ddi[it][S*ri:(S*ri+S),S*ci:(S*ci+S)])
        
    for ri in range(sq[0]/S):
        for ci in range(sq[1]/S):
            meanE=sqrt(sum((ddi['IERR'][S*ri:(S*ri+S),S*ci:(S*ci+S)])**2))/S**2
            stdI=numpy.std(ddi['I'][S*ri:(S*ri+S),S*ci:(S*ci+S)])#sample standard deviation
            #stdI=0
            ddo['IERR'][ri,ci]=numpy.max((meanE,stdI))
    return ddo

def binning_1D(q,I,E=[],Nbins=200,Stats='STD'):
    """An unweighted binning routine.
    The intensities are sorted across bins of equal size.
    If error provided is empty, the standard deviation of the intensities in
    the bins are computed."""
    # Let's make sure the input is consistent
    if size(q)!=size(I):
        print "Incompatible sizes of q and I"
        return
    elif (size(E)!=0) & (size(E)!=size(I)):
        print "Size of E is not identical to q and I"
        return

    #flatten q, I and E
    q=reshape(q,size(q),0)
    I=reshape(I,size(I),0)
    E=reshape(E,size(E),0)

    #define the bin edges and centres, and find out the stepsize while we're at it. Probably, there is no need for knowing the edges...
    qbin_edges=linspace(numpy.min(q),numpy.max(q),Nbins+1)
    stepsize=qbin_edges[1]-qbin_edges[0]
    qbin_centres=linspace(numpy.min(q)+0.5*stepsize,numpy.max(q)-0.5*stepsize,Nbins)
    
    #sort q, let I and E follow sort
    sort_ind=numpy.argsort(q,axis=None)
    q=q[sort_ind]
    I=I[sort_ind]
    Ibin=numpy.zeros(Nbins)
    Ebin=numpy.zeros(Nbins)    
    SDbin=numpy.zeros(Nbins)    
    SEbin=numpy.zeros(Nbins)    
    if (size(E)!=0):
        E=E[sort_ind]

    #now we can fill the bins
    for Bini in range(Nbins):
        #limit ourselves to only the bits we're interested in:
        lim_bool_array=((q>(qbin_centres[Bini]-stepsize))&(q<=(qbin_centres[Bini]+stepsize)))
        I_to_bin=I[lim_bool_array]
        if (size(E)!=0):
            E_to_bin=sum(E[lim_bool_array])

        #find out the weighting factors for each q,I,E-pair in the array  
        weight_fact=ones(size(I_to_bin))

        #sum the intensities in one bin and normalize by number of pixels
        Ibin[Bini]=sum(I_to_bin)/sum(weight_fact); 

        

        #now we deal with the Errors:
        if (size(E)!=0): #if we have errors supplied from outside
            #standard error calculation:
            SEbin[Bini]=sqrt(sum(E_to_bin**2*weight_fact))/sum(weight_fact)
            #Ebin2[Bini]=sqrt(sum((E_to_bin*weight_fact)**2))/sum(weight_fact) old, incorrect
            if Stats=='auto':
                SDbin[Bini]=sqrt(sum((I_to_bin-Ibin[Bini])**2*weight_fact)/(sum(weight_fact)-1)) #according to the definition of sample-standard deviation
                SEbin[Bini]=numpy.max(array([SEbin[Bini],SDbin[Bini]/sqrt(sum(weight_fact))])) #maximum between standard error and Poisson statistics

        else:           
            #calculate the standard deviation of the intensity in the bin both for samples with supplied error as well as for those where the error is supposed to be calculated
            SDbin[Bini]=sqrt(sum((I_to_bin-Ibin[Bini])**2*weight_fact)/(sum(weight_fact)-1)) #according to the definition of sample-standard deviation
            #calculate standard error by dividing the standard error by the square root of the number of measurements
            SEbin[Bini]=SDbin[Bini]/sqrt(sum(weight_fact))

    return qbin_centres,Ibin,SEbin

def binning_weighted_1D(q,I,E=[],Nbins=200,Stats='SE'):
    """Implementation of the binning routine written in Matlab.
    The intensities are divided across the q-range in bins of equal size.
    The intensities of a pixel are divided between the two neighbouring bins
    depending on the distances to the centres. If error provided is empty,
    the standard deviation of the intensities in the bins are computed.

    **Usage**::

        qbin, Ibin, Ebin = binning_weighted_1D(q, I, E = [],
                                               Nbins = 200, Stats = 'SE')

    **Optional input arguments**:

    *Nbins*: integer indicating the number of bins to divide the intensity
        over. Alternatively, this can be an array of equidistant bin centres.
        If you go this route, depending on the range, not all intensity may be
        counted.
    *Stats*: Can be set to 'auto'. This takes the maximum error between
        supplied Poisson statistics error-based errors or the standard error.

    Written by Brian R. Pauw, 2011, released under BSD open source license.
    """
    #let's make sure the input is consistent
    if size(q)!=size(I):
        print "Incompatible lengths of q and I, q must be of the same number of elements as I"
        return
    elif (size(E)!=0) & (size(E)!=size(I)):
        print "Size of E is not identical to q and I"
        return
    if (Stats.lower()!='std')&(Stats.lower()!='poisson')&(Stats.lower()!='se')&(Stats.lower()!='auto'):
        print "Statistics can only be set to 'SE' (default), or 'auto'. \n"
        print "Only use 'auto' for photon-counting detectors, selects largest error between SE and Poisson.\n"
        print "If errors are supplied, standard errors are not calculated except in the case of 'auto' "
        return
    if (size(Nbins)==1):
        if (Nbins<1):
            print "number of bins, Nbins, is smaller than one. I need at least one bin to fill"
            return
    if size(Nbins)>1:
        print "Nbins is larger than one value. Assuming that an equidistant list of bin centres has been supplied"

    #flatten q, I and E
    q=reshape(q,size(q),0)
    I=reshape(I,size(I),0)
    E=reshape(E,size(E),0)
    
    if size(Nbins)==1:
        #define the bin edges and centres, and find out the stepsize while we're at it. Probably, there is no need for knowing the edges...
        qbin_edges,stepsize=linspace(numpy.min(q),numpy.max(q),Nbins+1,retstep=True)#do not use qbin_edges!
        #stepsize=qbin_edges[1]-qbin_edges[0]
        qbin_centres=linspace(numpy.min(q)+0.5*stepsize,numpy.max(q)-0.5*stepsize,Nbins)
    else:
        if (numpy.min(q)>numpy.max(Nbins))|(numpy.max(q)<numpy.min(Nbins)):
            print "Bin centres supplied do not overlap with the q-range, cannot continue"
            return
        qbin_centres=sort(Nbins)
        stepsize=mean(diff(qbin_centres))
        Nbins=size(qbin_centres)

    #initialize output matrices
    Ibin=numpy.zeros(Nbins)
    SDbin=numpy.zeros(Nbins)    
    SEbin=numpy.zeros(Nbins)    

    #now we can fill the bins
    for Bini in range(Nbins):
        #limit ourselves to only the bits we're interested in:
        lim_bool_array=((q>(qbin_centres[Bini]-stepsize))&(q<=(qbin_centres[Bini]+stepsize)))
        q_to_bin=q[lim_bool_array]
        I_to_bin=I[lim_bool_array]
        if (size(E)!=0):
            E_to_bin=E[lim_bool_array]

        #find out the weighting factors for each q,I,E-pair in the array  
        q_dist=abs(q_to_bin-qbin_centres[Bini])
        weight_fact=(1-q_dist/stepsize)

        #sum the intensities in one bin
        Ibin[Bini]=sum(I_to_bin*weight_fact)/sum(weight_fact); 

        #now we deal with the Errors:
        if (size(E)!=0): #if we have errors supplied from outside
            #standard error calculation:
            SEbin[Bini]=sqrt(sum(E_to_bin**2*weight_fact))/sum(weight_fact)
            #Ebin2[Bini]=sqrt(sum((E_to_bin*weight_fact)**2))/sum(weight_fact) old, incorrect
            if Stats=='auto':
                SDbin[Bini]=sqrt(sum((I_to_bin-Ibin[Bini])**2*weight_fact)/(sum(weight_fact)-1)) #according to the definition of sample-standard deviation
                SEbin[Bini]=numpy.max(array([SEbin[Bini],SDbin[Bini]/sqrt(sum(weight_fact))])) #maximum between standard error and Poisson statistics
        else:           
            #calculate the standard deviation of the intensity in the bin both for samples with supplied error as well as for those where the error is supposed to be calculated
            SDbin[Bini]=sqrt(sum((I_to_bin-Ibin[Bini])**2*weight_fact)/(sum(weight_fact)-1)) #according to the definition of sample-standard deviation
            #calculate standard error by dividing the standard error by the square root of the number of measurements
            SEbin[Bini]=SDbin[Bini]/sqrt(sum(weight_fact))
    return qbin_centres,Ibin,SEbin
        
    
##general functions
#def csqr(Sc,I,Ic,E):
#    """Least-squares error for use with scipy.optimize.leastsq"""
#    cs=(I-Sc[0]*Ic-Sc[1])/E
#    #print "Size E",size(E)
#    #print "Sc: %f, %f" %(Sc[0],Sc[1])
#    return cs
#
#def Iopt(I,Ic,E,Sc,OutputI=False):
#    """Optimizes the scaling and background factor to match Ic closest to I.
#    Returns an array with scaling factors. Input Sc has to be a two-element
#    array with the scaling and background.
#
#    New version, using leastsq and csqr, speed improvement of over
#    a factor of 10 w.r.t. V1's use in the MC algorithm
#    """
#    # initial guesses. No bounds at the moment are applied,
#    # except that the intensity scaling has to be positive.
#    Sc,success=scipy.optimize.leastsq(csqr,Sc,args=(I.flatten(),Ic.flatten(),E.flatten()),full_output=0)
#    #print "Sc: %f, %f" %(Sc[0],Sc[1])
#    cval=csqr_v1(I,Sc[0]*Ic+Sc[1],E)
#    if OutputI:
#        return Sc,cval,Sc[0]*Ic+Sc[1]
#    else:
#        return Sc,cval
#
#def csqr_v1(I,Ic,E):
#    """Least-squares for data with known error,
#    size of parameter-space not taken into account."""
#    cs=sum(((I-Ic)/E)**2)/(size(I))
#    return cs
#
#def Iopt_v1(I,Ic,E,Sc,OutputI=False,Background=True):
#    """Optimizes the scaling and background factor to match Ic closest to I.
#    Returns an array with scaling factors. Input Sc has to be a two-element
#    array with the scaling and background.
#
#    Old version, using fmin and csqr_v1.
#    """
#    # initial guesses. No bounds at the moment are applied,
#    # except that the intensity scaling has to be positive.
#    # Background can be set to False to just find the scaling factor.
#    if Background:
#        Sc=scipy.optimize.fmin(lambda Sc: csqr_v1(I,Sc[0]*Ic+Sc[1],E),Sc,full_output=0, disp=0)
#        cval=csqr_v1(I,Sc[0]*Ic+Sc[1],E)
#    else:
#        Sc[1]=0.
#        Sc=scipy.optimize.fmin(lambda Sc: csqr_v1(I,Sc[0]*Ic,E),Sc,full_output=0, disp=0)
#        Sc[1]=0.
#        cval=csqr_v1(I,Sc[0]*Ic,E)
#    if OutputI:
#        return Sc,cval,Sc[0]*Ic+Sc[1]
#    else:
#        return Sc,cval


# vim: set ts=4 sts=4 sw=4 tw=0:
