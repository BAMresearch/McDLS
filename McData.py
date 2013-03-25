# -*- coding: utf-8 -*-
# McData.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

r'''
A class and supplementary functions for Monte-Carlo fitting of SAXS patterns.
                    *Data and parameter read-in procedures*

It is released under a `Creative Commons CC-BY-SA license
<http://creativecommons.org/licenses/by-sa/3.0/>`_.
Please cite as::

    Brian R. Pauw et al., J. Appl. Cryst. 46, (2013), pp. 365--371
        doi: http://dx.doi.org/10.1107/S0021889813001295

Classes and Functions Defined in This File
==========================================

 - **McData**: class for constructing the required set of information for use
   with the rest of the McSAS suite for Monte Carlo analysis of 
   small-angle scattering data.
 - **binning_array**: Can be used to do n-by-n pixel binning of 2D detector
   images. The returned uncertainty is the larger of either the binned
   uncertainty or the sample standard deviation in the bin.
 - **binning_1D**: bins the data and propagates errors, or calculates errors
   if not initially provided
 - **binning_weighted_1D**: Weighted binning, where the intensities of a
   pixel are divided between the two neighbouring bins depending on the
   distances to the centres. If error provided is empty, the standard
   deviation of the intensities in the bins are computed.
 - **pickle_read**: Reads in pickled data from a file (by filename)
 - **pickle_write**: write a block or dictionary to a file (by filename)
    
Made possible with help from (amongst others)
=============================================

 - | Samuel Tardif
   | Derivations (mostly observability) and checking of mathematics
 - | Jan Skov Pedersen
   | checking of mathematics
 - | Pawel Kwasniewski <kwasniew@esrf.fr>
   | Code cleanup and documentation
 - | Ingo Bressler <ingo.bressler@bam.de>
   | Code cleanup, modification and documentation

A Note on Units
===============

Internally, all length units are in meters, all angle units in degrees
clockwise from top. *Intensity* is in
:math:`\left[ 1 \over {m \cdot sr} \right]`,
*q* in :math:`\left[ 1 \over m \right]`.
The electron density contrast squared,
*drhosqr* is in :math:`\left[ m^{-4} \right]`.
Other units may be used, but if absolute units are supplied and absolute
volume fractions required, meters are necessitated.

Example Usage
=============


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

'''

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

class McData(object):
    r"""
    Main class containing all functions required to prepare the settings and 
    data for the Monte Carlo fit procedure.

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
            | Can be set to 'number' to force plotting of number-weighted
              distributions
            | 2013-03-19 issue remains that the observability in this case is
              incorrect.
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

    A McData object with the settings stored in its functions. Data should be 
    set and gotten using the set* and get* functions in order to do input 
    checking.

    """        

    dataset=None #where Q, PSI, I and IERR is stored, original dataset
    fitdata=None #may be populated with a subset of the aforementioned dataset, limited to q-limits or psi limits and to positive I values alone
    parameters=None #where the fitting and binning settings are stored
    functions=None #where the used functions are defined, this is where shape changes, smearing, and various forms of integrations are placed.

    def __init__(self,**kwargs):
        """
        intialization function, takes keyword-value input parameters. 
        input arguments can be one of the aforementioned parameter keyword-value pairs.
        This does the following::
            1. Initialises the variables to the right type
            2. Parses the input
            3. Stores the supplied data twice, one original and one for fitting 
                (the latter of which may be clipped or otherwise adjusted)
            4. Applies Q- and optional PSI- limits to the data
            5. Reshapes fitdata to 1-by-n dimensions
            6. Sets the function references
            7. Calculates the shape parameter bounds if not supplied
            8. Peforms simple checks on validity of input

        """
        #initialize
        self.dataset=dict() #where Q, PSI, I and IERR is stored, original dataset
        self.fitdata=dict() #may be populated with a subset of the aforementioned dataset, limited to q-limits or psi limits and to positive I values alone
        self.parameters=dict() #where the fitting and binning settings are stored
        #self.analysisresult=dict() #The place where the raw analysis details are stored, along with qfit and Ifit.
        #self.interpretresult=dict() #The place for the histogrammed data.
        self.functions=dict() #where the used functions are defined, this is where shape changes, smearing, and various forms of integrations are placed.
        #self.handles=dict() #I don't think handles can be stored in a dict.

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
        self.setfunctions(**kwargs)
        #calculate parameter bounds and store
        self.functions['BOUNDS']()
        #check and fix parameters where necessary
        self.check_parameters() #this is only a very simple function now in need of expansion

        return self

    def setfunctions(self,**kwargs):
        '''functions are defined here. 
        In particular here the following is specified:
        1. The parameter bounds estimation function 'BOUNDS'. Should be able to take input argument Bounds to update, should set the parameter bounds in self.parameter['Bounds']
        2. The random number generation function 'RAND' This must take its parameters from self, and have an optional input argument specifying the number of sets to return (for MC initialization). It should return a set of Nsets-by-nvalues to be used directly in 'FF'. This may be depreciated soon as is can be generated from within.
        3. The Form-factor function 'FF'. If called, this should get the required information from self and a supplied Nsets-by-nvalues shape-specifying parameter array. It should return an Nsets-by-q array. Orientational averaging should be part of the form-factor function (as it is most efficiently calculated there), so several form factor functions can exist for non-spherical objects.
        4. The shape volume calculation function 'VOL', which must be able to deal with input argument "Rpfactor", ranging from 0 to 1. Should accept an Nsets-by-nvalues array returning an Nsets number of (Rpfactor-compensated)-volumes. 
        5. The smearing function 'SMEAR'. Should take information from self, and an input Icalc, to output an Ismear of the same length.

        This function will actually cast the supplied function name into a function pointer.
        '''
        for kw in kwargs:
            if kw in self.functions.keys():
                if callable(kwargs[kw]):
                    self.functions[kw]=kwargs[kw]
                else:
                    #make it into a function handle/pointer
                    self.functions[kw]=getattr(self,kwargs[kw])

    def getfunctions(self,fname=None)
        '''
        returns the function handle or all handles (if no function name supplied).
        function name can be one of the following
        '''
        if fname==None:
            return self.functions
        else:
            return self.functions[fname]

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


    def _passthrough(self,In):
        """a passthrough mechanism returning the input unchanged"""
        return In

    def reshape_fitdata(self):
        """This ensures that q, I, PSI and E are in 1-by-n shape"""
        for key in self.fitdata.keys():
            self.fitdata[key]=reshape(self.fitdata[key],(1,prod(shape(self.fitdata[key]))))

    def clip_dataset(self):
        '''if q and/or psi limits are supplied in self.parameters, clips the dataset to within the supplied limits. Copies data to self.fitdata if no limits are set.'''

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
        '''gets the value of a parameter, so simple it is probably superfluous'''
        if parname==[]:
            return self.parameters
        else:
            return self.parameters[parname]

    def getdata(self,parname=[],dataset='fit'):
        '''gets the values of a dataset, retrieves from fitdata (clipped) by default. If the original data is wanted, use "dataset='original'" as kwarg. '''
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

    def setpar(self,**kwargs):
        '''
        Sets the supplied parameters given in keyword-value pairs for known setting keywords (unknown key-value pairs are skipped)
        If a supplied parameter is one of the function names, it is stored in the self.functions dict.
        '''
        for kw in kwargs:
            if kw in self.parameters.keys():
                self.parameters[kw]=kwargs[kw]
            else:
                pass #no need to store unknown keywords.

    def setdata(self,**kwargs):
        '''
        Sets the supplied data in the proper location. Optional argument "dataset" can be set to "fit" or "original" to define which dataset is set. Default is "original"
        '''
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

    def check_parameters(self):
        '''in here we can place checks for the parameters, for example to make sure histbins is defined for all, or to check if all parameters fall within their limits'''
        #for now, all I need is a check that Histbins is a 1D vector with n values, where n is the number of parameters specifying a shape. 

        testR=self.functions['RAND']()
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

        
        ############## Bounds calculation functions ##############    
    def EllBounds_2D(self):
        '''
        This function will take the q and psi input bounds and outputs properly formatted two-element size bounds for ellipsoids. Ellipsoids are defined by their equatorial radius, meridional radius and axis misalignment (default -45 to 45 degrees in PSI).
        '''
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
        '''
        This function will take the q and input bounds and outputs properly formatted two-element size bounds.
        '''
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

        ############## Random generator functions ##############    
    def random_uniform_ell(self,Nell=1):
        """
        Random number generator for generating uniformly-sampled size- and orientation parameters for ellipsoids.
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
        "Random number generator which behaves like its uniform counterpart, but with a higher likelihood of sampling smaller sizes. May speed up some fitting procedures."
        #get parameters from self
        Bounds=self.getpar('Bounds') 
        #generate Nsph random numbers
        Rset=zeros((Nell,3))
        Rset[:,0]=10**(numpy.random.uniform(log10(numpy.min(Bounds[0])),log10(numpy.max(Bounds[1])),Nell))
        Rset[:,1]=10**(numpy.random.uniform(log10(numpy.min(Bounds[2])),log10(numpy.max(Bounds[3])),Nell))
        Rset[:,2]=numpy.random.uniform(numpy.min(Bounds[4]),numpy.max(Bounds[5]),Nell)

        #output Nsph-by-3 array
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

    ############## Volume calculation functions ##############    
    def vol_ell(self,Rset,Rpfactor=[]):
        '''
        calculates the volume of an ellipsoid, taking Rpfactor from input or preset parameters
        '''
        if Rpfactor==[]:
            Rpfactor=self.getpar('Rpfactor')
            
        return ((4.0/3*pi)*Rset[:,0]**(2*Rpfactor)*Rset[:,1]**(Rpfactor))[:,newaxis]

    def vol_sph(self,Rset,Rpfactor=[]):
        '''
        calculates the volume of a sphere, taking Rpfactor from input or preset parameters
        '''
        if Rpfactor==[]:
            Rpfactor=self.getpar('Rpfactor')
            
        return (4.0/3*pi)*Rset**(3*Rpfactor)

    ############## Form factor calculation functions ##############    
    def FF_sph_1D(self,Rset):
        '''
        Calculate the Rayleigh function for a sphere
        '''
        q=self.getdata('Q')
        if size(Rset,0)>1: # multimensional matrices required, input Rsph has to be Nsph-by-1. q has to be 1-by-N
            qR=(q+0*Rset)*(Rset+0*q)
        else:
            qR=(q)*(Rset)

        Fsph=3*(sin(qR)-qR*cos(qR))/(qR**3)
        return Fsph

    def FF_ell_2D(self,Rset=[],Q=[],PSI=[]):
        """
        Calculates the form factor for oriented ellipsoids, normalized to 1 at Q=0.

        **all 2D functions should be able to potentially take externally supplied Q and PSI vectors**
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


####################################### End Object #########################################
#some quick pickle functions to make my life easier

def pickle_read(filename):
    #nargs can be 1-4, indicates number of output variables, if it is even possible to extract more from pickle
    fh=open(filename)
    O=pickle.load(fh)
    fh.close()
    return O

def pickle_write(filename,DBlock):
    #writes DBlock to a file
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
    #python implementation of an unweighted binning routine. The intensities are sorted across bins of equal size. If error provided is empty, the standard deviation of the intensities in the bins are computed.
    #let's make sure the input is consistent
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
    #USAGE: qbin,Ibin,Ebin=binning_weighted_1D(q,I,E=[],Nbins=200,Stats='SE'):
    #python implementation of the binning routine written in Matlab. The intensities are divided across the q-range in bins of equal size. The intensities of a pixel are divided between the two neighbouring bins depending on the distances to the centres. If error provided is empty, the standard deviation of the intensities in the bins are computed.
    #optional input arguments:
    #   Nbins: integer indicating the number of bins to divide the intensity over. Alternatively, this can be an array of equidistant bin centres. If you go this route, depending on the range, not all intensity may be counted.
    #   Stats: can be set to 'auto'. This takes the maximum error between supplied Poisson statistics error-based errors or the standard error. 
    #Written by Brian R. Pauw, 2011, released under BSD open source license.
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
        
# vim: set ts=4 sts=4 sw=4 tw=0:
