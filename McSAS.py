'''
This is a file with small programs for Monte-Carlo fitting of SAXS patterns. It is 
released under a Creative Commons CC-BY-SA license. Please cite as:

Brian R. Pauw, 2012, http://arxiv.org/abs/1210.5304 arXiv:1210.5304. 

(May be replaced later with a J. Appl. Cryst. reference)

Contents (updated 2013-01-16):
    Analyze_1D: A wrapper for the MC code which repeatedly runs the optimization and 
        computes the final volume-weighted sphere distribution.
    FF_sph_1D: computes the rayleigh form factor for a sphere radius or array of sphere 
        radii
    observability3: histogramming and uncertainty calculation for the Analyze_1D code.
        also calculates the minimum number of required contributions and the volume 
        fractions in absolute units.
    MCFit_sph: Monte-carlo fitting of the data to extract polydispersity assuming spheres, 
        by changing the radius of a fixed number of spheres
    binning_1D: bins the data and propagates errors, or calculates errors if not initially 
        provided
    binning_weighted_1D: Weighted binning, where the intensities of a pixel are divided 
        between the two neighbouring bins depending on the distances to the centres. If 
        error provided is empty, the standard deviation of the intensities in the bins are 
        computed. 
    MCFit_sph: Monte-carlo fitting of the data to extract polydispersity assuming spheres, 
        by changing the radius of a fixed number of spheres
    McPlot: A procedure for generating a data-fit plot and size histogram based on 
        Analyze_1D's results
    McCSV: Function to write an arbitrary number of semicolon-separated values to a file
    FixBounds: Internal function for estimating minimum and maximum size bounds based on 
        q values.
    csqr: least-squares error to use with scipy.optimize.leastsq
    Iopt: Optimize the scaling factor and background level of modeled data vs. intensity
    csqr_v1: least-squares for data with known error, size of parameter-space not taken 
        into account
    Iopt_v1: old intensity scaling factor optimisation, more robust but slower than Iopt
    pickle_read: Reads in pickled data from a file (by filename)
    pickle_write: write a block or dictionary to a file (by filename)
    
    (asterisk * indicates code to be depreciated)

Made possible with help from (amongst others):
    Samuel Tardif - Derivations (mostly observability) and checking of mathematics
    Jan Skov Pedersen - checking of mathematics
    Pawel Kwasniewski <kwasniew@esrf.fr> - Code cleanup and documentation
'''

import scipy # For many important functions
from scipy import optimize # For the leastsq optimization function
import numpy # For arrays
import os # Miscellaneous operating system interfaces
import time # Timekeeping and timing of objects
import sys # For printing of slightly more advanced messages to stdout
from numpy import *
import pickle #for pickle_read and pickle_write

# 1D functions:

###########################################################################################
###################################### Analyze_1D #########################################
###########################################################################################

def Analyze_1D(q,I,E,Bounds=[],Nsph=200,Maxiter=1e5,Rpfactor=1.5/3,Nreps=100,qlims=[0,inf],Histbins=50,Histscale='log',drhosqr=1,Convcrit=1.,StartFromMin=False,Maxntry=5,SimpleOutput=False,Plot=False):
    '''
    This function runs the monte-carlo fit MCFit_sph() several times, and returns 
    bin centres, means, standard deviations and observability limits for the bins. 
    Eventually, this may also plot. The drhosqr value can be set to the contrast 
    (if known, approximately on the order of 6x10^30 m^-2) so a volume fraction 
    for each histogram bin is calculated.
    Be aware that all length units are in meters, including those for q, I and 
    the scattering contrast.

    Parameters:
    -----------
    q : 1D array 
        q vector values, units: m^-1
    I : 1D array 
        I(q) values, units: (m sr)^-1, size should be identical to q
    E : 1D array
        Uncertainty (one standard deviation) of the measured I(q) values, these 
        should only reflect the relative uncertainties, absolute uncertainty 
        contributions should be indicated on the final histogram in addition to 
        the relative uncertainty of the result. Size should be identical to q, 
        values should be positive and larger than zero
    Bounds : list
        Two-element vector or list indicating upper and lower size bounds of the 
        particle radii used in the fitting procedure. If not provided, these will
        be estimated as: Rmax=pi/qmin and Rmin=pi/qmax. units: m
    Nsph : int, default: 200
        Number of spheres used for the MC simulation
    Maxiter : int, default: 1e5
        Maximum number of iterations for the MCfit_sph() function
    Rpfactor : float, default: 1.5/3
        Parameter used to compensate the volume^2 scaling of each sphere 
        contribution to the simulated I(q)
    Nreps : int, default: 100
        Number of repetitions of the MC fit for determination of final histogram 
        uncertainty.
    qlims : list, default: [0,inf]
        Limits on the fitting range in q. units: m^-1
    Histbins : int, default: 50
        Number of bins used for the histogramming procedure.
    Histscale : string, default: 'log'
        Can be set to 'log' for histogramming on a logarithmic size scale, 
        recommended for q- and/or size-ranges spanning more than a decade.
    drhosqr : float, default: 1
        Scattering contrast - when known it will be used to calculate the absolute
        volume fraction of each contribution, units: m^-2
    Convcrit : float, default: 1
        Convergence criterium for the least-squares fit. The fit converges once 
        the normalized chi squared < Convcrit. If convergence is reached with 
        Convcrit = 1, the model describes the data (on average) to within the 
        uncertainty, and thus all information has been extracted from the 
        scattering pattern.
    StartFromMin : bool, default: False
        If set to False, the starting configuration is a set of spheres with radii
        uniformly sampled between the given or estimated bounds.
        If set to True, the starting configuration is a set of spheres with radii
        set to the lower given or estimated Bound (if not zero). Practically, this
        makes little difference and this feature might be depreciated.
    Maxntry : int, default: 5
        If a single MC optimization fails to reach convergence within Maxiter, it
        may just be due to bad luck. The Analyze_1D procedure will try to redo 
        that MC optimization for a maximum of Maxntry tries before concluding that
        it is not bad luck but bad input.
    Plot : Bool, default: False
        If set to True, will generate a plot showing the data and fit, as well as
        the resulting size histogram.

    Returns:
    --------
    A : Output dictionary containing the following elements:
        'Imean' : 1D array 
            The fitted intensity, given as the mean of all Nreps results.
        'q' : 1D array
            Corresponding q values (may be different than the input q if qlims was 
            used) 
        'Istd' : array
            Standard deviation of the fitted I(q), calculated as the standard 
            deviation of all Nreps results.
        'Hx' : array
            Histogram bin left edge position
        'Hmid' : array
            Center positions for the size histogram bins
        'Hwidth' : array
            Histogram bin width
        'Hmean' : array
            Volume-weighted particle size distribution values for all Nreps results
        'Hstd' : array
            Standard deviations of the corresponding size distribution bins, calculated
            from Nreps repetitions of the MCfit_sph() function
        'Hy' : size (Histbins x Nreps) array
            Volume-weighted particle size distribution values for each MC fit repetition
        'Niter' : int
            Average number of MC iterations required for convergence
        'Rrep' : size (Nsph x Nreps) array
            Collection of Nsph sphere radii fitted to best represent the provided I(q) data.
            Contains the results of each of Nreps iterations. This can be used for
            rebinning without having to re-optimize.
        'Vf' : size (Nsph x Nreps) array
            Volume fractions for each of Nsph spheres 
            in each of Nreps iterations
        'Vft' : size (Nreps) array
            Total scatterer volume fraction for each of the Nreps iterations 
        'vfmin' : size (Nsph x Nreps) array
            minimum required volube fraction for each contribution to become statistically
            significant.
        'vfminbins' : size (Hmid) array 
            array with the minimum required volume fraction per bin to become statistically 
            significant. Used to display minimum required level in histogram.
        'Screps' : size (2 x Nreps) array
            Scaling and background values for each repetition. Used to display background 
            level in data and fit plot.


    See also:
    ---------
    MCfit_sph, observability3

    Usage:
    ------
    To fit I(q):
    A=Analyze_1D(q,I,numpy.maximum(0.01*I,E),Nsph=200,Convcrit=1,Bounds=array([pi/numpy.max(q),pi/numpy.min(q)]),Rpfactor=1.5/3,Maxiter=1e5,Histscale='lin',drhosqr=1,Nreps=100)
    Or, simplified:
    A=Analyze_1D(q,I,numpy.maximum(0.01*I,E))

    Plotting the data fit and histogram:
    McPlot(q,I,numpy.maximum(0.01*I,E),A)

    '''
    #for volume weighting, Rpfactor = 1.5/3
    #store initial values for q, I and E for plotting
    initialq=q
    initialI=I
    initialE=E
    # q and psi limits
    validbools = (q>qlims[0])&(q<=qlims[1]) #excluding the lower q limit may prevent q=0 from appearing
    I = I[validbools]
    if size(E)!=0:
        E = E[validbools]
    q = q[validbools]
    Rrep = zeros([Nsph,Nreps]) #Nsph needs to be set for this, hence the default.
    Niters = zeros([Nreps])
    Irep = zeros([len(I),Nreps])
    bignow = time.time() #for time estimation and reporting
    #fix bounds
    Bounds=FixBounds(q,Bounds=Bounds)

    #This is the loop that repeats the MC optimization Nreps times, after which we can calculate an uncertainty on the results.
    for nr in arange(0,Nreps):
        nt = 0 #keep track of how many failed attempts there have been 
        # do that MC thing! 
        Rrep[:,nr],Irep[:,nr],ConVal,Details = MCFit_sph(q,I,E,Bounds=Bounds,Nsph=Nsph,Maxiter=Maxiter,Rpfactor=Rpfactor,OutputI=True,Convcrit=Convcrit,StartFromMin=StartFromMin,OutputDetails=True)
        while ConVal>Convcrit:
            #retry in the case we were unlucky in reaching convergence within Maxiter.
            nt+=1
            Rrep[:,nr],Irep[:,nr],ConVal,Details=MCFit_sph(q,I,E,Bounds=Bounds,Nsph=Nsph,Maxiter=Maxiter,Rpfactor=Rpfactor,OutputI=True,Convcrit=Convcrit,StartFromMin=StartFromMin,OutputDetails=True)
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
    Imean = numpy.mean(Irep,axis=1) #mean fitted intensity
    Istd = numpy.std(Irep,axis=1) #standard deviation on the fitted intensity, usually not plotted for clarity
    # store in output dict
    A = dict()
    A['Rrep'] = Rrep
    A['Imean'] = Imean
    A['Istd'] = Istd
    A['q'] = q
    A['Niter'] = numpy.mean(Niters) #average number of iterations for all repetitions

    #the observability3 function histograms the results and can be used independently of the MC code for rebinning the result. 
    if SimpleOutput:
        return A
    
    print "histogramming..."
    B = observability3(q,I=I,E=E,Rrep=Rrep,Histbins=Histbins,Histscale=Histscale,Bounds=Bounds,Rpfactor=Rpfactor,drhosqr=drhosqr)
    #copy all content of the result of observability3 to the output matrix
    for keyname in B.keys():
        A[keyname] = B[keyname]

    if Plot:
        McPlot(initialq,initialI,initialE,A,Histscale=Histscale)

    return A

######################################## end ###############################################
def FixBounds(q,Bounds=[]):
    '''
    This function will take the q and input bounds and outputs properly formatted two-element size bounds.
    '''
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
    return Bounds

def FF_sph_1D(q,Rsph):
    #this function calculates the rayleigh function for a sphere
    if size(Rsph)>1: #multimensional matrices required
        Rsph=Rsph[:,newaxis] #change the dimension of Rsph array        
        qR=(q+0*Rsph)*(Rsph+0*q)
    else:
        qR=(q)*(Rsph)

    Fsph=3*(sin(qR)-qR*cos(qR))/(qR**3)
    return Fsph

###########################################################################################
#################################### observability3 #######################################
###########################################################################################
def observability3(q,I=[],E=[],Rrep=[],Histbins=30,Histscale='lin',Bounds=[],Rpfactor=1.,drhosqr=1.):
    '''
    Observability calculation for a series of spheres, over a range of q. 
    Additional intensity and errors may be supplied for error-weighted observability. 
    Intensity is used for determining the intesity scaling and background levels.
    
    Now with rebinning as well, so we can keep track of which contribution ends up in 
    which bin and calculate the correct minimum required contribution accordingly.
    '''
    #fix bounds in case not supplied
    Bounds=FixBounds(q,Bounds)
    #set the bin edges for our radius bins either based on a linear division or on a logarithmic division of radii.
    if Histscale == 'lin':
        Hx = linspace(Bounds[0],Bounds[1],Histbins+1) #Hx contains the Histbins+1 bin edges, or class limits.
    else:
        Hx = 10**(linspace(log10(Bounds[0]),log10(Bounds[1]),Histbins+1))

    nreps = size(Rrep,1) #this binning and observability calculation process can be applied to a number nreps of identical but independent MC optimizations, with the nreps repetitions allowign for the estimation of the uncertainty on the final result.
    #ov = zeros(shape(Rrep)) #observability
    Vf = zeros(shape(Rrep)) #volume fraction for each contribution
    qm = zeros(shape(Rrep)) #q-value at which the observability contribution for this component is largest
    vfmin = zeros(shape(Rrep)) #number of required spheres of that size to make a measurable impact
    vfminbin = zeros([Histbins,nreps]) #minimum required number of contributions /in a bin/ to make a measurable impact
    Vft = zeros([nreps]) #total volume fractions
    Hy = zeros([Histbins,nreps]) #total volume fraction contribution in a bin
    Screps = zeros([2,nreps]) #Intensity scaling factors for matching to the experimental scattering pattern (Amplitude A and flat background term b, defined in the paper)

    #loop over each repetition
    for ri in range(nreps):
        Rset = Rrep[:,ri] #the single set of R for this calculation
        FFset = FF_sph_1D(q,Rset) #Form factors, all normalized to 1 at q=0.
        Vset = (4./3*pi)*Rset**(3*Rpfactor) #compensated volume for each sphere in the set
        # Calculate the intensities
        Iset = FFset**2*(Vset[:,newaxis]+0*FFset)**2 # Intensity for each contribution as used in the MC calculation
        Vst = sum(Vset**2) # total compensated volume squared 
        It = sum(Iset,0) # the total intensity of the scattering pattern
        if I==[]:
            I = It #in this case we are not comparing to a measured intensity (untested)
        if E==[]:
            E = 0.01*I #in this case we assume 1% standard deviation uncertainty on the data
        
        # Now for each sphere, calculate its volume fraction (p_c compensated):
        Vsa = 4./3*pi*Rset**(3*Rpfactor)
        # And the real particle volume:
        Vpa = 4./3*pi*Rset**(3)

        Sci = numpy.max(I)/numpy.max(It) #initial guess for the scaling factor.
        Sc,Cv = Iopt(I,It,E,[Sci,1]) #optimize scaling and background for this repetition
        Screps[:,ri]=Sc #scaling and background for this repetition.
        Vf[:,ri] = Sc[0]*Vsa**2/(Vpa*drhosqr) # a set of volume fractions
        Vft[ri] = sum(Vf[:,ri]) # total volume squared
        for isi in range(size(Iset,0)): #For each sphere
            #ov[isi,ri] = (Iset[isi,:]/(It)).max() #calculate the observability (the maximum contribution for that sphere to the total scattering pattern) NOTE: no need to compensate for p_c here, we work with volume fraction later which is compensated by default. additionally, we actually do not use this value.
            qmi = numpy.argmax(Iset[isi,:]/(It)) #determine where this maximum observability is of contribution isi (index)
            qm[isi,ri] = q[qmi] #point where the contribution of isi is maximum
            #close approximation:
            #vfmin[isi,ri] = (E[qmi]*Vf[isi,ri]/(Sc[0]*Iset[isi,qmi]))
            #or more precice but slower:
            vfmin[isi,ri] = numpy.min(E*Vf[isi,ri]/(Sc[0]*Iset[isi,:]))

        # Now bin whilst keeping track of which contribution ends up in which bin:
        for bini in range(Histbins):
            findi = ((Rset>=Hx[bini])*(Rset<Hx[bini+1])) #indexing which contributions fall into the radius bin
            Hy[bini,ri] = sum(Vf[findi,ri]) #y contains the volume fraction for that radius bin
            if sum(findi)==0:
                vfminbin[bini,ri] = 0
            else:
                vfminbin[bini,ri] = numpy.max(vfmin[findi,ri])
                vfminbin[bini,ri] = numpy.mean(vfmin[findi,ri])
            if isnan(Hy[bini,ri]):
                Hy[bini,ri] = 0.
    Hmid = zeros(Histbins)
    vfminbins = zeros(Histbins)
    for bini in range(Histbins):
        Hmid[bini] = numpy.mean(Hx[bini:bini+2])
        vb = vfminbin[bini,:]
        vfminbins[bini] = numpy.max(vb[vb<inf])
    Hmean = numpy.mean(Hy,axis=1)
    Hstd = numpy.std(Hy,axis=1)
    B = dict()
    #B['ov'] = ov
    #B['qm'] = qm
    B['Hx'] = Hx
    B['Hy'] = Hy
    B['Hmid'] = Hmid
    B['Hmean'] = Hmean
    B['Hstd'] = Hstd
    B['Hwidth'] = diff(Hx)
    B['vfminbins'] = vfminbins
    B['vfmin'] = vfmin
    B['Vf'] = Vf
    B['Vft'] = Vft
    B['Screps'] = Screps
    return B
######################################## end ###############################################

###########################################################################################
################################### Monte-carlo procedure #################################
###########################################################################################
def MCFit_sph(q,I,E,Nsph=200,Bounds=[],Convcrit=1.,Rpfactor=1.5/3,Maxiter=1e5,Prior=[],Qlimits=numpy.array([]),MaskNegI=False,OutputI=False,StartFromMin=False,OutputDetails=False,OutputIterations=False):
    '''
    Rewrite of the monte-carlo method previously implemented in Matlab. 
    Simpler form, but open source might mean slight improvements.
    '''
    # Initialise parameters
    #OutputI can be set to True only if MaskNeg is not True
    #Convcrit=1 #reasonable value for poisson weighting. any lower than this and we would be fitting noise
    Weighting = 'Poisson' # only one implemented
    Method = 'Randmove' # only one implemented
    if MaskNegI == True:
        OutputI = False
    if (size(Qlimits) == 0): # no q-limits supplied
        #Qlimits = numpy.min(q)
        Qlimits = 0.
    if (size(Qlimits) == 1): # only lower q limit supplied
        #I do not think people supply integers as limits. This code should be removed for clarity:
        #if isinstance(Qlimits,int): # integer supplied, removing a number of values
        #    Qlimits = q[Qlimits]
        Qlimits = append(Qlimits,numpy.max(q))
    if size(Qlimits) == 2: # make sure they are in the right order
        #I do not think people supply integers as limits. This code should be removed for clarity:
        #if isinstance(Qlimits[0],int): # integer supplied, removing a number of values
        #    Qlimits[0] = q[Qlimits[0]]
        #if isinstance(Qlimits[1],int): # integer supplied, removing a number of values
        #    Qlimits[1] = q[Qlimits[1]]
        Qlimits = numpy.array([numpy.min(Qlimits),numpy.max(Qlimits)])
    # Apply limits
    I = I[(q > Qlimits[0])&(q <= Qlimits[1])] #this should exclude q=0 (you never know what users do)
    E = E[(q > Qlimits[0])&(q <= Qlimits[1])]
    q = q[(q > Qlimits[0])&(q <= Qlimits[1])]
    # Optional masking of negative intensity
    if MaskNegI == True:
        q = q[I >= 0]
        E = E[I >= 0]
        I = I[I >= 0]
    #fix bounds:
    Bounds=FixBounds(q,Bounds)

    #Rpower=2 #squared results in most reasonable end result. none=3, linear=2.5, squared=2, cubed=1.5, fourth=1. This reduces the volume dependency in the calculation leading to quicker convergence
    #Maxiter=1e6
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
            if numpy.min(Bounds)==0:
                #mb=1./Nsph*pi/numpy.max(q)
                mb=pi/numpy.max(q)
            else:
                #mb=numpy.min(Bounds)/numpy.float(Nsph)
                mb=numpy.min(Bounds)
            Rset=numpy.ones(Nsph)*mb
        else:
            Rset=numpy.random.uniform(numpy.min(Bounds),numpy.max(Bounds),Nsph)
    elif (size(Prior)!=0)&(size(Nsph)==0):
        Nsph=size(Prior)
        Rset=Prior
    elif size(Prior)==Nsph:
        Rset=Prior
    elif size(Prior)<Nsph:
        print "size of prior is smaller than Nsph. duplicating random prior values"
        #while size(Prior)<Nsph:
        Addi=numpy.random.randint(size(Prior),size=Nsph-size(Prior))
        Rset=concatenate((Prior,Prior[Addi]))
        print "size now:", size(Rset)
    elif size(Prior)>Nsph:
        print "Size of prior is larger than Nsph. removing random prior values"
        Remi=numpy.random.randint(size(Prior),size=Nsph) #remaining choices
        Rset=Prior[Remi]
        print "size now:", size(Rset)
    
    Arange=array(range(Nsph)) #indices to array values, we'll need this later for some logical indexing operations

    #calculate their form factors
    FFset=FF_sph_1D(q,Rset)
    Vset=(4.0/3*pi)*Rset**(3*Rpfactor)
    #calculate the intensities
    Iset=FFset**2*(Vset[:,newaxis]+0*FFset)**2 #a set of intensities
    Vst=sum(Vset**2) # total volume squared
    It=sum(Iset,0) # the total intensity - eq. (1)
    # Optimize the intensities and calculate convergence criterium
    Sc,Conval1=Iopt_v1(I,It/Vst,E,[1,1]) # V1 is more robust w.r.t. a poor initial guess
    Sc,Conval=Iopt(I,It/Vst,E,Sc) # reoptimize with V2, there might be a slight discrepancy in the residual definitions of V1 and V2 which would prevent optimization.
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
        Rt=numpy.random.uniform(numpy.min(Bounds),numpy.max(Bounds),1)
        Ft=FF_sph_1D(q,Rt)
        Vtt=(4.0/3*pi)*Rt**(3*Rpfactor)
        Itt=(Ft**2*Vtt**2)
        # Calculate new total intensity
        Itest=(It-Iset[Ri,:]+Itt) # we do subtractions and additions, which give us another factor 2 improvement in speed over summation and is much more scalable
        Vstest = (sqrt(Vst)-Vset[Ri])**2+Vtt**2
        # optimize intensity and calculate convergence criterium
        Sct,Convalt = Iopt(I,Itest/Vstest,E,Sc) # using version two here for a >10 times speed improvement
        # test if the radius change is an improvement:
        if Convalt<Conval: # it's better
            Rset[Ri],Iset[Ri,:],It,Vset[Ri],Vst,Sc,Conval=(Rt,Itt,Itest,Vtt,Vstest,Sct,Convalt)
            print "Improvement in iteration number %i, Chi-squared value %f of %f\r" %(Niter,Conval,Convcrit),
            Nmoves+=1
            if OutputIterations:
                # output each iteration, starting with number 0. 
                # Iterations will be stored in Details['Itersph'], Details['IterIfit'], 
                # Details['IterConval'], Details['IterSc'] and 
                # Details['IterPriorUnaccepted'] listing the unaccepted 
                # number of moves before the recorded accepted move.
                Details['Itersph']=concatenate((Details['Itersph'],Rset[:,newaxis]),axis=1) #new iterations will (have to) be appended to this, cannot be zero-padded due to size constraints
                Details['IterIfit']=concatenate((Details['IterIfit'],(Itest/Vstest*Sct[0]+Sct[1])[:,newaxis]),axis=1) #ibid.
                Details['IterConVal']=concatenate((Details['IterConVal'],numpy.array(Convalt)[newaxis]))
                Details['IterSc']=concatenate((Details['IterSc'],Sct[:,newaxis]),axis=1)
                Details['IterPriorUnaccepted']=concatenate((Details['IterPriorUnaccepted'],numpy.array(Nnotaccepted)[newaxis]))
            Nnotaccepted=-1
        # else nothing to do
        Ri+=1 # move to next sphere in list
        Ri=Ri%(Nsph) # loop if last sphere
        Nnotaccepted+=1 # number of non-accepted moves, resets to zero after accepted move.
        Niter+=1 # add one to the iteration number           
    if Niter>=Maxiter:
        print "exited due to max. number of iterations (%i) reached" %(Niter)
    else:
        print "Normal exit"
    print "Number of iterations per second",Niter/(time.time()-Now)
    print "Number of valid moves",Nmoves
    print "final Chi-squared value %f" %(Conval)
    Details['Niter']=Niter
    Details['Nmoves']=Nmoves
    Details['elapsed']=(time.time()-Now)

    Ifinal=sum(Iset,0)/sum(Vset**2)
    Sc,Conval=Iopt(I,Ifinal,E,Sc)    
    if OutputI:
        if OutputDetails:
            return Rset,(Ifinal*Sc[0]+Sc[1]),Conval,Details
        else:
            return Rset,(Ifinal*Sc[0]+Sc[1]),Conval
    else:
        if OutputDetails:
            return Rset,Conval,Details # ifinal cannot be output with variable length intensity outputs (in case of masked negative intensities or q limits)
        else:
            return Rset,Conval # ifinal cannot be output with variable length intensity outputs (in case of masked negative intensities or q limits)

######################################## end ###############################################

def McPlot(q,I,E,A,Histscale='log',AxisMargin=0.3):
    '''
    This function plots the output of the Monte-Carlo procedure in two windows, with the left window the measured signal versus the fitted intensity (on double-log scale), and the righthand window the size distribution
    '''

    #set plot font
    import matplotlib.font_manager as fm
    plotfont = fm.FontProperties(
                family = 'Courier New Bold', fname = '/Library/Fonts/Courier New Bold.ttf')
    textfont = fm.FontProperties( #Baskerville.ttc does not work when saving to eps
                family = 'Times New Roman', fname = '/Library/Fonts/Times New Roman.ttf')

    #initialize figure and axes
    fig=figure(figsize=(14,7),dpi=80,facecolor='w',edgecolor='k')
    q_ax=fig.add_subplot(121,axisbg=(.95,.95,.95),xlim=(numpy.min(q)*(1-AxisMargin),numpy.max(q)*(1+AxisMargin)),ylim=(numpy.min(I)*(1-AxisMargin),numpy.max(I)*(1+AxisMargin)))
    xlabel('q, 1/m',fontproperties=plotfont)
    ylabel('I, 1/(m sr)',fontproperties=plotfont)
    R_ax=fig.add_subplot(122,axisbg=(.95,.95,.95),xlim=(numpy.min(A['Hx'])*(1-AxisMargin),numpy.max(A['Hx'])*(1+AxisMargin)),ylim=(0,numpy.max(A['Hmean'])*(1+AxisMargin)))
    if Histscale=='log':
        xscale('log')
    xlabel('Radius, m',fontproperties=plotfont)
    ylabel('[Rel.] Volume Fraction',fontproperties=plotfont)
    fig.subplots_adjust(left=0.06,bottom=0.11,right=0.96,top=0.95,wspace=0.23,hspace=0.13)

    #set axis font and ticks
    R_ax.set_yticklabels(R_ax.get_yticks(), fontproperties = plotfont)
    R_ax.set_xticklabels(R_ax.get_xticks(), fontproperties = plotfont)
    q_ax.set_yticklabels(q_ax.get_yticks(), fontproperties = plotfont)
    q_ax.set_xticklabels(q_ax.get_xticks(), fontproperties = plotfont)
    #R_ax.spines['bottom'].set_color('black')
    R_ax.spines['bottom'].set_lw(2)
    R_ax.spines['top'].set_lw(2)
    R_ax.spines['left'].set_lw(2)
    R_ax.spines['right'].set_lw(2)
    R_ax.tick_params(axis='both',colors='black',width=2,which='major',direction='in',length=6)
    R_ax.tick_params(axis='x',colors='black',width=2,which='minor',direction='in',length=3)
    R_ax.tick_params(axis='y',colors='black',width=2,which='minor',direction='in',length=3)
    q_ax.spines['bottom'].set_lw(2)
    q_ax.spines['top'].set_lw(2)
    q_ax.spines['left'].set_lw(2)
    q_ax.spines['right'].set_lw(2)
    q_ax.tick_params(axis='both',colors='black',width=2,which='major',direction='in',length=6)
    q_ax.tick_params(axis='x',colors='black',width=2,which='minor',direction='in',length=3)
    q_ax.tick_params(axis='y',colors='black',width=2,which='minor',direction='in',length=3)

    #fill R axes
    axes(R_ax)
    bar(A['Hx'][0:-1],A['Hmean'],width=A['Hwidth'],color='orange',edgecolor='black',linewidth=1,zorder=2,label='MC size histogram')
    plot(A['Hmid'],A['vfminbins'],'r--',lw=5,label='Minimum visibility limit',zorder=3)
    R_eb=errorbar(A['Hmid'],A['Hmean'],A['Hstd'],zorder=4,fmt='k.',ecolor='k',elinewidth=2,capsize=4,ms=0,lw=2,solid_capstyle='round',solid_joinstyle='miter')
    legend(loc=1,fancybox=True,prop=textfont)
    title('Radius size histogram',fontproperties=textfont)
    #not quite sure if I can set the error bar linewidth on the bar plot in "bar" or whether I have to replot the errors themselves separately to get that level of control 

    #fill q axes
    axes(q_ax)
    eb=errorbar(q,I,E,zorder=2,fmt='k.',ecolor='k',elinewidth=2,capsize=4,ms=5,label='Measured intensity',lw=2,solid_capstyle='round',solid_joinstyle='miter')
    grid(lw=2,color='black',alpha=.5,dashes=[1,6],dash_capstyle='round',zorder=-1)
    xscale('log')
    yscale('log')
    aq=sort(A['q'])
    aI=A['Imean'][argsort(A['q'])]
    plot(aq,aI,'r-',lw=3,label='MC Fit intensity',zorder=4)
    plot(aq,numpy.mean(A['Screps'][1,:])+0*aq,'g-',linewidth=3,label='MC Background level:\n\t ({0:03.3g})'.format(numpy.mean(A['Screps'][1,:])),zorder=3)
    legend(loc=1,fancybox=True,prop=textfont)
    title('Measured vs. Fitted intensity',fontproperties=textfont)
    
    #Handles=dict
    #Handles['fig']=fig
    #Handles['R_ax']=R_ax
    #Handles['q_ax']=q_ax
    #return Handles
    return

def McCSV(filename,*args):
    '''
    This function writes a semicolon-separated csv file to [filename] containing an arbitrary number of columns of *args. in case of variable length columns, empty fields will contain ''.

    Input arguments should be lists or 1D arrays, which will be ouput as:
    arg1[0];arg2[0];arg3[0]
    arg1[1];arg2[1];arg3[1]
    etc.
    
    existing files with the same filename will be overwritten

    '''
    #uses sprintf rather than csv for flexibility
    ncol=len(args)
    #make format string used for every line, don't need this
    #linestr=''
    #for coli in range(ncol):
    #    linestr=linestr+'{'+'};'
    #linestr=linestr[0:-1]+'\n' #strip the last semicolon, add a newline


    #find out the longest row
    nrow=0
    for argi in range(len(args)):
        nrow=numpy.max((nrow,len(args[argi])))

    #now we can open the file:
    fh=open(filename,'w')
    emptyfields=0
    for rowi in range(nrow):
        linestr=''
        for coli in range(ncol):
            #print 'rowi {} coli {} len(args[coli]) {}'.format(rowi,coli,len(args[coli]))
            if len(args[coli])<=rowi: #we ran out of numbers for this arg
                linestr=linestr+';' #add empty field
                emptyfields+=1
            else:
                linestr=linestr+'{};'.format(args[coli][rowi])
        linestr=linestr[0:-1]+'\n'

        fh.write(linestr)

    fh.close()
    print '{} lines written with {} columns per line, and {} empty fields'.format(rowi,ncol,emptyfields)


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
        
    
#general functions
def csqr(Sc,I,Ic,E):
    #least-squares error for use with scipy.optimize.leastsq
    cs=(I-Sc[0]*Ic-Sc[1])/E
    #print "Size E",size(E)
    #print "Sc: %f, %f" %(Sc[0],Sc[1])
    return cs

def Iopt(I,Ic,E,Sc,OutputI=False):
    #new version, using leastsq and csqr, speed improvement of over a factor of 10 w.r.t. V1's use in the MC algorithm
    #optimizes the scaling and background factor to match Ic closest to I. returns an array with scaling factors. Input Sc has to be a two-element array wiht the scaling and background   
    #initial guesses. No bounds at the moment are applied, except that the intensity scaling has to be positive.
    Sc,success=scipy.optimize.leastsq(csqr,Sc,args=(I.flatten(),Ic.flatten(),E.flatten()),full_output=0)
    #print "Sc: %f, %f" %(Sc[0],Sc[1])
    cval=csqr_v1(I,Sc[0]*Ic+Sc[1],E)
    if OutputI:
        return Sc,cval,Sc[0]*Ic+Sc[1]
    else:
        return Sc,cval

def csqr_v1(I,Ic,E):
    #least-squares for data with known error, size of parameter-space not taken into account
    cs=sum(((I-Ic)/E)**2)/(size(I))
    return cs

def Iopt_v1(I,Ic,E,Sc,OutputI=False,Background=True):
    #old version, using fmin and csqr_v1
    #optimizes the scaling and background factor to match Ic closest to I. returns an array with scaling factors. Input Sc has to be a two-element array wiht the scaling and background   
    #initial guesses. No bounds at the moment are applied, except that the intensity scaling has to be positive.
    #Background can be set to False to just find the scaling factor.
    if Background:
        Sc=scipy.optimize.fmin(lambda Sc: csqr_v1(I,Sc[0]*Ic+Sc[1],E),Sc,full_output=0, disp=0)
        cval=csqr_v1(I,Sc[0]*Ic+Sc[1],E)
    else:
        Sc[1]=0.
        Sc=scipy.optimize.fmin(lambda Sc: csqr_v1(I,Sc[0]*Ic,E),Sc,full_output=0, disp=0)
        Sc[1]=0.
        cval=csqr_v1(I,Sc[0]*Ic,E)
    if OutputI:
        return Sc,cval,Sc[0]*Ic+Sc[1]
    else:
        return Sc,cval

#########################################DEAD CODE?################################################

#older code, may need to be updated before it can or should be used again.

#Monte-carlo procedure for cylinders
def MCFit_cyl_RadialIsotropic(q,I,E,Ncyl=200,Bounds=[],Convcrit=1,Rpfactor=1.,Maxiter=1e5,R1Prior=[],R2Prior=[],Oprior=[],Qlimits=numpy.array([]),MaskNegI=False,OutputI=False,R2IsAspect=False,OutputDetails=False):
    #Bound has to be a 4-element vector, with bounds for R1, R2
    #rewrite of the monte-carlo method previously implemented in Matlab. Simpler form, but open source might mean slight improvements.
    #initialise parameters
    #Convcrit=1 #reasonable value for poisson weighting. any lower than this and we would be fitting noise
    #OutputI can be set to True only if MaskNeg is not True
    Weighting='Poisson' #only one implemented
    Method='Randmove' #only one implemented
    if MaskNegI==True:
        OutputI=False
    if (size(Qlimits)==0): #no q-limits supplied
        Qlimits=numpy.min(q)
    if (size(Qlimits)==1): #only lower q limit supplied
        if isinstance(Qlimits,int): #integer supplied, removing a number of values
            Qlimits=q[Qlimits]
        Qlimits=append(Qlimits,numpy.max(q))
    if size(Qlimits)==2: #make sure they are in the right order
        if isinstance(Qlimits[0],int): #integer supplied, removing a number of values
            Qlimits[0]=q[Qlimits[0]]
        if isinstance(Qlimits[1],int): #integer supplied, removing a number of values
            Qlimits[1]=q[Qlimits[1]]
        Qlimits=numpy.array([numpy.min(Qlimits),numpy.max(Qlimits)])
    #apply limits
    I=I[(q>=Qlimits[0])&(q<=Qlimits[1])]
    E=E[(q>=Qlimits[0])&(q<=Qlimits[1])]
    q=q[(q>=Qlimits[0])&(q<=Qlimits[1])]
    #optional masking of negative intensity
    if MaskNegI==True:
        q=q[I>=0]
        E=E[I>=0]
        I=I[I>=0]
        
    if size(Bounds)==0: #if the bounds are not supplied, make a good guess
        Bounds[0:2]=array([pi/numpy.max(q),pi/numpy.min(q)]) #reasonable, but not necessarily correct, parameters
        Bounds[2:4]=Bounds[0:2]
        print 'Bounds not provided, so set related to minimum and maximum q. Lower and upper bounds are {0} and {1} for R1 as well as R2'.format(Bounds[0],Bounds[1])

    #Rpower=2 #squared results in most reasonable end result. none=3, linear=2.5, squared=2, cubed=1.5, fourth=1. This reduces the volume dependency in the calculation leading to quicker convergence
    #Maxiter=1e6

    #intialise variables
    FFsqrset=[]
    Vset=[]
    Niter=0
    Conval=inf
    Ri=0 #index of sphere to change. We'll sequentially change spheres, which is perfectly random since they are in random order.
    
    #generate initial set of spheres
    if size(R1Prior)==0:
        R1set=numpy.random.uniform(numpy.min(Bounds[0:2]),numpy.max(Bounds[0:2]),Ncyl)
        R2set=numpy.random.uniform(numpy.min(Bounds[2:4]),numpy.max(Bounds[2:4]),Ncyl)
        Oset=numpy.random.uniform(numpy.min(Bounds[4:6]),numpy.max(Bounds[4:6]),Ncyl)
    elif (size(R1Prior)!=0)&(size(Ncyl)==0):
        Ncyl=size(R1Prior)
        R1set=R1Prior
        R2set=R2Prior
        Oset=OPrior
    elif size(R1Prior)==Ncyl:
        R1set=R1Prior
        R2set=R2Prior
        Oset=OPrior
    elif size(R1Prior)<Ncyl:
        print "size of prior is smaller than Ncyl. duplicating random prior values"
        #while size(Prior)<Nsph:
        Addi=numpy.random.randint(size(R1Prior),size=Ncyl-size(R1Prior))
        R1set=concatenate((R1Prior,R1Prior[Addi]))
        R2set=concatenate((R2Prior,R2Prior[Addi]))
        Oset=concatenate((OPrior,OPrior[Addi]))
        print "size now:", size(R1set)
    elif size(R1Prior)>Ncyl:
        print "Size of prior is larger than Ncyl. removing random prior values"
        Remi=numpy.random.randint(size(R1Prior),size=Ncyl) #remaining choices
        R1set=R1Prior[Remi]
        R2set=R2Prior[Remi]
        Oset=OPrior[Remi]
        print "size now:", size(R1set)
    
    Arange=array(range(Ncyl)) #indices to array values, we'll need this later for some logical indexing operations

    #calculate their form factors
    #FFset=FF_sph_1D(q,Rset)
    FFsqrset=zeros((size(R1set),size(q)))
    print "Generating initial set of cylinders..."
    now=time.time()
    for Ri in range(size(R1set)):
        if R2IsAspect:
            FFsqrset[Ri,:]=Integrate_Shape_as_Radial_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1set[Ri],R2set[Ri]*R1set[Ri],Oset[Ri]]))
        else:
            FFsqrset[Ri,:]=Integrate_Shape_as_Radial_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1set[Ri],R2set[Ri],Oset[Ri]]))
    print "average time (s) per cylinder: {0}".format((time.time()-now)/size(R1set))

    if R2IsAspect:
        Vset=pi*R1set**(2*Rpfactor)*2*(R2set*R1set)**(Rpfactor)
    else:
        Vset=pi*R1set**(2*Rpfactor)*2*R2set**(Rpfactor)
    #calculate the intensities
    Iset=FFsqrset*(Vset[:,newaxis]+0*FFsqrset)**2 #a set of intensities
    Vst=sum(Vset**2) #total volume squared
    It=sum(Iset,0) #the total intensity
    #optimize the intensities and calculate convergence criterium
    Sc,Conval1=Iopt_v1(I,It/Vst,E,[1,1]) #V1 is more robust w.r.t. a poor initial guess
    Sc,Conval=Iopt(I,It/Vst,E,Sc) #reoptimize with V2, there might be a slight discrepancy in the residual definitions of V1 and V2 which would prevent optimization.
    #print "Initial conval V1",Conval1
    print "Initial conval V2",Conval
    #start the MC procedure
    Now=time.time()
    Nmoves=0 #tracking the number of moves
    while (Conval>Convcrit) &(Niter<Maxiter):
        R1t=numpy.random.uniform(numpy.min(Bounds[0:2]),numpy.max(Bounds[0:2]),1)
        R2t=numpy.random.uniform(numpy.min(Bounds[2:4]),numpy.max(Bounds[2:4]),1)
        Ot=numpy.random.uniform(numpy.min(Bounds[4:6]),numpy.max(Bounds[4:6]),1)
        if R2IsAspect:
            Fst=Integrate_Shape_as_Radial_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1t,R2t*R1t,Ot]))
            Vtt=pi*R1t**(2*Rpfactor)*2*(R1t*R2t)**(Rpfactor)
        else:
            Fst=Integrate_Shape_as_Radial_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1t,R2t,Ot]))
            Vtt=pi*R1t**(2*Rpfactor)*2*R2t**(Rpfactor)
        Itt=(Fst*Vtt**2)
        #calculate new total intensity
        Itest=(It-Iset[Ri,:]+Itt) #we do subtractions and additions, which give us another factor 2 improvement in speed over summation and is much more scalable
        Vstest=(sqrt(Vst)-Vset[Ri])**2+Vtt**2
        #optimize intensity and calculate convergence criterium
        Sct,Convalt=Iopt(I,Itest/Vstest,E,Sc)#using version two here for a >10 times speed improvement
        #test if the radius change is an improvement:
        if Convalt<Conval: #it's better
            R1set[Ri],R2set[Ri],Oset[Ri],Iset[Ri,:],It,Vset[Ri],Vst,Sc,Conval=(R1t,R2t,Ot,Itt,Itest,Vtt,Vstest,Sct,Convalt)
            print "Improvement in iteration number %i, convergence value %f of %f\r" %(Niter,Conval,Convcrit),
            Nmoves+=1
        #else nothing to do
        Ri+=1 #move to next sphere in list
        Ri=Ri%(Ncyl) #loop if last sphere
        Niter+=1 #add one to the iteration number           
    if Niter>=Maxiter:
        print "exited due to max. number of iterations (%i) reached" %(Niter)
    else:
        print "Normal exit"
    print "Number of iterations per second",Niter/(time.time()-Now)
    print "Number of valid moves",Nmoves
    print "final convergence value %f" %(Conval)
    Ifinal=sum(Iset,0)/sum(Vset**2)
    Sc,Conval=Iopt(I,Ifinal,E,Sc)    
    Details=dict()
    Details['Niter']=Niter
    Details['Nmoves']=Nmoves
    Details['elapsed']=(time.time()-Now)

    if OutputI:
        if OutputDetails:
            return R1set,R2set,Oset,(Ifinal*Sc[0]+Sc[1]),Conval,Details
        else:
            return R1set,R2set,Oset,(Ifinal*Sc[0]+Sc[1]),Conval
    else:
        if OutputDetails:
            return R1set,R2set,Oset,Conval,Details #ifinal cannot be output with variable length intensity outputs (in case of masked negative intensities or q limits)
        else:
            return R1set,R2set,Oset,Conval #ifinal cannot be output with variable length intensity outputs (in case of masked negative intensities or q limits)

def Analyze_1D_Cylinder_RadiallyIsotropic(q,I,E,Bounds=[],Ncyl=[],Maxiter=[],Rpfactor=1.,Nreps=100,qlims=[0,inf],Histbins=50,Histscale='log',drhosqr=1,Convcrit=1.,Maxntry=7,R2IsAspect=False):
    #this function runs the monte-carlo fit several times, and returns bin centres, means, standard deviations and observability limits for the bins. Eventually, this may also plot. The drhosqr value can be set to the contrast (if known) so a volume fraction for each histogram bin is calculated.
    #Rpower = 3 when Rpfactor = 1. Rpower = 3*Rpfactor
    #for volume weighting, Rpfactor = 1.5/3
    #Rpower = Rpfactor*3
    if Histscale=='lin':
        Hx1 = linspace(Bounds[0],Bounds[1],Histbins+1)
        if R2IsAspect:
            Hx2 = linspace(Bounds[2]*Bounds[0],Bounds[3]*Bounds[1],Histbins+1)
        else:
            Hx2 = linspace(Bounds[2],Bounds[3],Histbins+1)
    else:
        Hx1 = 10**(linspace(log10(Bounds[0]),log10(Bounds[1]),Histbins+1))
        if R2IsAspect:
            Hx2 = 10**(linspace(log10(Bounds[2]*Bounds[0]),log10(Bounds[3]*Bounds[1]),Histbins+1))
        else:
            Hx2 = 10**(linspace(log10(Bounds[2]),log10(Bounds[3]),Histbins+1))
    #Rmc,ConVal = MCFit_sph(q,I,E,Bounds=Bounds,Nsph=Nsph,Maxiter=Maxiter,Rpower=Rpower)
    #if ConVal>1:
    #    print "test optimisation not converging to below 1! fix first, then ask again"
    #    return

    #q and psi limits
    validbools = (q>=qlims[0])&(q<=qlims[1])
    I = I[validbools]
    if size(E)!=0:
        E = E[validbools]
    q = q[validbools]

    R1rep = zeros([Ncyl,Nreps])
    R2rep=zeros([Ncyl,Nreps])
    Orep=zeros([Ncyl,Nreps])
    Vf=zeros([Ncyl,Nreps])
    Vft=zeros([Nreps])
    #Rrep[:,0]=Rmc
    Hy1=zeros([Histbins,Nreps])
    Hy2=zeros([Histbins,Nreps])
    #Hr=hist(Rmc,Hx,normed=True,weights=1/Rmc**((3-Rpower)*2))
    #Hy[:,0]=Hr[0]
    Irep=zeros([len(I),Nreps])
    #try optimizing Nreps times
    bignow=time.time()
    Niters=zeros([Nreps])

    for nr in arange(0,Nreps):
        nt=0
        #do the thing!
        R1rep[:,nr],R2rep[:,nr],Orep[:,nr],Irep[:,nr],ConVal,Details=MCFit_cyl_RadialIsotropic(q,I,E,Ncyl=Ncyl,Bounds=Bounds,Convcrit=Convcrit,Rpfactor=Rpfactor,Maxiter=Maxiter,OutputI=True,R2IsAspect=R2IsAspect,OutputDetails=True)
        #if the optimization fails, perhaps it was an uncommon failure, so we retry max 10 times to get convergence.
        while ConVal>Convcrit:
            nt+=1
            R1rep[:,nr],R2rep[:,nr],Orep[:,nr],Irep[:,nr],ConVal,Details=MCFit_cyl_RadialIsotropic(q,I,E,Ncyl=Ncyl,Bounds=Bounds,Convcrit=Convcrit,Rpfactor=Rpfactor,Maxiter=Maxiter,OutputI=True,R2IsAspect=R2IsAspect,OutputDetails=True)
            if nt>Maxntry:
                print "could not reach optimization criterion within {0} attempts, exiting...".format(Maxntry+2)
                return

        Niters[nr]=Details['Niter']
        #absolute intensity calculation-->
        #calculate volume contribution
        if R2IsAspect:
            Vsa=pi*R1rep[:,nr]**(2*Rpfactor)*2*(R1rep[:,nr]*R2rep[:,nr])**(Rpfactor)
            #real particle volume:
            Vpa=pi*R1rep[:,nr]**(2)*2*(R1rep[:,nr]*R2rep[:,nr])
        else:
            Vsa=pi*R1rep[:,nr]**(2*Rpfactor)*2*R2rep[:,nr]**(Rpfactor)
            #real particle volume:
            Vpa=pi*R1rep[:,nr]**(2)*2*R2rep[:,nr]
        #calculate intensities without scaling:

        FFsqrset=zeros((size(R1rep[:,nr]),size(q)))
        print "Generating intensity of set of cylinders..."
        now=time.time()
        for Ri in range(size(R1rep[:,nr])):
            if R2IsAspect:
                FFsqrset[Ri,:]=Integrate_Shape_as_Radial_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1rep[Ri,nr],R1rep[Ri,nr]*R2rep[Ri,nr],Orep[Ri,nr]]))
            else:
                FFsqrset[Ri,:]=Integrate_Shape_as_Radial_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1rep[Ri,nr],R2rep[Ri,nr],Orep[Ri,nr]]))
        #calculate the intensities
        Ics=FFsqrset*(Vsa[:,newaxis]+0*FFsqrset)**2 #a set of intensities
        Icyls=sum(Ics,0) #the total intensity

        #determine scaling factor for each contribution:
        Sci=numpy.max(I)/numpy.max(Icyls)
        Sc,Cv=Iopt(I,Icyls,E,[Sci,1])
        print "Scaling: {0}, background {1}".format(Sc[0],Sc[1])
        #now for each sphere, calculate its volume fraction:
        sR=size(R1rep,0)
        Vf[:,nr]=Sc[0]*Vsa**2/(Vpa*drhosqr) #a set of intensities
        Vft[nr]=sum(Vf[:,nr]) #total volume squared
        print "Calculating volfrac for this iteration: {0}".format(Vft[nr])
        #<-- absolute intensity calculation
        #Hr1=hist(R1rep[:,nr],Hx1,normed=False,weights=Vf[:,nr])
        Hr1=numpy.histogram(R1rep[:,nr],Hx1,density=False,weights=Vf[:,nr])
        #Hr2=hist(R2rep[:,nr],Hx2,normed=False,weights=Vf[:,nr])
        Hr2=numpy.histogram(R2rep[:,nr],Hx2,density=False,weights=Vf[:,nr])
        #old, non-volume fraction: Hr=hist(Rrep[:,nr],Hx,normed=False,weights=1/Rrep[:,nr]**((3-Rpower)*2))
        Hy1[:,nr]=Hr1[0]
        Hy2[:,nr]=Hr2[0]
        biglap=time.time()
        #in minutes:
        tottime=(biglap-bignow)/60.
        avetime=(tottime/(nr+1))
        remtime=(avetime*Nreps-tottime)
        print "\t*finished optimization number {0} of {1} \r\n\t*total elapsed time: {2} minutes \r\n\t*average time per optimization {3} minutes \r\n\t*total time remaining {4} minutes".format(nr+1,Nreps,tottime,avetime,remtime)

    Imean=numpy.mean(Irep,axis=1)
    Istd=numpy.std(Irep,axis=1)
    H1mid,H2mid=zeros(Histbins),zeros(Histbins)
    for bini in range(Histbins):
        H1mid[bini]=numpy.mean(Hx1[bini:bini+2])
        H2mid[bini]=numpy.mean(Hx2[bini:bini+2])
    H1mean=numpy.mean(Hy1,axis=1)
    H2mean=numpy.mean(Hy2,axis=1)
    H1std=numpy.std(Hy1,axis=1)
    H2std=numpy.std(Hy2,axis=1)
    obs,qmo=observability(q,H1mid)
    iobs=1/obs/obs[-1]*H1mean[-1] #normalize inverse observability to amount in the last bin


    #store in output dict
    A=dict()
    A['H1width']=diff(Hx1)
    A['H1mid']=H1mid
    A['H1mean']=H1mean
    A['H1std']=H1std
    A['R1rep']=R1rep
    A['Hy1']=Hy1

    A['H2width']=diff(Hx2)
    A['H2mid']=H2mid
    A['H2mean']=H2mean
    A['H2std']=H2std
    A['R2rep']=R2rep
    A['Hy2']=Hy2

    A['iobs']=iobs
    A['Imean']=Imean
    A['Istd']=Istd
    A['Vf']=Vf
    A['Vft']=Vft
    A['q']=q

    A['Niter']=numpy.mean(Niters)

    return A

def Icyl_InPlaneAverage_RandomSampling(q,psi=array([0.,90.]),nsamples='auto',R1=1.,R2=1.):
    """
    calculates the in-plane average of a shape (rotational average, but only rotated around the beam axis, perpendicular to the detector plane).
    This calculation works by uniform random number generation for q and psi, within the bounds dictated by the input. 
    Set "nsamples" to a number of samples. 'auto' sounds cool but has not been implemented yet. 
    input 'q' is supposed to be a vector of q values for which the intensity is requested
    """

    Qs=numpy.random.uniform(low=numpy.min(abs(q)),high=numpy.max(abs(q)),size=nsamples)
    PSIs=numpy.random.uniform(low=numpy.min(abs(psi)),high=numpy.max(abs(psi)),size=nsamples)
    Is=Icyl_2D(Qs,PSIs,R1,R2)

    qbin,Ibin,SEbin=binning_weighted_1D(Qs,Is,Nbins=q,Stats='SE')
    return qbin,Ibin,SEbin
    #this is not really an efficient way of doing it, requires many points


def Integrate_Shape_as_Spherical_Isotropic(q,psirange=numpy.array([0.1,360.1]),Func='',parameters=numpy.array([0]),psidiv=303,randfact=0):
    """
    Generate a 1D intensity of a spherically isotropic shape. 
    This is identical to a fully rotationally isotropic averaged shape, for radially isotropic shapes, use Integrate_Shape_as_Radial_Isotropic, which averages the shape around the axis perpendicular to the detector.
    """
    d_to_r=1./180*pi
    psi=linspace(numpy.min(psirange),numpy.max(psirange),psidiv)
    Q=(q+0*psi[:,newaxis])
    PSI=(0*q+psi[:,newaxis])
    if randfact>0:
        PSI=PSI+numpy.random.uniform(high=diff(psirange)/(randfact*psidiv),size=shape(PSI))
    if Func=='FF_ell_2D':
        FFc=FF_ell_2D(Q,PSI,parameters[0],parameters[1],rot=parameters[2])
        FFcsqr=numpy.mean(FFc**2*abs(sin(PSI*d_to_r)),axis=0)
        #FFcsqr=1/diff(psirange)*numpy.trapz(FFc**2,axis=0,x=psi)
        return FFcsqr
    elif Func=='FF_cyl_2D':
        FFc=FF_cyl_2D(Q,PSI,parameters[0],parameters[1],rot=parameters[2])
        FFcsqr=numpy.mean(FFc**2*abs(sin(PSI*d_to_r)),axis=0)
        #FFcsqr=1/diff(psirange)*numpy.trapz(FFc**2,axis=0,x=psi)
        return FFcsqr

def Integrate_Shape_as_Radial_Isotropic(q,psirange=numpy.array([0.1,360.1]),Func='',parameters=numpy.array([0]),psidiv=303,randfact=0):
    """
    Generate a 1D intensity of a radially isotropic shape. 
    This is not identical to a fully rotationally isotropic shape, but merely isotropic as if the shape 
    is rotated around the axis perpendicular to the detector.
    """
    psi=linspace(numpy.min(psirange),numpy.max(psirange),psidiv)
    Q=(q+0*psi[:,newaxis])
    PSI=(0*q+psi[:,newaxis])
    if randfact>0:
        PSI=PSI+numpy.random.uniform(high=diff(psirange)/(randfact*psidiv),size=shape(PSI))
    if Func=='FF_cyl_2D':
        FFc=FF_cyl_2D(Q,PSI,parameters[0],parameters[1],rot=parameters[2])
        FFcsqr=numpy.mean(FFc**2,axis=0)
        #FFcsqr=1/diff(psirange)*numpy.trapz(FFc**2,axis=0,x=psi)
        return FFcsqr
    elif Func=='Icyl_2D':
        Ic=Icyl_2D(Q,PSI,parameters=parameters)
        Ic=numpy.mean(Ic,axis=0)
        return Ic

def Icyl_2D(Q,PSI,R1cyl=0,R2cyl=0,parameters=''):
    #This function calculates the scattering intensity of a cylinder with radius R1cyl and height 2*R2cyl
    if parameters!='':
        R1cyl=parameters[0]
        R2cyl=parameters[1]
        rot=parameters[2]
    PSIr=(PSI-rot)/180*pi
    qRsina=Q*R1cyl*sin(PSIr)
    qLcosa=Q*R2cyl*cos(PSIr)
    FF=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )
    Vc=pi*R1cyl**2*2*R2cyl
    return abs(FF)**2 * Vc**2

def FF_ell_2D(q,psi,R1,R2,rot):
    #R1<R2, prolate ellipsoid (cigar-shaped), R1>R2, oblate ellipsoid (disk-shaped), rotation is offset from perfect orientation (psi-rot)
    d_to_r=1./360*2*pi #degrees to radians, forget the dot and get yourself into a non-floating point mess, even though pi is floating point...
    if size(R1)==1:
        sda=sin((psi-rot)*d_to_r)
        cda=cos((psi-rot)*d_to_r)
        r=sqrt(R1**2*sda**2+R2**2*cda**2)
        qr=q*r
        Fell=3*(sin(qr)-qr*cos(qr))/(qr**3)
    else: #calculate a series
        if len(q)==1: #one-dimensional q and psi, add another dimension
            q=q[:,newaxis]
            psi=psi[:,newaxis]
        Fell=zeros([size(q,0),size(q,1),size(R1)])
        for Ri in range(size(R1)):
            sda=sin((psi-rot[Ri])*d_to_r)
            cda=cos((psi-rot[Ri])*d_to_r)
            r=sqrt(R1[Ri]**2*sda**2+R2[Ri]**2*cda**2)
            qr=q*r
            Fell[:,:,Ri]=3*(sin(qr)-qr*cos(qr))/(qr**3)

    return Fell #this will be a three-dimensional result for ranges of Ri, and one- to two-dimensional for single radii sets
    
def FF_cyl_2D(Q,PSI,R1,R2,rot=0.):
    #Cylinder form factor, rotation is offset from perfect orientation (psi-rot)
    d_to_r=1./180*pi #degrees to radians, forget the dot and get yourself into a non-floating point mess, even though pi is floating point...
    if size(R1)==1:
        PSIr=(PSI-rot)*d_to_r
        qRsina=Q*R1*sin(PSIr)
        qLcosa=Q*R2*cos(PSIr)
        Fcyl=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )
    else: #calculate a series
        if ndim(q)==1: #one-dimensional q and psi, add another dimension
            Q=Q[:,newaxis]
            PSI=PSI[:,newaxis]
        Fell=zeros([size(Q,0),size(Q,1),size(R1)])
        for Ri in range(size(R1)):
            PSIr=(PSI-rot[Ri])*d_to_r
            qRsina=Q*R1[Ri]*sin(PSIr)
            qLcosa=Q*R2[Ri]*cos(PSIr)
            Fcyl[:,:,Ri]=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )

    return Fcyl #this will be a three-dimensional result for ranges of Ri, and one- to two-dimensional for single radii sets

#dead code, may be deleted.
    
#useful, but not in this set of functions:

def Icyl_1D_SphericallyIsotropic(q,R1set=[],R2set=[],OBounds=numpy.array([0,0.01]),Rpfactor=1.,R2IsAspect=False):
    #generates intensity of a set of randomly oriented cylinders with radius R1set and length R2set*2.
        
    #intialise variables
    FFsqrset=[]
    Vset=[]
    
    Ncyl=len(R1set)
    #generate set of spheres
    Oset=numpy.random.uniform(numpy.min(OBounds),numpy.max(OBounds),Ncyl)

    #calculate their form factors
    #FFset=FF_sph_1D(q,Rset)
    FFsqrset=zeros((size(R1set),size(q)))
    print "Generating initial set of cylinders..."
    now=time.time()
    for Ri in range(size(R1set)):
        if R2IsAspect:
            FFsqrset[Ri,:]=Integrate_Shape_as_Spherical_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1set[Ri],R2set[Ri]*R1set[Ri],Oset[Ri]]))
        else:
            FFsqrset[Ri,:]=Integrate_Shape_as_Spherical_Isotropic(q,Func='FF_cyl_2D',parameters=numpy.array([R1set[Ri],R2set[Ri],Oset[Ri]]))
    print "average time (s) per cylinder: {0}".format((time.time()-now)/size(R1set))

    if R2IsAspect:
        Vset=pi*R1set**(2*Rpfactor)*2*(R2set*R1set)**(Rpfactor)
    else:
        Vset=pi*R1set**(2*Rpfactor)*2*R2set**(Rpfactor)
    #calculate the intensities
    Iset=FFsqrset*(Vset[:,newaxis]+0*FFsqrset)**2 #a set of intensities
    It=sum(Iset,0) #the total intensity
    return It

def Iell_1D_SphericallyIsotropic(q,R1set=[],R2set=[],OBounds=numpy.array([0,0.01]),Rpfactor=1.,R2IsAspect=False):
    #generates intensity of a set of randomly oriented cylinders with radius R1set and length R2set*2.
        
    #intialise variables
    FFsqrset=[]
    Vset=[]
    
    Nell=len(R1set)
    #generate set of spheres
    Oset=numpy.random.uniform(numpy.min(OBounds),numpy.max(OBounds),Nell)

    #calculate their form factors
    #FFset=FF_sph_1D(q,Rset)
    FFsqrset=zeros((size(R1set),size(q)))
    print "Generating initial set of ellipsoids..."
    now=time.time()
    for Ri in range(size(R1set)):
        if R2IsAspect:
            FFsqrset[Ri,:]=Integrate_Shape_as_Spherical_Isotropic(q,Func='FF_ell_2D',parameters=numpy.array([R1set[Ri],R2set[Ri]*R1set[Ri],Oset[Ri]]))
        else:
            FFsqrset[Ri,:]=Integrate_Shape_as_Spherical_Isotropic(q,Func='FF_ell_2D',parameters=numpy.array([R1set[Ri],R2set[Ri],Oset[Ri]]))
    print "average time (s) per cylinder: {0}".format((time.time()-now)/size(R1set))

    if R2IsAspect:
        Vset=4./3*pi*R1set**(2*Rpfactor)*(R2set*R1set)**(Rpfactor)
    else:
        Vset=4./3*pi*R1set**(2*Rpfactor)*R2set**(Rpfactor)
    #calculate the intensities
    Iset=FFsqrset*(Vset[:,newaxis]+0*FFsqrset)**2 #a set of intensities
    Vst=sum(Vset**2) #total volume squared
    It=sum(Iset,0) #the total intensity
    return It

def observability(q,Rset,I=[],E=[],Rpfactor=3):
    #observability calculation for a series of spheres, over a range of q. Additional intensity and errors may be supplied for error-weighted observability. Intensity is used for determining the intesity scaling and background levels.
    FFset=FF_sph_1D(q,Rset)
    Vset=(4.0/3*pi)*Rset**(3*Rpfactor)
    #calculate the intensities
    Iset=FFset**2*(Vset[:,newaxis]+0*FFset)**2 #a set of intensities
    Vst=sum(Vset**2) #total volume squared
    It=sum(Iset,0) #the total intensity
    #optimize the intensities and calculate convergence criterium
    ov=zeros(size(Rset))
    qm=zeros(size(Rset))
    if (size(I)==size(E))&(size(I)!=0):
        Sc,Conval1=Iopt_v1(I,It/Vst,E,[1,1]) #V1 is more robust w.r.t. a poor initial guess
        Sc,Conval=Iopt(I,It/Vst,E,Sc) #reoptimize with V2, there might be a slight discrepancy in the residual definitions of V1 and V2 which would prevent optimization.
        print "Conval V1",Conval1
        print "Conval V2",Conval
        print Sc
        for isi in range(size(Iset,0)):
            ov[isi]=sum(Sc[0]*Iset[isi,:]/(size(I)*E))
            #ov[isi]=numpy.max(Iset[isi,:]/(size(I)*Sc[0]*E))
            qm[isi]=q[numpy.argmax(Sc[0]*Iset[isi,:]/(size(I)*E))]
    else:
        for isi in range(size(Iset,0)):
            ov[isi]=numpy.max(Iset[isi,:]/It)
            qm[isi]=q[numpy.argmax(Iset[isi,:]/It)]
    return ov,qm
