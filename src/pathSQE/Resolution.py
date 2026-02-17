import numpy as np
import matplotlib.pyplot as plt


hbar = 6.58211899e-16*1e12*1e3 #eV*s * ps/s * meV/eV = meV*ps
kB = 8.61734315e-5*1e3  #ev/K * meV/eV
hbar_meVs = 6.582119*1e-13 # meV * s
VasptoTHz=15.633302   #convert eV-Angstroms to THz
THz2meV=4.1357   #convert THz to millielectronvolts
factors=14.400   #have no idea what this is for
m_neut=1.674927e-27 #kg
hbarJoules=1.05457173e-34 #m^2 kg/s
joulesTomeV = 6.24150647996E+21
epsilon = 1e-10 #very small number to use in numerical comparisons

class Instrument :
    """Class to define useful instrumental parameters for resolution convolution informationp.
    Takes instName (HERIX, CNCS, ARCS, SEQ)
    FWHM (estimated FWHM resolution at elastic line)
    IncidentEnergy (incident energy used for the experiment).
    """
    
    
    def __init__(self,instName,FWHM,IncidentEnergy) :
        self.instName = instName,
        self.FWHM = FWHM
        self.ei = IncidentEnergy
        
        
        if instName == 'HERIX' :
            self.type = 'IXS'
    

        if instName == 'CNCS' :

            self.type='INS'
            self.L1=28.3 #m "L1" mod face to sample chopper
            self.L2=1.487 #m "L2" sample chopper to Sample
            self.L3=3.5 #m variable to 3.5 m.... "L3" Sample to detector
            self.dtp=0.00000103 #pulse width at pulse shaping chopper; use FWHM at mon nearest
            self.dtd=0.000014  #sample spatial uncertainty in m
            self.dtc=0.00000026 #chopper uncertainty in s; use FWHM at mon nearest

        if instName == 'SEQ' :

            self.type='INS'
            self.L1=18.0 #m "L1" mod face to Fermi
            self.L2=2.0 #m "L2" Fermi to Sample
            self.L3=5.5 #m variable to 5.5 m.... "L3" Sample to detector
            self.dtp=0.000012 #pulse width at pulse shaping chopper; use FWHM at mon nearest
            self.dtd=0.000024  #sample spatial uncertainty in m
            self.dtc=0.000036 #chopper uncertainty in s; use FWHM at mon nearest


        if instName == 'ARCS' : 
            self.type = 'INS'
            self.L1=11.6 #m "L1" mod face to Fermi
            self.L2=1.99 #m "L2" Fermi to Sample
            self.L3=3.0 #m variable to 3.5 m.... "L3" Sample to detector
            self.dtp=0.0000012 #pulse width at pulse shaping chopper in seconds
            self.dtd=0.0000024 #sample spatial uncertainty in meters
            self.dtc=0.0000036 #chopper uncertainty in seconds
            

    
    def setPulseWidth(self,pw):
        """setPulseWidth(pw) : setter function for neutron beam pulse width at pulse shaping chopper. 
         Value should be given in seconds."""
        self.dtp = pw
    def setChopperUncertainty(self,ci):
        """setChopperUncertainty(pw) : setter function for neutron beam pulse width at due to Fermi chopper. 
         Value should be given in seconds."""
        self.dtc = ci
    def setSampleUncertainty(self,si):
        """setSampleUncertainty(pw) : setter function for timing uncertainty due to spatial extent. 
        Value should be given in seconds."""
        self.dtd = si


class Resolution :
    def __init__(self,type,inst,sigmaq=None,poly=None) :
        """Class to perform resolution function calculations.  
        Wrapper to set of subclasses which hold the correct calculating functions.
        Determination of which class is used is based on the inputs passed inp.
        'constant' - calculates a constant energy resolution using the given instrument and FWHM.
        'polynomial' -- uses a fourth order polynomial fit to a windsor-type calculation provided by the insrument scientist.
        
        Other inputs:
        inst - an instants of the Instrument object which contains instrumental parameters that might be needed for the calculationp.
        sigmaq - estimated Q resolution provided in the SimphoniesInput.py file (constant values only supported at this point.)
        poly - polynomial coefficients as described above.
        """
        self._name = type
        self.inst = inst 
        self.sigmaq = sigmaq
        self.poly = poly
        
        if (type == 'constant' and sigmaq!=None):
            self.__class__ = ConstantRes

        if (type == 'constant' and sigmaq==None):
            self.__class__ = ConstantResNoQResolution


        if (type == 'polynomial' and sigmaq!=None):
            self.__class__ = PolynomialRes

        if (type == 'polynomial' and sigmaq==None):
            self.__class__ = PolynomialResNoQResolution

            
        if (type == 'instrument' and sigmaq!=None):
            self.__class__ = InstrumentResolution

        if (type == 'instrument' and sigmaq==None):
            self.__class__ = InstrumentResolutionNoQResolution
          

class ConstantRes(Resolution):
    """creates constant resolution subclass of Resolution Functionp.  
    Takes a FWHM value provided either from user input 
    as ElasticFWHM or taken from the Instrument Class."""
    def sigma(self):
        sig = self.inst.FWHM/2.355
        return sig
    
    def Gauss(self,evalues,qvalues,ecenter,qcenter):
        Ne=1.0/np.sqrt(2*np.pi*self.sigma()**2)
        Nq=1.0/np.sqrt(2*np.pi*self.sigmaq**2)
        delta_q=qvalues-qcenter
        
       
        q_square=np.dot(delta_q,delta_q)
        
        
        Res= Ne*np.exp(-(evalues-ecenter)**2/(2*self.sigma()**2))*Nq*np.exp(-(q_square)/(2*self.sigmaq**2))
        return Res

class ConstantResNoQResolution(Resolution):
    """creates constant resolution subclass of Resolution Function with no Q resolutionp..  
    Takes a FWHM value provided either from user input 
    as ElasticFWHM or taken from the Instrument Class."""
    
    
    def sigma(self):
        sig = self.inst.FWHM/2.355
        return sig
    
    def Gauss(self,evalues,qvalues,ecenter,qcenter):
        Ne=1.0/np.sqrt(2*np.pi*self.sigma()**2)
        #Nq=1.0/np.sqrt(2*np.pi*self.sigmaq**2)
        #delta_q=qvalues-qcenter
        #q_square=np.dot(delta_q,delta_q)
        Res= Ne*np.exp(-(evalues-ecenter)**2/(2*self.sigma()**2))*Nq*np.exp(-(q_square)/(2*self.sigmaq**2))
        return Res


class PolynomialRes(Resolution):
    """creates polynomial resolution subclass of Resolution Functionp.  
    Takes a fourth order polynomial coefficient values provided either from user input 
    or taken from the Instrument Class."""
    def sigma(self,E):
        """ Energy resolution function for an instrument using coefficients of a characteristic polynomial.
         These values are typically supplied by an instrument scientist.
         """
        sig = (self.poly[0] + self.poly[1]*E +self.poly[2]*E**2 + self.poly[3]*E**3 +self.poly[4]*E**4)/2.35 #convert to unit to c
        return sig
    
    def Gauss(self,evalues,qvalues,ecenter,qcenter):
        Ne=1.0/np.sqrt(2*np.pi*self.sigma(ecenter)**2)
        Nq=1.0/np.sqrt(2*np.pi*self.sigmaq**2)
        delta_q=qvalues-qcenter
        q_square=np.dot(delta_q,delta_q)
        
        Res= Ne*np.exp(-(evalues-ecenter)**2/(2*self.sigma(ecenter)**2))#*Nq*np.exp(-(q_square)/(2*self.sigmaq**2))
        return Res
    
class PolynomialResNoQResolution(Resolution):
    """creates polynomial resolution subclass of Resolution Function with no Q resolutionp.  
    Takes a fourth order polynomial coefficient values provided either from user input 
    or taken from the Instrument Class."""
    def sigma(self,E):
        """ Energy resolution function for an instrument using coefficients of a characteristic polynomial.
         These values are typically supplied by an instrument scientist.
         """
        sig = (self.poly[0] + self.poly[1]*E +self.poly[2]*E**2 + self.poly[3]*E**3 +self.poly[4]*E**4)/2.35 #convert to unit to c
        return sig
    
    def Gauss(self,evalues,qvalues,ecenter,qcenter):
        Ne=1.0/np.sqrt(2*np.pi*self.sigma(ecenter)**2)
        #Nq=1.0/np.sqrt(2*np.pi*self.sigmaq**2)
        #delta_q=qvalues-qcenter
        #q_square=np.dot(delta_q,delta_q)
        Res= Ne*np.exp(-(evalues-ecenter)**2/(2*self.sigma(ecenter)**2)) #*Nq*np.exp(-(q_square)/(2*self.sigmaq**2))
        return Res



class InstrumentResolution(Resolution):
    #print('here Resolution')
    """creates constant resolution subclass of Resolution Functionp.  
    Takes a FWHM value provided either from user input 
    as ElasticFWHM or taken from the Instrument Class."""
    def sigmaconst(self):
        sig = self.inst.FWHM/2.355
        return sig
    
    def Gausselastic(self,evalues,ecenter,escale):
        Ne=1.0/np.sqrt(2*np.pi*(self.sigmaconst())**2)
        Res= Ne*np.exp(-0.5*(evalues-ecenter)**2/(self.sigmaconst()**2))#*Nq*np.exp(-0.5*(q_square)/((fwhmq*self.sigmaq)**2))
        #print(np.trapz(Res))
        return Res
    def GaussQelastic(self,qvalues,qcenter,qscale):
        #Nq=1.0/np.sqrt(2*qscale*np.pi*(self.sigmaq)**2)
        sigma = self.sigmaq*np.linalg.norm(qcenter)
        q_diff = qvalues - qcenter ## This is an array of vectors.
        q_square = np.hstack([np.linalg.norm(q)**2 for q in q_diff])
        Res= np.exp(-(0.5/qscale)*(q_square)/((sigma)**2))/np.sum(np.exp(-(0.5/qscale)*(q_square)/((sigma)**2)))
        #plt.plot(Res)
        #plt.savefig("Qres.png")
        #plt.close()
        #print(np.trapz(Res))
        return Res

    """creates instrument resolution subclass of Resolution Function with Q resolutionp.  
    Uses instrument parameters to estimate the energy resolution as a function of energy transfer
    using the Windsor approximationp."""

    def energyToVelocity(self,energy):
        epsilon = 1.0e-4
        #if energy is less than zero here, an invalid value has been passed inp.
        
        if energy < epsilon :
            energy = 1e-3
        
        velocity=630.2*np.sqrt(energy/2.072) #meter per second
        return velocity 

    def sigma(self,ecenter,escale):

        L1 = self.inst.L1
        L2 = self.inst.L2
        L3 = self.inst.L3
        fwhm_mon1 = self.inst.dtp
        fwhm_mon2 = self.inst.dtc
        tsamp = self.inst.dtd
        vi = self.energyToVelocity(self.inst.ei)
        ef = self.inst.ei - ecenter  #calculate the energy transfer
        vf = self.energyToVelocity(ef)
    
        dE = 1.5*escale*m_neut*np.sqrt(((vi**3/L1) + L2*vf**3/(L1*L3))**2*fwhm_mon1**2 +
                     ((vi**3/L1) + (L1+L2)*vf**3/(L1*L3))**2*fwhm_mon2**2 +
                     ((vf**3/L3)**2*(tsamp/vf)**2)) ## BANDER modified 1.5*fwhm*
        dE = dE* joulesTomeV
        return dE    
    def McVinesigma(self,ecenter,escale):
        #from https://rez.mcvine.ornl.gov/
        mcvinesigma = np.load('/path/to/Interpolated_arcs_res_ARCS-100-1.5-AST_420.0_Ei_50.0.npy')
        Sigma = mcvinesigma[:,1][np.abs(mcvinesigma[:,0]-ecenter).argmin()]
        return Sigma
    
    def Gauss(self,evalues,qvalues,ecenter,qcenter,escale,qscale):
        Ne=0.5/np.sqrt(2*np.pi*(self.sigma(ecenter,escale))**2)
        #Nq=1.0/np.sqrt(2*np.pi*self.sigmaq**2)
        #delta_q=qvalues-qcenter
        #q_square=np.dot(delta_q,delta_q)
        Res= Ne*np.exp(-(evalues-ecenter)**2/(2*self.sigma(ecenter,escale)**2))#*Nq*np.exp(-(q_square)/(2*self.sigmaq**2))

        return Res

    def GaussQ(self,qvalues,qcenter,qscale):
        q_diff = qvalues - qcenter ## This is an array of vectors.
        q_square = np.hstack([np.linalg.norm(q)**2 for q in q_diff])
        Res= np.exp(-(0.5/qscale)*(q_square)/((self.sigmaq)**2))/np.trapz(np.exp(-(0.5/qscale)*(q_square)/((self.sigmaq)**2)))
        #print(np.trapz(Res))
        return Res

    def Laurentz(self,evalues,ecenter,escale):
        gamma = self.McVinesigma(ecenter,escale)*2.2355
        Res =  1/np.pi * 0.5 * gamma / ((evalues-ecenter)**2+(0.5*gamma)**2)
        return (Res)/np.trapz(Res)
    def Voigt(self,evalues,ecenter,escale):
        gamma = self.McVinesigma(ecenter,escale)*2.2355
        Res = 1/np.pi * 0.5 * (gamma)*evalues / ((evalues**2-(ecenter)**2)**2+(0.5*(gamma)*evalues)**2)
        return Res/np.trapz(Res)
    def GaussMcVine(self,evalues,ecenter,escale):
        sigma = self.McVinesigma(ecenter,escale)
        Res= np.exp(-(evalues-ecenter)**2/(2*sigma**2))                                           
        return Res/np.trapz(Res)
    def PseudoVoigt(self,evalues,ecenter,escale):
        b = 0.5
        return b*self.GaussMcVine(evalues,ecenter,escale) + (1-b)*self.Laurentz(evalues,ecenter,escale)

class InstrumentResolutionNoQResolution(Resolution):
    """creates instrument resolution subclass of Resolution Function with no Q resolutionp.  
        Uses instrument parameters to estimate the energy resolution as a function of energy transfer
    using the Windsor approximationp."""
    
    def energyToVelocity(self,energy):
        energy=float(energy)
        epsilon = 1.0e-4
        #if energy is less than zero here, an invalid value has been passed in, make this something small but finite.
        if energy < epsilon :
            energy = 1e-3
        
        velocity=630.2*np.sqrt(energy/2.072) #meter per second
        return velocity 
    
    def sigma(self,ef) :

        L1 = self.inst.L1
        L2 = self.inst.L2
        L3 = self.inst.L3
        fwhm_mon1 = self.inst.dtc
        fwhm_mon2 = self.inst.dtp
        tsamp = self.inst.dtd
        vi = self.energyToVelocity(self.inst.ei)
        vf = self.energyToVelocity(ef)

    
        dE = m_neut*np.sqrt(((vi**3/L1) + L2*vf**3/(L1*L3))**2*fwhm_mon2**2 +
                     ((vi**3/L1) + (L1+L2)*vf**3/(L1*L3))**2*fwhm_mon1**2 +
                     ((vf**3/L3)**2*(tsamp/vf)**2)) 
        dE = dE* joulesTomeV
        return dE    
    
    def Gauss(self,evalues,qvalues,ecenter,qcenter):
        print("InstrumentResNoQRes")
        Ne=1.0/np.sqrt(2*np.pi*self.sigma(ecenter)**2)
       
        Res= Ne*np.exp(-(evalues-ecenter)**2/(2*self.sigma(ecenter)**2))
        return Res

