# this consists of the class TurbulenceRANS, to compute RANS turbulence boundary conditions

import numpy as np

# available turbulence models: kOmegaSST, kEpsilon, SpalartAllmaras
class turbulenceRANS:
    def __init__(self,U=1.0,nu=1.0,rho=1.0,L=1.0):
        self.mesh = None
        self.nu = None
        self.k = None
        self.epsilon = None
        self.omega = None
        self.U = U
        self.rho = rho
        self.L = L
        self.Re = U*L/nu


        self.I = None # turbulence intensity
        # Constants
        self.Cmu = 0.09
        self.sigma_k = 1.0
        self.sigma_epsilon = 1.3
        self.C1 = 1.44
        self.C2 = 1.92
        self.kappa = 0.41
    
    def set_intensity(self, intensity):
        self.I = intensity

    def calc_intensity(self):
        self.I = 0.16*self.Re**(-1./8.)


    def calc_k(self):
        self.k = 1.5*(self.U*self.I)**2

    


