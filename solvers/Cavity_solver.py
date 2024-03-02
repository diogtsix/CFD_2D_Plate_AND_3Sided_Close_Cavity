import numpy as np
from utilities.thomasAlgorithm import Thomas

class cavity_solver:
    
    def __init__(self, preprocessor):
        
        self.preprocessor = preprocessor
        self.diagonal1 = None
        self.lower_diagonal1 = None
        self.upper_diagonal1 = None
        self.diagonal2 = None
        self.lower_diagonal2 = None
        self.upper_diagonal2 = None
    
    
    def ADI_method(self):
        
        pre = self.preprocessor
        
        # Build const and const vectors
        
        B = 1 + pre.dt / (pre.reynolds * (pre.dx ** 2))
        E = 1 - pre.dt / (pre.reynolds * (pre.dy ** 2))
        BB = 1 + pre.dt / (pre.reynolds * (pre.dy ** 2))
        EE = 1 - pre.dt / (pre.reynolds * (pre.dx ** 2))
        
        self.diagonal1 = ((pre.dt / (pre.dx**2)) + 1) * np.ones((pre.nodx - 2))
        self.lower_diagonal1 = -(pre.dt / (2 * pre.dx**2)) * np.ones((pre.nodx - 2))
        self.upper_diagonal1 = -(pre.dt / (2 * pre.dx**2)) * np.ones((pre.nodx - 2))
        
        self.diagonal2 = ((pre.dt / (pre.dy**2)) + 1) * np.ones((pre.nody - 2))
        self.lower_diagonal2 = -(pre.dt / (2 * pre.dy**2)) * np.ones((pre.nody - 2))
        self.upper_diagonal2 = -(pre.dt / (2 * pre.dy**2)) * np.ones((pre.nody - 2))
        
        self.lower_diagonal1[0] = 0
        self.upper_diagonal1[-1] = 0
        self.lower_diagonal2[0] = 0
        self.upper_diagonal2[-1] = 0
        
        # Stream Function Clalculations : Solving the 2 equation system
        
        self.solveStreamFunction(pre, B, E, BB, EE)
        
        
    def solveStreamFunction(self, pre, B, E, BB, EE):
        
        # constant values
        ddx = pre.dx**2
        ddy = pre.dy**2
        dt = pre.dt
        dx = pre.dx
        dy = pre.dy
        # Relaxation technique for Zhta and Psi matrices
        pre.zhta1[:,:] = (1 - pre.relax_factor) * pre.zhta1[:,:] + pre.relax_factor * pre.zhta3[:,:]
        pre.psi1[:,:] = (1 - pre.relax_factor) * pre.psi1[:,:] + pre.relax_factor * pre.psi3[:,:]
        
        # First TimeStep
        for jj in range(1, pre.nody - 1):
            
            RHS_constants =  np.zeros((pre.nodx - 2))
            
            for ii in range(1, pre.nodx - 1):
                
                a1 = (dt / 2) * pre.zhta1[jj, ii]   
                a2 = (dt /(2*ddy)) * pre.psi1[jj+1,ii]
                a3 = (1-(dt/(dy**2))) * pre.psi1[jj,ii]
                a4 = (dt/(2*(dy**2))) * pre.psi1[jj-1,ii]
                              
                RHS_constants[ii - 1] = a1 + a2 + a3 + a4
                
                if ii ==1 : 
                    
                    RHS_constants[ii - 1] = RHS_constants[ii - 1] + (dt/(2*ddx)) * pre.psi2[jj, 0]
                elif ii == pre.nodx - 2:
                    RHS_constants[ii - 1] = RHS_constants[ii - 1] + (dt/(2*ddx)) * pre.psi2[jj, pre.nodx - 1]
        
            pre.psi2[jj, 1: pre.nodx - 1] = Thomas(self.lower_diagonal1, self.diagonal1, 
                                                    self.upper_diagonal1, RHS_constants)
    
