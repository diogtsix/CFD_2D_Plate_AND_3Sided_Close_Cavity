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
        
        self.diagonal3 = None
        self.lower_diagonal3 = None
        self.upper_diagonal3 = None
        
        self.diagonal4 = None
        self.lower_diagonal4 = None
        self.upper_diagonal4 = None
    
    
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
        
        self.diagonal3 = np.zeros((pre.nodx-2))
        self.lower_diagonal3 = np.zeros((pre.nodx-2))
        self.upper_diagonal3 = np.zeros((pre.nodx-2))
        
        self.diagonal4 = np.zeros((pre.nody-2))
        self.lower_diagonal4 = np.zeros((pre.nody-2))
        self.upper_diagonal4 = np.zeros((pre.nody-2))
        
        ITER = []
        errorMat = []
        error = 1
        iter = 1
        
        while error>(10**(-1)):
            # Stream Function Clalculations : Solving the 2 equation system
        
            pre = self.solveStreamFunction(pre)
        
            pre =  self.solveVorticityFunction(pre, B, E, BB, EE)
        
            pre = self.velocityFieldCalculation(pre)
        
            # Display error and iterations while running 
        
            iter, error =  self.algorithmCharacteristics(pre, ITER, error, errorMat, iter)
        

        
    def solveStreamFunction(self, pre):
        
        # constant values
        ddx = pre.dx**2
        ddy = pre.dy**2
        dt = pre.dt
        dx = pre.dx
        dy = pre.dy
        # Relaxation technique for Zhta and Psi matrices
        pre.zhta1[:,:] = (1 - pre.relax_factor) * pre.zhta1[:,:] + pre.relax_factor * pre.zhta3[:,:]
        pre.psi1[:,:] = (1 - pre.relax_factor) * pre.psi1[:,:] + pre.relax_factor * pre.psi3[:,:]
        
        # First TimeStep (half - timestep) ---------------------------------------------------------------------
        # Vectorized computation for RHS_constants

        a1 = (dt / 2) * pre.zhta1[1:pre.nody-1, 1:pre.nodx-1]
        a2 = (dt / (2*ddy)) * pre.psi1[2:pre.nody, 1:pre.nodx-1]
        a3 = (1 - (dt / (ddy))) * pre.psi1[1:pre.nody-1, 1:pre.nodx-1]
        a4 = (dt / (2*ddy)) * pre.psi1[:pre.nody-2, 1:pre.nodx-1]

        RHS_constants = a1 + a2 + a3 + a4

        # Handling boundary conditions
        RHS_constants[:, 0] += (dt / (2*ddx)) * pre.psi2[1:pre.nody-1, 0]
        RHS_constants[:, -1] += (dt / (2*ddx)) * pre.psi2[1:pre.nody-1, -1]    
        
        for jj in range(1, pre.nody - 1):     
            pre.psi2[jj, 1: pre.nodx - 1] = Thomas(self.lower_diagonal1, self.diagonal1, 
                                                    self.upper_diagonal1, RHS_constants[jj-1, :])
    

        # Second Time Step ( second half time step) ----------------------------------------
        
        a1 = (dt / 2) * pre.zhta1[1:pre.nody-1, 1:pre.nodx-1]
        a2 = (dt / (2*ddx)) * pre.psi2[1:pre.nody - 1, 2:pre.nodx]
        a3 = (1 - (dt / (ddx))) * pre.psi2[1:pre.nody-1, 1:pre.nodx-1]
        a4 = (dt / (2*ddx)) * pre.psi2[1:pre.nody-1, :pre.nodx-2]

        RHS_constants = a1 + a2 + a3 + a4       
        
        
        # Handling boundary conditions
        RHS_constants[0 , :] += (dt / (2*ddy)) * pre.psi3[0, 1:pre.nodx-1]
        RHS_constants[-1, :] += (dt / (2*ddy)) * pre.psi3[ - 1, 1:pre.nodx-1]    
        
        
        for jj in range(1, pre.nodx - 1):     
            pre.psi3[1 : pre.nody - 1, jj] = Thomas(self.lower_diagonal2, self.diagonal2, 
                                                    self.upper_diagonal2, RHS_constants[:, jj-1])
    
        return pre
    
    def solveVorticityFunction(self, pre, B, E, BB, EE):
        
        pre = self.vorticityBoundaryConditions(pre)
        ddx = pre.dx**2
        ddy = pre.dy**2
        dx = pre.dx
        dy = pre.dy
        dt = pre.dt
        
        # First Step ( first half time step)
        for jj in range(1, pre.nody - 1):
            
            RHS_constants = np.zeros((pre.nodx-2))
            lk = 1
            lp = 0
            
            for ii in range(1, pre.nodx - 1): 
                
                psy =(pre.psi3[jj + 1, ii] - pre.psi3[jj - 1, ii]) / (2 *dy)
                psx =(pre.psi3[jj, ii + 1] - pre.psi3[jj,ii - 1]) / (2 *dx)
    
                A = - (dt * psy) / (4*dx) - dt / (2*pre.reynolds*ddx)
                C = (dt * psy) / (4*dx) - dt/ (2 * pre.reynolds*ddx)
                D = dt/ (2 * pre.reynolds*ddy) - (dt * psx) / (4 *dy)
                F = (dt/(2 * pre.reynolds*ddy)) + (dt * psx)/ ( 4 * dy)
                
                if ii == 1:
                    
                    k = D*pre.zhta1[jj - 1, ii] + E* pre.zhta1[jj,ii] + F*pre.zhta1[jj+1, ii]
                    RHS_constants[ii-1] = k - A * pre.zhta2[jj, ii - 1]
                    
                    self.upper_diagonal3[lp] = C
                    self.diagonal3[ii - 1] = B
                    lp = lp + 1
                    
                elif ii == (pre.nodx - 2) : 
                    
                    k = D * pre.zhta1[jj - 1, ii] + E * pre.zhta1[jj, ii] + F * pre.zhta1[jj+1, ii]
                    RHS_constants[ii - 1] = k - C * pre.zhta2[jj, ii+1]
                    
                    self.lower_diagonal3[lk] = A
                    self.diagonal3[ii- 1] = B 
                    lk = lk + 1
                else :
                    
                    k = D * pre.zhta1[ jj - 1, ii] + E *pre.zhta1[jj,ii] + F *pre.zhta1[jj+1, ii]
                    RHS_constants[ii-1] = k
                    
                    self.upper_diagonal3[lp] = C
                    self.lower_diagonal3[lk] = A
                    self.diagonal3[ii - 1] = B
                    lp = lp + 1
                    lk = lk + 1
            
            pre.zhta2[jj, 1: pre.nodx - 1] = Thomas(self.lower_diagonal3, self.diagonal3, 
                                                     self.upper_diagonal3, RHS_constants)
            
        # Second step (second half time step)
        
        for ii in range(1, pre.nodx - 1) :
            
            RHS_constants = np.zeros((pre.nodx-2))
            lk = 1
            lp = 0
            
            for jj in range(1, pre.nody - 1):
                
        
                psy =(pre.psi3[jj + 1, ii] - pre.psi3[jj - 1, ii]) / (2 *dy)
                psx =(pre.psi3[jj, ii + 1] - pre.psi3[jj,ii - 1]) / (2 *dx)
    
                AA =  (dt * psx) / (4*dy) - dt / (2*pre.reynolds*ddy)
                CC = - (dt * psx) / (4 * dy) - dt/ (2 * pre.reynolds*ddy)
                DD = dt/ (2 * pre.reynolds*ddx) + (dt * psy) / (4 *dx)
                FF = (dt/(2 * pre.reynolds*ddx)) - (dt * psy)/ ( 4 * dx)
        
                if jj == 1:
                    
                    k = DD*pre.zhta2[jj, ii - 1] + EE* pre.zhta2[jj,ii] + FF*pre.zhta2[jj, ii+1]
                    RHS_constants[jj-1] = k - AA * pre.zhta3[jj - 1, ii]
                    
                    self.upper_diagonal4[lp] = CC
                    self.diagonal4[jj - 1] = BB
                    lp = lp + 1
                    
                elif jj == (pre.nody - 2) : 
                    
                    k = DD * pre.zhta2[jj, ii - 1] + EE * pre.zhta2[jj, ii] + FF * pre.zhta2[jj, ii + 1]
                    RHS_constants[jj - 1] = k - CC * pre.zhta3[jj+1, ii]
                    
                    self.lower_diagonal4[lk] = AA
                    self.diagonal4[jj- 1] = BB 
                    lk = lk + 1
                else :
                    
                    k = DD * pre.zhta2[ jj, ii - 1] + EE *pre.zhta2[jj,ii] + FF *pre.zhta2[jj, ii+1]
                    RHS_constants[jj-1] = k
                    
                    self.upper_diagonal4[lp] = CC
                    self.lower_diagonal4[lk] = AA
                    self.diagonal4[jj - 1] = BB
                    lp = lp + 1
                    lk = lk + 1       
                 
            pre.zhta3[1:pre.nody - 1, ii] = Thomas(self.lower_diagonal4, self.diagonal4, 
                                                     self.upper_diagonal4, RHS_constants)   
            
        return pre
                    
    def vorticityBoundaryConditions(self, pre): 
        
        ddx = pre.dx**2
        ddy = pre.dy**2
        
        # Left Wall BC
        val1 = (2/(ddx)) * (pre.psi3[:,0] - pre.psi3[:,1])
        pre.zhta1[:, 0] = val1
        pre.zhta2[:, 0] = val1
        pre.zhta3[:, 0] = val1
        
        # Right Wall BC
        val2 = (2/(ddx)) * (pre.psi3[:,pre.nodx - 1] - pre.psi3[:,pre.nodx - 2])
        pre.zhta1[:, pre.nodx - 1] = val2
        pre.zhta2[:, pre.nodx - 1] = val2
        pre.zhta3[:, pre.nodx - 1] = val2
        
        # Bottom Wall BC
        val3= (2/(ddy)) * (pre.psi3[0, :] - pre.psi3[1, :])
        pre.zhta1[0, :] = val3
        pre.zhta2[0, :] = val3
        pre.zhta3[0, :] = val3    
        
        # Top Wall BC
        val4= (2/(ddy)) * (pre.psi3[pre.nody-2,:] - pre.psi3[pre.nody-1,:] + pre.dy + pre.wall_vel)
        pre.zhta1[pre.nody-1,:] = val4
        pre.zhta2[pre.nody-1,:] = val4
        pre.zhta3[pre.nody-1,:] = val4   
        
        return pre        
          
          
    def velocityFieldCalculation(self, pre):
        
        # Calculate u_vel using slices for the entire array at once
        a1 = pre.psi3[2:pre.nody, 1:pre.nodx-1]
        a2 = pre.psi3[:pre.nody-2, 1:pre.nodx-1]
        pre.u_vel[1:pre.nody-1, 1:pre.nodx-1] = (a1 - a2) / (2 * pre.dy)
    
        # Calculate v_vel using slices for the entire array at once
        b1 = pre.psi3[1:pre.nody-1, 2:pre.nodx]
        b2 = pre.psi3[1:pre.nody-1, :pre.nodx-2]
        pre.v_vel[1:pre.nody-1, 1:pre.nodx-1] = (b1 - b2) / (2 * pre.dx)
        
        return pre

    def algorithmCharacteristics(self, pre, ITER, error, errorMat, iter):
        error1 = np.square(np.subtract(pre.zhta1, pre.zhta3))
        
        ITER.append(iter)
        print(error)
        error = error1.sum()
        errorMat. append(error)
        
        print(iter)
        iter += 1
        
        return iter, error
        
        