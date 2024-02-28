import numpy as np


class BlasiusSolution:
    
    air_dynamic_viscosity = 1.789*(10**(-5))
    density = 1.225 #kg/m3
    
    def __init__(self, number_of_columms,dx,dy):
        
        self.number_of_columms = number_of_columms
        self.dx = dx
        self.dy = dy        
        self.Reynolds = np.zeros((number_of_columms,1))
        self.delta = np.zeros((number_of_columms,1))
        self.delta1 = np.zeros((number_of_columms,1))
        self.delta2 = np.zeros((number_of_columms,1))
        self.cf_friction_coeff = np.zeros((number_of_columms,1))
        self.x_cordinate = np.zeros((number_of_columms,1))
        



    def blasius_exact_solution(self):
        
        for ii in range(1,self.number_of_columms):
            self.Reynolds[ii,0] = (self.density*ii*self.dx)/self.air_dynamic_viscosity
            self.delta[ii,0] = 5*ii*self.dx/(self.Reynolds[ii,0]**0.5)
            self.delta1[ii,0] = 1.72*ii*self.dx/(self.Reynolds[ii,0]**0.5)
            self.delta2[ii,0] = 0.664*ii*self.dx/(self.Reynolds[ii,0]**0.5)
            self.cf_friction_coeff[ii,0] = 0.664/(self.Reynolds[ii,0]**0.5)
            self.x_cordinate[ii,0] = ii*self.dx+1
        self.x_cordinate[0,0] = self.x_cordinate[0,0]+1
            
        

        