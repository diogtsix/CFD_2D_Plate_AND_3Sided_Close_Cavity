import numpy as np
# We use finite differences to discritize the Prandlt equations for 
# our 2D flow field
class Solver_explicit: 
    
    air_dynamic_viscosity = 1.789*(10**(-5))
    density = 1.225 #kg/m3
                                   
    def __init__(self,grid_nodes_x,grid_nodes_y,grid_u_velocity,grid_v_velocity,dx,dy):
        
        self.grid_nodes_x = grid_nodes_x
        self.grid_nodes_y = grid_nodes_y
        self.grid_u_velocity = grid_u_velocity
        self.grid_v_velocity = grid_v_velocity
        self.dx = dx
        self.dy = dy
    

        
    
    def u_i_plus_1_j(self,Dx,Dy,viscosity,AN,PN,KN,N_u,N_v):
            
        first_term = (2*Dx*viscosity/(Dy**2))*(PN-2*N_u+KN)/N_u
        second_term = (Dx/Dy)*(N_v*(PN-KN)/N_u)
        velocity = AN+first_term-second_term
            
        return velocity
        
    def v_i_j_plus_1(self,Dx,Dy,DN,AN,KN):
            
        first_term=KN.vel[1]
        second_term=(Dy/Dx)*(DN.vel[0]-AN.vel[0])
        velocity=first_term-second_term
            
        return velocity
        
    def v_i_1(self,Dy,PN):
        
        velocity = 0.5*Dy*PN.vel[1]
     
        return velocity
            
            
    
    def solve(self):
        shape = np.shape(self.grid_nodes_x)
        elements_in_x = shape[1]
        elements_in_y = shape[0]
        
        for i in range(1, elements_in_x + 1):
            for j in range(1, elements_in_y + 1):
            

                Node = (i,j)
                left_node = (i-1,j)
                right_node = (i+1,j)
                top_node = (i,j+1)
                bottom_node = (i,j-1) 
                
               
                
                first_term_vel = self.u_i_plus_1_j(self.dx,self.dy, 
                                                   self.air_dynamic_viscosity,
                                                   self.grid_u_velocity[left_node],
                                                   self.grid_u_velocity[top_node],
                                                   self.grid_u_velocity[bottom_node],
                                                   self.grid_u_velocity[Node],
                                                   self.grid_v_velocity[Node],
                                                    )
                
                  
        
a = np.zeros([3,4])
print(a) 
b = (2,2)
print("\n")
print(a[b])