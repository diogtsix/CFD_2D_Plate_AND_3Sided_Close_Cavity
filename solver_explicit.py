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
            
        first_term=KN
        second_term=(Dy/Dx)*(DN-AN)
        velocity=first_term-second_term
            
        return velocity
        
    def v_i_1(self,Dy,PN):
        
        velocity = 0.5*Dy*PN
     
        return velocity
            
            
    
    def solve(self):
        shape = np.shape(self.grid_nodes_x)
        elements_in_x = shape[1]
        elements_in_y = shape[0]
        
        for i in range(1, elements_in_x-1):
            for j in range(1, elements_in_y-1):
            

                Node = (j,1)
                left_node = (j,i-1)
                right_node = (j,i+1)
                top_node = (j+1,i)
                bottom_node = (j-1,i) 
                
               
                
                first_term_vel = self.u_i_plus_1_j(self.dx,
                                                   self.dy, 
                                                   self.air_dynamic_viscosity,
                                                   self.grid_u_velocity[left_node],
                                                   self.grid_u_velocity[top_node],
                                                   self.grid_u_velocity[bottom_node],
                                                   self.grid_u_velocity[Node],
                                                   self.grid_v_velocity[Node],
                                                    )
                

                
                self.grid_u_velocity[right_node] = first_term_vel
                
                second_term_vel = self.v_i_j_plus_1(self.dx,
                                                    self.dy,
                                                    self.grid_u_velocity[right_node],
                                                    self.grid_u_velocity[left_node],
                                                    self.grid_v_velocity[bottom_node],
                                                    )
                                
                self.grid_v_velocity[top_node] = second_term_vel
                
                if self.grid_nodes_y[Node] == self.dy:
                    third_term_vel = self.v_i_1(self.dy,
                                                self.grid_v_velocity[top_node]
                                                )
                    
                    self.grid_v_velocity[Node] = third_term_vel