import numpy as np
# We use finite differences to discritize the Prandlt equations for 
# our 2D flow field
class Solver_explicit: 
    
    air_dynamic_viscosity = 1.789*(10**(-5))
    density = 1.225 #kg/m3
                                   
    def __init__(self,grid_nodes_x,grid_nodes_y,
                 grid_u_velocity,grid_v_velocity,dx,dy):
        
        self.grid_nodes_x = grid_nodes_x
        self.grid_nodes_y = grid_nodes_y
        self.grid_u_velocity = grid_u_velocity
        self.grid_v_velocity = grid_v_velocity
        self.dx = dx
        self.dy = dy
        self.x_delta_position = []
        self.y_delta_position = []
        self.x_delta1_position = []
        self.y_delta1_position = []
        self.x_delta2_position = []
        self.y_delta2_position = []
        self.t_wall = [] 
        self.Cf_friction_coeff = []
        self.Cf_x_cordinate = []
    

        
    
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
            

                Node = (j,i)
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
    
    def Blasius_delta(self):
        dd = 0
        ddd = 1
        for i in range(1,self.grid_u_velocity.shape[1]):
            for j in range(0,self.grid_u_velocity.shape[0]):
                if self.grid_nodes_x[j,i] < 10: # this is vald only for plate_length = 10
                    if self.grid_u_velocity[j,i] >= 0.99 and dd==0:
                        self.x_delta_position.append(self.grid_nodes_x[j,i])
                        self.y_delta_position.append(self.grid_nodes_y[j,i])
                    
                        dd = 1
                        ddd = 0
                    if self.grid_nodes_y[j,i] == 0 and ddd == 0:
                        dd = 0
                        ddd = 1    
        return self.x_delta_position , self.y_delta_position
                    
    def Blasius_delta1_and_delta2(self):
        delta1 = 0
        delta2 = 0
        dd =0        
        for i in range(1,self.grid_u_velocity.shape[1]):
            for j in range(0,self.grid_u_velocity.shape[0]):
                
                if self.grid_nodes_x[j,i] > 0.01*self.dx and self.grid_nodes_x[j,i]< 10:
                    
                    delta1 = delta1 + (1 - self.grid_u_velocity[j,i])*self.dy
                
                    delta2 = delta2 +self.grid_u_velocity[j,i]*(1 - self.grid_u_velocity[j,i])*self.dy
                    if j == self.grid_u_velocity.shape[0]-1 and dd == 0:
                    
                        self.x_delta1_position.append(self.grid_nodes_x[j,i])
                        self.y_delta1_position.append(delta1)
                        self.x_delta2_position.append(self.grid_nodes_x[j,i])
                        self.y_delta2_position.append(delta2)
                        dd = 1
                        delta1 = 0
                        delta2 = 0
                
                    if dd == 1 and self.grid_nodes_y[j,i] ==0:
                        dd = 0
                
    def Blasius_t_wall_and_Cf(self):
        t = 0 
        for i in range(1,self.grid_u_velocity.shape[1]):
            for j in range(0,self.grid_v_velocity.shape[0]):
                
                if self.grid_nodes_y[j,i] == 0:
                    u0 = self.grid_u_velocity[j,i]
                if self.grid_nodes_y[j,i] == self.dy:
                    u1 = self.grid_u_velocity[j,i]
                    t = 1
                if t == 1:
                    self.t_wall.append(self.air_dynamic_viscosity*((u1-u0)/self.dy))
                    ll = 2*self.air_dynamic_viscosity*((u1-u0)/self.dy)
                    ll = ll/self.density
                    self.Cf_friction_coeff.append(ll)
                    self.Cf_x_cordinate.append(self.grid_nodes_x[j,i])
    

                
    