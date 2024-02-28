import numpy as np

class Solver_implicit: 
  
    air_dynamic_viscosity = 1.789*(10**(-5))
    density = 1.225 #kg/m3
      
    def __init__(self, grid_nodes_x, grid_nodes_y, grid_u_velocity, grid_v_velocity, dx, dy):
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
        

    def u_i_plus_1_system_eq(self):
        nodes_on_y_axis = round(self.grid_nodes_y[-1,0]/self.dy)
        nodes_on_x_axis =round(self.grid_nodes_x[0,-1]/self.dx)

        # Fix the size of LHS_constants to (nodes_on_x_axis - 2, 1)
        LHS_constants = np.zeros((nodes_on_y_axis - 2, 1))
        lower_diagonal = np.ones((nodes_on_y_axis - 3, 1))
        diagonal = np.ones((nodes_on_y_axis - 2, 1))
        upper_diagonal = np.ones((nodes_on_y_axis - 3, 1))


        xx = self.dx**2
        yy = self.dy**2
        xy = self.dx * self.dy
        s1 = 0

        # Correct the range of iteration here
        for i in range(1, nodes_on_x_axis - 2):
            for j in range(1, nodes_on_y_axis - 2):
                A2 = self.grid_v_velocity[j, i] * xy - 2 * xx
                B2 = 4 * self.grid_u_velocity[j, i] * yy + 4 * xx
                C2 = -self.grid_v_velocity[j, i] * xy - 2 * xx

                diagonal[s1, 0] = B2

                if j != 1:
                    lower_diagonal[s1 - 1, 0] = C2
                if j != (nodes_on_y_axis - 1):
                    upper_diagonal[s1, 0] = A2

                D2 = 4 * yy * self.grid_v_velocity[j, i] - 4 * xx

                K1 = D2 * self.grid_u_velocity[j, i]
                K2 = A2 * self.grid_u_velocity[j + 1, i]
                K3 = C2 * self.grid_u_velocity[j - 1, i]

                LHS_constants[s1, 0] = K1 - K2 - K3

                if j == 1:
                    l = C2 * self.grid_u_velocity[j - 1, i + 1]
                elif j == nodes_on_y_axis - 1:
                    l = -A2 * self.grid_u_velocity[j + 1, i + 1]
                else:
                    l = 0
     
                LHS_constants[s1, 0] = LHS_constants[s1, 0] + l 
            
                s1 += 1   
            s1 = 0
        
            # Apply Thomas algorithm to solve the tridiagonal system
            u = Solver_implicit.Thomas(lower_diagonal, diagonal, upper_diagonal, LHS_constants)
            self.grid_u_velocity[1:nodes_on_y_axis - 1,i+1:i+2] = u
        
        # Update v-velocity component for next time step
            self.v_i_plus_1(nodes_on_x_axis)
            
    def v_i_plus_1(self, nodes_on_x_axis):
        for jj in range(1, self.grid_nodes_y.shape[0] - 1):
            Lamda = -(2 * (self.dy / self.dx)) * (self.grid_u_velocity[jj, nodes_on_x_axis - 1] - self.grid_u_velocity[jj, nodes_on_x_axis - 2]) - (self.grid_v_velocity[jj + 1, nodes_on_x_axis - 1] - self.grid_v_velocity[jj, nodes_on_x_axis - 1])
            self.grid_v_velocity[jj + 1, nodes_on_x_axis - 1] = Lamda + self.grid_v_velocity[jj, nodes_on_x_axis - 1]
    
    def solve(self):
        self.u_i_plus_1_system_eq()
                
    @staticmethod
    def Thomas(a, b, c, d):
    
        nf = len(d) # number of equations
        ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
        for it in range(1, nf):
  
            mc = ac[it-1]/bc[it-1]

            bc[it] = bc[it] - mc*cc[it-1] 
            dc[it] = dc[it] - mc*dc[it-1]
 
        xc = bc
        xc[-1] = dc[-1]/bc[-1]

        for il in range(nf-2, -1, -1):
            xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

        return xc
    
    
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

                
