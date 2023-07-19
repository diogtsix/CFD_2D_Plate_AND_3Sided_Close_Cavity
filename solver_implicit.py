import numpy as np
# We use finite differences to discritize the Prandlt equations for 
# our 2D flow field
class Solver_implicit: 
    
    
    def __init__(self,grid_nodes_x,grid_nodes_y,grid_u_velocity,grid_v_velocity,dx,dy):
        
        self.grid_nodes_x = grid_nodes_x
        self.grid_nodes_y = grid_nodes_y
        self.grid_u_velocity = grid_u_velocity
        self.grid_v_velocity = grid_v_velocity
        self.dx = dx
        self.dy = dy
        
        
    def u_i_plus_1_system_eq(self):
        
        #Create trisdiagonal coefficient matrix for our system
        nodes_on_y_axis = np.round(self.grid_nodes_y.shape[0])
        nodes_on_x_axis = np.round(self.grid_nodes_x.shape[1])
        LHS_constants = np.zeros((nodes_on_y_axis - 1, 1))
        lower_diagonal = np.ones((nodes_on_y_axis - 2, 1))

        diagonal = np.ones((nodes_on_y_axis - 1,1))
        upper_diagonal = np.ones((nodes_on_y_axis - 2, 1))
        
        xx = self.dx**2
        yy = self.dy**2
        xy = self.dx*self.dy
        s1 =0
        for i in range(1,nodes_on_x_axis-1):
            for j in range(1,nodes_on_y_axis-1):
                A2 = self.grid_v_velocity[j,i]*xy-2*xx
                B2 = 4*self.grid_u_velocity[j,i]*yy+4*xx
                C2 = -self.grid_v_velocity[j,i]*xy-2*xx
                
                diagonal[s1,0] = B2
                
                if j!=1:
                    lower_diagonal[s1-1,0] = C2
                if j!=(nodes_on_y_axis-1):
                    upper_diagonal[s1,0] = A2
                    
                D2 = 4*yy*self.grid_v_velocity[j,i] - 4*xx
                
                K1 = D2*self.grid_u_velocity[j,i]
                K2 = A2*self.grid_u_velocity[j+1,i]
                K3 = C2*self.grid_u_velocity[j-1,i]
                
                LHS_constants[s1,0] = K1 - K2 - K3
                
                
                if j==1 :
                    l=C2*self.grid_u_velocity[j-1,i+1]
                elif j==nodes_on_y_axis-1 :
                    l=-A2*self.grid_u_velocity[j+1,i+1]
                else: 
                    l=0
     
                LHS_constants[s1,0]=LHS_constants[s1,0]+l 
            
            
                s1+=1  
        
            s1=0
            u=Solver_implicit.Thomas(self,lower_diagonal,diagonal,upper_diagonal,LHS_constants)
            u=np.transpose(u)
            self.grid_u_velocity[1:nodes_on_y_axis,i+1]=u
            
            Solver_implicit.v_i_plus_1(self,i,nodes_on_y_axis)
            
            
    def v_i_plus_1(self,i,nodes_on_y_axis):
        for jj in range(1,nodes_on_y_axis-1):
            Lamda = -(2*(self.dy/self.dx))*(self.grid_u_velocity[jj,i+1]-self.grid_u_velocity[jj,i])-(self.grid_v_velocity[jj+1,i]-self.grid_v_velocity[jj,i])
            self.grid_v_velocity[jj+1,i+1]=Lamda+self.grid_v_velocity[jj,i+1]
            
        return self.grid_v_velocity
            
    def solve(self):
        Solver_implicit.u_i_plus_1_system_eq(self)
        
                
    @staticmethod
    def Thomas(self, a, b, c, d):
        nf = len(d)  # number of equations
        ac, bc, cc, dc = map(np.array, (a, b, c, d))  # copy arrays
        for it in range(1, nf):
            mc = np.divide(ac[it - 1], bc[it - 1])

            bc[it] -= mc * cc[it - 1]
            dc[it] -= mc * dc[it - 1]

        xc = np.zeros(nf)
        xc[-1] = dc[-1] / bc[-1]

        for il in range(nf - 2, -1, -1):
            xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

        return xc