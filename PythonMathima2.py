import matplotlib.pyplot as plt 
import numpy as np
import time


TOTAL_COLLISIONS=0




class Ball:
    
    def __init__ (self, pos=[0, 0], velocity=[1, 0], mass=1, radius=1):
        self.pos = pos
        self.velocity = velocity
        self.mass=mass
        self.radius=radius
        
    def draw(self,ax):
        ax.plot(self.pos[0],self.pos[1], "bo")
        
    def update_state(self, board, dt):
        
        self.pos[0]=self.pos[0]+self.velocity[0]*dt
        self.pos[1]=self.pos[1]+self.velocity[1]*dt
        
        self.draw(board)
        

    def collide(self,ball):
        temp_vel_x=self.velocity[0]
        temp_vel_y=self.velocity[1]
        self.velocity[0]=ball.velocity[0]
        self.velocity[1]=ball.velocity[1]
        ball.velocity[0]=temp_vel_x
        ball.velocity[1]=temp_vel_y
        

        
def main_loop(board,ax2, balls, board_size=100):
    global TOTAL_COLLISIONS
    t = 0
    dt=0.1
    plot_time=[]
    plot_collisions=[]
    while (t<1000):
        
        t0=time.time()
        
        
        board.clear()
        board.grid()
        board.axis("equal")
        
        
        
        board.plot([0, board_size, board_size,0 ,0],[0,0,board_size, board_size,0],'r')
        
        
        
        for ball in balls:
            ball.update_state(board,dt)
        
        
        
        do_collisions(balls, board_size)
        
        ax2.plot(plot_time, plot_collisions, 'r')
       
        
        t+=dt
        plt.draw()
        
        plt.pause(0.1)
        
        dt=time.time()-t0
        plot_time.append(t)
        plot_collisions.append(TOTAL_COLLISIONS)

def check_collision_balls(balls1,balls2):
    dist_check=np.sqrt((balls1.pos[0]-balls2.pos[0])**2+(balls1.pos[1]-balls2.pos[1])**2) <=(balls1.radius+balls2.radius)
    
    if dist_check:
        balls1.collide(balls2)


def do_collisions(balls, board_size):
    
    global TOTAL_COLLISIONS
    
    for i in range(len(balls)):
        
        collision = check_collision_walls(balls[i], board_size)
        
        for j in range(i+1,len(balls)):
            check_collision_balls(balls[j],balls[i])
            
            if collision:
                TOTAL_COLLISIONS += 1
                
                
        
      
def check_collision_walls(ball, board_size):
    
    collision = False
    
    if (ball.pos[0]-ball.radius<=0) or (ball.pos[0]+ball.radius>=board_size):
        ball.velocity[0]=ball.velocity[0]*(-1)
        collision = True
        
    if (ball.pos[1]-ball.radius<=0) or (ball.pos[1]+ball.radius>=board_size):
        ball.velocity[1]=ball.velocity[1]*(-1)
        
        collision = True
    return collision 

   
def initialize_balls(number_of_balls=40,board_size=100,velocity=3,avg_mass=1, ball_radius=1):
    balls = []
    
    for i in range(number_of_balls):
        ball=Ball(pos=[np.random.random()*board_size, np.random.random()*board_size],
                  velocity=[np.random.random()*2*velocity-velocity, np.random.random()*2*velocity-velocity],
                  mass=np.random.random()*avg_mass,
                  radius=ball_radius)
        balls.append(ball)
    return balls
     
          
        
        
if __name__== "__main__" :
   
   fig, (ax, ax2)=plt.subplots(1,2)
   ax.grid()
   ax.axis("equal")
   
   balls = initialize_balls()
   
   main_loop(ax,ax2, balls,100)
   






    