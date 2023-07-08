import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate



class Nodes:
    def __init__(self,id=[0,0], coord=[0,0],vel=[0,0]):
        self.coord=coord
        self.id=id
        self.vel=vel
            

def create_nodes(Dx,Dy):
   
    j_num=round(Heigh/dy)
    i_num=round(Width/dx)
    a=0
    b=0
    Nodes_num=i_num*j_num
    Nod=[]
    
    for i in range(i_num+1):
         for j in range(j_num+1):
            v1=1
            
            Nodd=Nodes(id=[j,i],coord=[b,a],vel=[v1,0])
            a+=Dy
            a=round(a,5)
            
            if Nodd.coord[0]>=0 and Nodd.coord[0]<=10:
                if Nodd.id[0]==0:
                   Nodd.vel[0]=0
            
            
            Nod.append(Nodd)
         a=0
         b+=Dx
         b=round(b,5)
         
    return Nod 

def eq3(self,Dx,Dy,viscosity,AN,PN,KN,N):
    
    AA=(2*Dx*viscosity/(Dy**2))*(PN.vel[0]-2*N.vel[0]+KN.vel[0])/N.vel[0]
    BB=(Dx/Dy)*(N.vel[1]*(PN.vel[0]-KN.vel[0])/N.vel[0])
    CC=AN.vel[0]+AA-BB
    self.vel[0]=CC
    

def eq2(self,Dx,Dy,DN,AN,KN):
    
    AA=KN.vel[1]
    BB=(Dy/Dx)*(DN.vel[0]-AN.vel[0])
    self.vel[1]=AA-BB

def eq4(self,Dy,PN):
  
    self.vel[1]=0.25*PN.vel[1]
    

def main_loop(Nod, Dx, Dy):
    
    j_num=round(Heigh/Dy)
    i_num=round(Width/Dx)
    
    
    for i in range(2*j_num+1, len(Nod)-j_num-1):
        
        
        N=Nod[i]
        if N.id[0]>=1 and N.id[0]<=j_num-1:
             DN=Nod[i+j_num+1]
             AN=Nod[i-j_num-1]
             PN=Nod[i+1]
             KN=Nod[i-1]

             eq3(DN,Dx,Dy,viscosity,AN,PN,KN,N)
             eq2(PN,Dx,Dy,DN,AN,KN)
        
             if PN.id[0]==2:
                 eq4(N,Dy,PN)
    
           
def draw1(self):
    
    
    for kk in [2,3,4,5,6,7,8,9,10]: 
   
        ccx=[]
        ccy=[]
    
    
        for nod in self:
         
          if nod.coord[0]==kk:
                
             ccx.append(nod.vel[0])
             ccy.append(nod.coord[1])
             
      
            
                 
    
    
        plt.plot(ccx,ccy,label=f"x={kk}")
        
    
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}    
    plt.grid()
    plt.axis([0,1,0,Heigh])       
    plt.legend()
    plt.xlabel("Ταχύτητα u [m/s]",**csfont)
    plt.ylabel("Απόσταση Y [m]",**csfont)
    plt.title("Ρητή Μέθοδος : Κατανομη ταχυτήτων u, σε διαφορετικες θέσεις x",**csfont)
    plt.show()
    

def CompareBlasius(Nod,dx,dy,i_num,j_num):
    
    # ----------Delta ------------------
    Deltay=[]
    Deltax=[]
    dd=0
    ddd=1
    for nod in Nod:
        
        if nod.coord[0]>0.01*dx:
            if nod.vel[0]>=0.99 and dd==0:
                Deltay.append(nod.coord[1])
                Deltax.append(nod.coord[0])

                dd=1
                ddd=0
            elif nod.coord[1]==0 and ddd==0:
                dd=0
                ddd=1
            
        
   
    

    # ----------Delta  1  &    Delta 2------------------
    Delta1y=[]
    Delta1x=[]
    Delta2y=[]
    Delta2x=[]
    delta1=0
    delta2=0
    i=0
    for nod in Nod:
        
       if nod.coord[0]>0.01*dx:
           delta1=delta1+(1-nod.vel[0])*dy
        
           delta2=delta2+nod.vel[0]*(1-(nod.vel[0]))*dy
        
           if nod.id[0]==j_num:
               Delta1y.append(delta1)
               Delta1x.append(nod.coord[0])
            
               Delta2y.append(delta2)
               Delta2x.append(nod.coord[0])
               
            
               delta1=0
               delta2=0
            
    # ----------t-wall & Cf------------------        
    
    AirDynamic = 1.789*(10**(-5))   
    density=1.225 
    t=0     
    tw=[] 
    Cf=[]  
    XX=[]
    for nod in Nod:
        if nod.coord[0]>0.01*dx:
           if nod.id[0]==0:
              u0=nod.vel[0]
           if nod.id[0]==1:
              u1=nod.vel[0]
              t=1
           if t==1:
              tw.append(AirDynamic*((u1-u0)/dy))
              ll=2*AirDynamic*((u1-u0)/dy)
              ll=ll/density
              Cf.append(ll)
              XX.append(nod.coord[0])
              t=0
            
            
    # --------Blasius exact Solutions - Equations------------
    
    Rex=np.zeros((i_num,1))
    Bd=np.zeros((i_num,1))
    Bd1=np.zeros((i_num,1))
    Bd2=np.zeros((i_num,1))
    Btw=np.zeros((i_num,1))
    Bcf=np.zeros((i_num,1))
    
    
    for ii in range(1,i_num):
        Rex[ii,0]=(density*ii*dx)/AirDynamic
        Bd[ii,0]=5*ii*dx/(Rex[ii,0]**0.5)
        Bd1[ii,0]=1.72*ii*dx/(Rex[ii,0]**0.5)
        Bd2[ii,0]=0.664*ii*dx/(Rex[ii,0]**0.5)
        Bcf[ii,0]=0.664/(Rex[ii,0]**0.5)    
    
    
    
    #---------------plots----------------------
    
    
    # ------------------------DELTA plots----------------
    plt.figure(1)
    plt.plot(Deltax,Deltay)
    plt.plot(Deltax,Bd)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('δ [m]',**csfont)
    plt.title('Πάχος οριακού Στρώματος δ',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    # plt.show()
    
    # -----------------DELTA 1 plots----------------
    
    
    plt.figure(2)
    plt.plot(Delta1x,Delta1y)
    plt.plot(Delta1x,Bd1)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('δ1 [m]',**csfont)
    plt.title('Πάχος Μετατόπισης οριακού Στρώματος δ1',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    # plt.show()
    
    
    
    
    # -----------------DELTA 2 plots----------------
    
    
    plt.figure(3)
    plt.plot(Delta2x,Delta2y)
    plt.plot(Delta2x,Bd2)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('δ2 [m]',**csfont)
    plt.title('Πάχος Απώλειας Ορμής οριακού Στρώματος δ2',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    # plt.show()
    
    
    
    # -----------------Cf plots----------------
    
    
    plt.figure(4)
    plt.plot(XX,Cf)
    plt.plot(XX,Bcf)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('Cf',**csfont)
    plt.title('Συντελεστής Τριβής Cf',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    # plt.show()
    
    
    
        # ----------------ERRROR---------------
    
    E=0
    
    for k in range(len(Deltay)):
       
        E=E+Deltay[k]-Bd[k]
    
    Er=E/len(Deltay)
    
    print(Er)
    

if __name__=="__main__":
    
    
    
    viscosity=1.46*10**(-5) 
    dx=0.001
  
    dy=0.01

    Heigh=0.09
    Width=12
    Width=10
    j_num=round(Heigh/dy)
    i_num=round(Width/dx) 
    
    
   
    Nod=create_nodes(dx,dy)
    main_loop(Nod,dx,dy)
  
    draw1(Nod)
 
    
    CompareBlasius(Nod,dx,dy,i_num,j_num)
    
   
    
   
