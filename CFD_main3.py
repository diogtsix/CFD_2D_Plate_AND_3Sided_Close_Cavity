from cmath import e

import matplotlib.pyplot as plt
import numpy as np

from math import e




class Nodes:
    def __init__(self,id=[0,0], coord=[0,0],vel=[0,0], Tcoord=[0,0]):
        self.coord=coord
        self.id=id
        self.vel=vel
        self.Tcoord=Tcoord
            

def create_nodes(Dx,Dy):
   
    j_num=round(Heigh/dy)
    i_num=round(Width/dx)
    a=0
    b=0
    Nodes_num=i_num*j_num
    Nod=[]
    aa=0
    AA=[]
    AAA=[]
    A4=[]
   
   
    aa=np.linspace(0,(0.25*Heigh)**(1/2),j_num+1)
    
    
   
    
   
    
    for i in range(i_num+1):
         for j in range(j_num+1):
            v1=1
            
            aaa=4*aa[j]**2
            aaa=np.round(aaa,8)
            # print(aaa)
            
          
            
         
            
         
            Nodd=Nodes(id=[j,i],coord=[b,aaa],vel=[v1,0],Tcoord=[b,aa[j]])
      
            if Nodd.coord[0]>=-1 and Nodd.coord[0]<=11:
                if Nodd.id[0]==0:
                   Nodd.vel[0]=0
            
            # AA.append(aaa)
            # AAA.append(b)
            # A4.append(b+0.1)
            # Nod.append(Nodd)
            aaa=0
         a=0
         
         b+=Dx
         b=round(b,5)
   
       
        #  plt.plot(AAA,AA,'o--',label=f"x-y Plane")
        #  plt.plot(A4,aa,'o--',label=f"ξ-η Plane")
         
        #  plt.xlabel('x Node Distance (m)')
        #  plt.ylabel('y Node Distance (m)')
        #  plt.legend()
        #  plt.title('Physical and Computational Plane Discretization')
        #  plt.show()
         
         

    return Nod 
   

def eq3(self,Dx,Dy,viscosity,AN,PN,KN,N,Dj,Dh):
    
    A1=N.vel[1]/N.vel[0]

    A2=(8*N.Tcoord[1])**(-1)
    
   
  
    A3=D
    A4=PN.vel[0]-KN.vel[0]
    
    A=A1*A2*A3*A4
    
  
  
    B1=1/N.vel[0]
 
    B2=(8*N.Tcoord[1])**(-3)
    B3=Dv
    B4=PN.vel[0]-KN.vel[0]
    B=B1*B2*B3*B4
 
    C1=1/N.vel[0]

    C2=(8*N.Tcoord[1])**(-2)

    C3=Dv2
   
    C4=PN.vel[0]-2*N.vel[0]+KN.vel[0]
    C=C1*C2*C3*C4
    
 
    CC=N.vel[0]-A-B+C
 
    
  
    
    
    self.vel[0]=CC

    

def eq2(self,Dx,Dy,DN,AN,KN,Dj,Dh,N):
    
    AA1=8*(N.Tcoord[1])
    AA2=D
    AA3=DN.vel[0]-AN.vel[0]
    AA4=AA1*AA2*AA3
  
    
    self.vel[1]=KN.vel[1]-AA4
    

def eq4(self,Dy,PN,Dj,Dh):
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
        

             eq3(DN,Dx,Dy,viscosity,AN,PN,KN,N,Dj,Dh)
          
             eq2(PN,Dx,Dy,DN,AN,KN,Dj,Dh,N)

        
             if PN.id[0]==2:
                 eq4(N,Dy,PN,Dj,Dh)
    
    
    
           
def draw1(self):
    
    
    for kk in [2,3,4,5,6,7,8,9,10]: 
    # for kk in [1.928, 3]: 
   
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
    
    YY=[]
    aa=np.linspace(0,(0.25*Heigh)**(1/2),j_num+1)
    for j in aa:
        YY.append(4*j**2)
  
    # print(YY)
    # input()
    
    # ----------Delta ------------------
    Deltay=[]
    Deltax=[]
    dd=0
    ddd=1
    for nod in Nod:
        
        if nod.coord[0]>0.01*dx:
            if nod.vel[0]>=0.998 and dd==0:
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
           

        
            DY=8*nod.Tcoord[1]*Dh
          
            
          
            
               
            delta1=delta1+(1-nod.vel[0])*DY
            delta2=delta2+nod.vel[0]*(1-(nod.vel[0]))*DY

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
        DY=8*nod.Tcoord[1]*Dh
        if (nod.coord[0]>0.01*dx):
           if nod.id[0]==0:
              u0=nod.vel[0]
           if nod.id[0]==1:
              u1=nod.vel[0]
              t=1
           if t==1:
              
              if nod.coord[0]>=1.3 and nod.coord[0]<=11:
                   
                  tw.append(AirDynamic*((u1-u0)/DY))
                  ll=2*AirDynamic*((u1-u0)/DY)
                  ll=ll/density
                  Cf.append(ll)
                  XX.append(nod.coord[0])
                  t=0
              else:
                  tw.append(0)
                 
                  Cf.append(0)
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
    
    
    
    viscosity=1.46*(10**(-5) )

    # dx=0.1
    dx=0.001
  
    dy=0.01
    # dy=0.01

    Heigh=0.09
    Width=10
    
    
    j_num=round(Heigh/dy)
    i_num=round(Width/dx) 
    
    Nod=create_nodes(dx,dy)
    N1=Nod[0]
    N2=Nod[1]
    
    Dj=dx
    Dh=N2.Tcoord[1]-N1.Tcoord[1]
    
    D=Dj/(2*Dh)
    Dv=viscosity*D
    D2=Dj/(Dh**2)
    Dv2=viscosity*D2
  
   
    
    
    
   
    
    main_loop(Nod,dx,dy)
    
    # for n in Nod:
    #     print(n.coord, n.vel)
    
    # for n in Nod:
    #     print(n.coord, n.Tcoord)
   
  
  
  
    draw1(Nod)
    
    CompareBlasius(Nod,dx,dy,i_num,j_num)