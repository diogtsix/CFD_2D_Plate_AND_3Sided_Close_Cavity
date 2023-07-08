from re import L
import matplotlib.pyplot as plt
import numpy as np




def main_loop2():
    
    # --------------------Basic Matrices - Boundary Conditions-----------
    s=(j_num+1,i_num+1)
    Cordx=np.zeros(s)
    Cordy=np.zeros(s)
    Velu=np.ones(s)
    Velv=np.zeros(s)
    x=0
    y=0
    for i in range(i_num+1):
      for j in range(j_num+1):
          Cordy[j,i]=y
          Cordx[j,i]=x
          
          if x>=1-dx and x<=11 and y==0:  
              Velu[j,i]=0  
 
          y+=dy
      y=0    
      x+=dx
 

    #-------------------Calculation of u(i+1) column----------------------
   
  
         
    ss=(j_num-1,1)
    K=np.zeros(ss)
   
    
    ss2=(j_num-2,1)
    KatoDiag=np.ones(ss2)
    
    Diag=np.ones(ss)
    PanoDiag=np.ones(ss2)
    s1=0

    for i in range(1, i_num):  
        for j in range(1, j_num):
            
            A2=Velv[j,i]*xy-2*xx
            B2=4*Velu[j,i]*(y2)+4*xx
            C2=-Velv[j,i]*xy-2*xx
            
            Diag[s1,0]=B2
   
            if j!=1:
                KatoDiag[s1-1,0]=C2  
            if j!=(j_num-1):   
                PanoDiag[s1,0]=A2

            D2=4*(y2)*Velu[j,i]-4*xx
 
            K1=D2*Velu[j,i]
            K2=A2*Velu[j+1,i]
            K3=C2*Velu[j-1,i]
  
            K[s1,0]=K1-K2-K3
            
            if j==1 :
                l=C2*Velu[j-1,i+1]
               
            elif j==j_num-1 :
                l=-A2*Velu[j+1,i+1]
           
            else: 
                l=0
     
            K[s1,0]=K[s1,0]+l 
            
            
            s1+=1  
        s1=0
        u=Thomas(KatoDiag,Diag,PanoDiag,K)
        u=np.transpose(u)
        Velu[1:j_num,i+1]=u

     #-------------------Calculation of v(i+1) column----------------------   
        for jj in range(1, j_num):
            Lamda=-(2*(dy/dx))*(Velu[jj,i+1]-Velu[jj,i])-(Velv[jj+1,i]-Velv[jj,i])
            Velv[jj+1,i+1]=Lamda+Velv[jj,i+1]

      #-------------------Plots------------------------------------------ 
       
    ccy=Cordy[:,10]      
    for kk in [20,30,40,50,60,70,80,90,100]:  
   
         
         ii=dx*kk    
         ccx=Velu[:,kk]
       
         plt.plot(ccx,ccy,label=f"x={ii}")
       

    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}    
    plt.grid()
    plt.axis([0,1,0,Heigh])       
    plt.legend()
    plt.xlabel("Ταχύτητα u [m/s]",**csfont)
    plt.ylabel("Απόσταση Y [m]",**csfont)
    plt.title("Πεπλεγμένη Μέθοδος: Κατανομη ταχυτήτων u, σε διαφορετικες θέσεις x",**csfont)
    plt.show()

    # Compare Blasius -  Numerical-----------------------
    
    
     # ----------Delta ------------------
     
     
    Deltay=np.zeros((i_num+1,1))
    Deltax=np.zeros((i_num+1,1))
   
  
    for i in range(0,i_num+1):
        for j in range(0,j_num+1):
            if Velu[j,i]>=0.99:
                Deltay[i,0]=Cordy[j,i]
                Deltax[i,0]=Cordx[j,i]
                break
    


    # ----------Delta  1  &    Delta 2------------------
    Delta1y=np.zeros((i_num+1,1))
    Delta1x=np.zeros((i_num+1,1))
    Delta2y=np.zeros((i_num+1,1))
    Delta2x=np.zeros((i_num+1,1))
  

    for i in range(0,i_num+1):
        
        Delta1x[i,0]=Cordx[j,i]
        Delta2x[i,0]=Cordx[j,i]
        for j in range(0,j_num+1):
            
            Delta1y[i,0]=Delta1y[i,0]+(1-Velu[j,i])*dy
            
            Delta2y[i,0]=Delta2y[i,0]+Velu[j,i]*(1-Velu[j,i])*dy
            
            
            
    # ----------t-wall & Cf------------------        
    
    AirDynamic = 1.789*(10**(-5))   
    density=1.225 
      
    tw=np.zeros((i_num+1,1))
    Cf=np.zeros((i_num+1,1))  
    XX=Cordx[0,:]
    for i in range(0,i_num+1):
            u0=Velu[0,i]
            u1=Velu[1,i]
            
            tw[i,0]=AirDynamic*(u1-u0)/dy
            
            Cf[i,0] =2*tw[i,0]/density
           
            
    # --------Blasius exact Solutions - Equations------------
    
    Rex=np.zeros((i_num,1))
    Bd=np.zeros((i_num,1))
    Bd1=np.zeros((i_num,1))
    Bd2=np.zeros((i_num,1))
    Btw=np.zeros((i_num,1))
    Bcf=np.zeros((i_num,1))
    XXX=np.zeros((i_num,1))
    
    for ii in range(1,i_num):
        Rex[ii,0]=(density*ii*dx)/AirDynamic
        Bd[ii,0]=5*ii*dx/(Rex[ii,0]**0.5)
        Bd1[ii,0]=1.72*ii*dx/(Rex[ii,0]**0.5)
        Bd2[ii,0]=0.664*ii*dx/(Rex[ii,0]**0.5)
        Bcf[ii,0]=0.664/(Rex[ii,0]**0.5)    
        XXX[ii,0]=ii*dx+1
    XXX[0,0]=XXX[0,0]+1
    


    
    #---------------plots----------------------
    
    
    # ------------------------DELTA plots----------------
    plt.figure(1)
    plt.plot(Deltax[9:i_num+1],Deltay[9:i_num+1])
    plt.plot(XXX,Bd)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('δ [m]',**csfont)
    plt.title('Πάχος οριακού Στρώματος δ',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    plt.show()
    
    # -----------------DELTA 1 plots----------------
    
    
    plt.figure(2)
    plt.plot(Delta1x[9:i_num+1],Delta1y[9:i_num+1])
    plt.plot(XXX,Bd1)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('δ1 [m]',**csfont)
    plt.title('Πάχος Μετατόπισης οριακού Στρώματος δ1',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    plt.show()
    
    
    
    
    # -----------------DELTA 2 plots----------------
    
    
    plt.figure(3)
    plt.plot(Delta2x[9:i_num+1],Delta2y[9:i_num+1])
    plt.plot(XXX,Bd2)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('δ2 [m]',**csfont)
    plt.title('Πάχος Απώλειας Ορμής οριακού Στρώματος δ2',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    plt.show()
    
    
    
    # -----------------Cf plots----------------
    
    
    plt.figure(4)
    plt.plot(XX[9:i_num+1],Cf[9:i_num+1])
    plt.plot(XXX,Bcf)
    csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
    plt.grid()
    plt.xlabel('Απόσταση x [m]',**csfont)
    plt.ylabel('Cf',**csfont)
    plt.title('Συντελεστής Τριβής Cf',**csfont)
    plt.legend(['Αριθμητική','Blasius'])   
    plt.show()
    
  
    
#-------------------Thomas Algorithm Function----------------------
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


if __name__=="__main__":
    
    
    AirDensity = 1.225                                                                      # Air density [kg/m^3]
    AirDynamic = 1.789*(10**(-5))                                                           # Dynamic viscosity [kg/ms]

    viscosity=AirDynamic/AirDensity
    n=viscosity
    dx=0.1
    dy=0.01      
    xx=n*dx
    xy=dx*dy
    y2=dy*dy
    Heigh=0.09
    
    Width=11
   
    j_num=round(Heigh/dy)
    i_num=round(Width/dx)

    main_loop2()
   