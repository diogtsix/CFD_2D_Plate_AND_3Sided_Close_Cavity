# ------------Diogenis Tsichlakis-----------------------------------
# -----------eksamhno Spoudon: 10o----------------------------------
# -----------AEM : 6216------------------------------------


import matplotlib.pyplot as plt
import numpy as np
import time 



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


#-------------Starting Values ----------

Re=1
Er=[]
ITER=[]
wmega=1
Uw=-1
L=1
dx=0.01
dy=0.01
dt=0.0001
ddx=dx**2
ddy=dy**2


nodx=round(L/dx)+1
nody=round(L/dy)+1


#---------STEP 1 : Grid, BoundaryCon, Starting values, Basic matrices

s=(nody,nodx)
X=np.linspace(0,L,nodx)
Y=np.linspace(0,L,nody)
uvel=np.zeros(s)
vvel=np.zeros(s)
Psi1=np.zeros(s)
Psi2=np.zeros(s)
Psi3=np.zeros(s)
Zhta1=np.zeros(s)
Zhta2=np.zeros(s)
Zhta3=np.zeros(s)

uvel[-1,:]=Uw
e=1
iter=0

#-------------------Constant Variables-----------

B=1+dt/(Re*(dx**2))
E=1-dt/(Re*(dy**2))
BB=1+dt/(Re*(dy**2))
EE=1-dt/(Re*(dx**2))

#------------------------------------------            
# s1=(nodx-2,1)
# s2=(nodx-3,1)
# ss1=(nody-2,1)
# ss2=(nody-3,1)

s2=(nodx-2)

ss2=(nody-2)


Diag1=((dt/ddx)+1)*np.ones(s2)
KDiag1=-(dt/(2*ddx))*np.ones(s2)
PDiag1=-(dt/(2*ddx))*np.ones(s2)


Diag2=((dt/ddy)+1)*np.ones(ss2)
KDiag2=-(dt/(2*ddy))*np.ones(ss2)
PDiag2=-(dt/(2*ddy))*np.ones(ss2)

KDiag1[0]=0
PDiag1[-1]=0
KDiag2[0]=0
PDiag2[-1]=0

               #---------------First Time Step for Psi---------------------
start=time.time()
             
while e>(10**(-4)):
       

    
    Zhta1[:,:]=(1-wmega)*Zhta1[:,:]+wmega*Zhta3[:,:]
    Psi1[:,:]=(1-wmega)*Psi1[:,:]+wmega*Psi3[:,:]

    for j in range(1,nody-1):
        L=np.zeros(s2)
      
        
        
        for i in range(1,nodx-1):
            
            L[i-1]=(dt/2)*Zhta1[j,i]+(dt/(2*(dy**2)))*Psi1[j+1,i]+(1-(dt/(dy**2)))*Psi1[j,i] + (dt/(2*(dy**2)))*Psi1[j-1,i]

            
            if i==1:
                
                
                L[i-1]=L[i-1]+(dt/(2*ddx))*Psi2[j,0]
                
            elif i==nodx-2:
                
                L[i-1]=L[i-1]+(dt/(2*ddx))*Psi2[j,nodx-1]
           
       
         
      
        Psi2[j,1:nodx-1]=Thomas(KDiag1,Diag1,PDiag1,L)
        
    for i in range(1,nodx-1):
        
        L=np.zeros(ss2)
    
        for j in range(1,nody-1):
            
            L[j-1]=(dt/2)*Zhta1[j,i]+(dt/(2*ddx))*Psi2[j,i+1]+(1-(dt/ddx))*Psi2[j,i]+(dt/(2*ddx))*Psi2[j,i-1]
                
                
            if j==1:
                
                L[j-1]= L[j-1]+(dt/(2*ddy))*Psi3[0,i]
            elif j==nody-2:
                L[j-1]= L[j-1]+(dt/(2*ddy))*Psi3[j+1,i]

       
        
        Psi3[1:nody-1,i]=(Thomas(KDiag2,Diag2,PDiag2,L))   
 
  
#---------STEP 3 : Zhta Boundary Conditions-------------------------


        #----------left--------
        
    M1=(2/(ddx))*(Psi3[:,0]-Psi3[:,1])  
    # M1=(1-wmega)*Zhta3[:,0]+wmega*M1
    Zhta1[:,0]=M1 
    Zhta2[:,0]=M1   
    Zhta3[:,0]=M1 
        
        #-------------Right----------
        
    M2=(2/(ddx))*(Psi3[:,nodx-1]-Psi3[:,nodx-2])
    # M2=(1-wmega)*Zhta3[:,nodx-1]+wmega*M2
    Zhta1[:,nodx-1]=M2
    Zhta2[:,nodx-1]=M2
    Zhta3[:,nodx-1]=M2
        
        #---------------bottom   
    M3=(2/ddy)*(Psi3[0,:]-Psi3[1,:])
    # M3=(1-wmega)*Zhta3[0,:]+wmega*M3
    Zhta1[0,:]=M3
    Zhta2[0,:]=M3
    Zhta3[0,:]=M3 
        #---------------top (first row)
        
    M4=-(2/(ddy))*(Psi3[nody-2,:]-Psi3[nody-1,:]+Uw*dy)
    # M4=(1-wmega)*Zhta3[nody-1,:]+wmega*M4
    Zhta1[nody-1,:]=M4
    Zhta2[nody-1,:]=M4
    Zhta3[nody-1,:]=M4    
    #--------Zhta Calculation-----------------

    Diag3=np.zeros(s2)
    KDiag3=np.zeros(s2)
    PDiag3=np.zeros(s2)
    
    for j in range(1,nody-1):
        
        K=np.zeros(s2)
        lk=1
        lp=0
        
        for i in range(1,nodx-1):
            psy=(Psi3[j+1,i]-Psi3[j-1,i])/(2*dy)
            psx=(Psi3[j,i+1]-Psi3[j,i-1])/(2*dx)
            
            A=-(dt*psy)/(4*dx)-dt/(2*Re*(ddx))
            C=(dt*psy)/(4*dx)-dt/(2*Re*(ddx))
            
            D=dt/(2*Re*(ddy))-(dt*psx)/(4*dy)
            
            F=(dt/(2*Re*(ddy)))+(dt*psx)/(4*dy)
            
            if i==1:
                k=D*Zhta1[j-1,i]+E*Zhta1[j,i]+F*Zhta1[j+1,i]
                K[i-1]=k-A*Zhta2[j,i-1]
                
                PDiag3[lp]=C
                Diag3[i-1]=B
                lp=lp+1
                
            elif i==(nodx-2):
                k=D*Zhta1[j-1,i]+E*Zhta1[j,i]+F*Zhta1[j+1,i]
                K[i-1]=k-C*Zhta2[j,i+1]
                
                KDiag3[lk]=A
                Diag3[i-1]=B
                lk=lk+1
                
            else : 
                
                k=D*Zhta1[j-1,i]+E*Zhta1[j,i]+F*Zhta1[j+1,i]
                K[i-1]=k
                
                PDiag3[lp]=C
                KDiag3[lk]=A
                Diag3[i-1]=B
                lp=lp+1
                lk=lk+1
    
        
        Zhta2[j,1:nodx-1]=(Thomas(KDiag3,Diag3,PDiag3,K))    

    Diag4=np.zeros(ss2)
    KDiag4=np.zeros(ss2)
    PDiag4=np.zeros(ss2)
 
    for i in range(1,nodx-1):
        
        
        K=np.zeros(ss2)
        lk=1
        lp=0
        
        for j in range(1,nody-1):
            
            
            psy=(Psi3[j+1,i]-Psi3[j-1,i])/(2*dy)
            psx=(Psi3[j,i+1]-Psi3[j,i-1])/(2*dx)
            
            AA=(dt*psx)/(4*dy)-dt/(2*Re*(ddy))
            
            CC=-(dt*psx)/(4*dy)-dt/(2*Re*(ddy))
            DD=dt/(2*Re*(ddx))+(dt*psy)/(4*dx)
            
            FF=dt/(2*Re*(ddx))-(dt*psy)/(4*dx)
            
            if j==1:
                
                k=DD*Zhta2[j,i-1]+EE*Zhta2[j,i]+FF*Zhta2[j,i+1]
                K[j-1]=k-AA*Zhta3[j-1,i]
                
                PDiag4[lp]=CC
                Diag4[j-1]=BB
                lp=lp+1
                
                
            elif j==(nody-2):
                
                k=DD*Zhta2[j,i-1]+EE*Zhta2[j,i]+FF*Zhta2[j,i+1]
                K[j-1]=k-CC*Zhta3[j+1,i]
                
                KDiag4[lk]=AA
                Diag4[j-1]=BB
                lk=lk+1
            
            else:
                
                k=DD*Zhta2[j,i-1]+EE*Zhta2[j,i]+FF*Zhta2[j,i+1]
                K[j-1]=k
                
                KDiag4[lk]=AA
                PDiag4[lp]=CC
                Diag4[j-1]=BB
                lk=lk+1
                lp=lp+1
        
        
        Zhta3[1:nody-1,i]=(Thomas(KDiag4,Diag4,PDiag4,K)) 
            
    #------------------Velocities-----------
    
    for j in range(1,nody-1):
        
        for i in range(1,nodx-1):
          #---------STEP 7 : Velocities u,v Calculation-----------------------     
            uvel[j,i]=(Psi3[j+1,i]-Psi3[j-1,i])/(2*dy)
    
            vvel[j,i]=-(Psi3[j,i+1]-Psi3[j,i-1])/(2*dx)
        
    Error1=np.subtract(Zhta1,Zhta3)
    Error1=np.square(Error1)
    ITER.append(iter)
    print(e)
    e = Error1.sum()
    Er.append(e)
    iter+=1
    print(iter)
end=time.time()
TIME=end-start
print(TIME)
    

x=X
y=Y
   
[X,Y]=np.meshgrid(x,y)

# Re=1--------------
# lvl1=np.linspace(-200,53,400)
# lvl2=np.linspace(-0.1,0,30)

#Re=100---------------

# lvl2=np.linspace(-0.07,0,30)
# lvl1=np.linspace(-400,62,400)

#Re=500-----------------
# lvl2=np.linspace(-0.0609,0,30)
# lvl1=np.linspace(-200,97,400)

# plt.figure(1)

# Xromata = plt.cm.get_cmap("rainbow").copy()
# AA1=plt.contour(X,Y,Zhta3,lvl1,cmap=Xromata)
# AA1.cmap.set_over('blue')
# AA1.cmap.set_under('red')
# AA1.changed()
# plt.colorbar()
# plt.clabel(AA1,lvl1,colors='black')
# plt.xlabel('X distance')
# plt.ylabel('Y distance')
# plt.title(f'Vorticity Contours [Re = {Re}]')
# csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}

# plt.figure(2)

       
# Xromata = plt.cm.get_cmap("rainbow").copy()
# AA2=plt.contour(X,Y,Psi3,lvl2,cmap=Xromata)
# AA2.cmap.set_over('blue')
# AA2.cmap.set_under('red')
# AA2.changed()
# plt.colorbar()
# plt.clabel(AA2,lvl2,colors='black')
# plt.xlabel('X distance')
# plt.ylabel('Y distance')
# plt.title(f'Stream Function Contours [Re ={Re}]')
# csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}


# plt.figure(3)

# plt.plot(ITER[10:ITER[-1]+1],Er[10:ITER[-1]+1])
# plt.grid()

# plt.xlabel('Iteration Number')
# plt.ylabel('Error')
# plt.title(f' Error Convergence [Re ={Re} , ω = {wmega}]')
# csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'} 