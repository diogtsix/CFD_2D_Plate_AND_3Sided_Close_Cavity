import matplotlib.pyplot as plt
import numpy as np

class cavity_postprocess:
    def __init__(self, preprocessorObj, solverObj):
        
        self.preprocessorObj = preprocessorObj
        self.solverObj = solverObj
        
        self.lvl1 = [] 
        self.lvl2 = []
        
        self.setLevels()


    def setLevels(self):
        
        if self.preprocessorObj.reynolds == 1:
            
            self.lvl1=np.linspace(-200,53,400)
            self.lvl2=np.linspace(-0.1,0,30)
            
        if self.preprocessorObj.reynolds == 100:
            self.lvl2=np.linspace(-0.07,0,30)
            self.lvl1=np.linspace(-400,62,400)   
                
        if self.preprocessorObj.reynolds == 500:
            self.lvl2=np.linspace(-0.0609,0,30)
            self.lvl1=np.linspace(-200,97,400)



    def streamFunctionContour(self):
        
        pre = self.preprocessorObj
        solver = self.solverObj
        plt.figure(1)
        
        Colors = plt.cm.get_cmap("rainbow").copy()
        AA2=plt.contour(pre.xCoords, pre.yCoords, pre.psi3, self.lvl2, cmap=Colors)
        AA2.cmap.set_over('blue')
        AA2.cmap.set_under('red')
        AA2.changed()
        plt.colorbar()
        plt.clabel(AA2, self.lvl2, colors='black')
        plt.xlabel('X distance')
        plt.ylabel('Y distance')
        plt.title(f'Stream Function Contours [Re ={pre.reynolds}]')
        csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
        
       # plt.show()
        
    def vorticyContour(self):
        
        pre = self.preprocessorObj
        solver = self.solverObj
        plt.figure(2)
        Color = plt.cm.get_cmap("rainbow").copy()
        AA1=plt.contour(pre.xCoords, pre.yCoords,pre.zhta3,self.lvl1,cmap=Color)
        AA1.cmap.set_over('blue')
        AA1.cmap.set_under('red')
        AA1.changed()
        plt.colorbar()
        plt.clabel(AA1, self.lvl1, colors='black')
        plt.xlabel('X distance')
        plt.ylabel('Y distance')
        plt.title(f'Vorticity Contours [Re = {pre.reynolds}]')
        csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'}
        
     #   plt.show()



    def errorConvergence(self):
        pre = self.preprocessorObj
        solver = self.solverObj
        
        plt.figure(3)

        plt.plot(solver.ITER[10:solver.ITER[-1]+1],solver.errorMat[10:solver.ITER[-1]+1])
        plt.grid()

        plt.xlabel('Iteration Number')
        plt.ylabel('Error')
        plt.title(f' Error Convergence [Re ={pre.reynolds} , ω = {pre.relax_factor}]')
        csfont = {'fontname':'Times New Roman','fontsize':15,'fontweight':'bold'} 
        
        #plt.show()