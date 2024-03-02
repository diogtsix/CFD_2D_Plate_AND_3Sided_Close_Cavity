import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog
from solvers.Cavity_solver import cavity_solver
from preprocessor.Cavity_preprocessor import cavity_preprocessor
import time 


start=time.time()
pre = cavity_preprocessor(1, 0.01, 0.01, 0.0001, 1.2, 1, 1, 0.0001)
pre.initialize_velocity_field()

solve = cavity_solver(pre)
solve.ADI_method()

end = time.time()
print( end - start)