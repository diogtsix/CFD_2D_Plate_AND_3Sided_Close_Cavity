import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog
from solvers.Cavity_solver import cavity_solver
from preprocessor.Cavity_preprocessor import cavity_preprocessor

pre = cavity_preprocessor(1, 0.01, 0.01, 0.001, 1, 1, 1, 0.001)
pre.initialize_velocity_field()

solve = cavity_solver(pre)
solve.ADI_method()

