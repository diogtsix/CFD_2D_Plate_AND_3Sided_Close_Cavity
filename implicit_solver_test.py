import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog
from solver_explicit import Solver_explicit
from solver_implicit import Solver_implicit
from preprocessor import Preprocessor
from blasius_solution import BlasiusSolution
grid_y_size = 0.09
free_flow_velocity = 1
dx = 0.01
dy = 0.01


PreProcess = Preprocessor(dx=dx, dy=dy, grid_y_size=grid_y_size, free_flow_velocity=free_flow_velocity)
PreProcess.create_grid()
PreProcess.initialize_velocity_field()
Result = Solver_implicit(PreProcess.grid_nodes_x,
                        PreProcess.grid_nodes_y,
                        PreProcess.grid_u_velocity,
                        PreProcess.grid_v_velocity,
                        PreProcess.dx,
                        PreProcess.dy)
Result.solve()

Result.Blasius_delta()
Result.Blasius_delta1_and_delta2()
Result.Blasius_t_wall_and_Cf()

Blasius = BlasiusSolution(Result.grid_nodes_x.shape[1],
                          Result.dx,
                          Result.dy)
Blasius.blasius_exact_solution()



