import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog
from solver_explicit import Solver_explicit
from solver_implicit import Solver_implicit
from preprocessor import Preprocessor
from bl import 

grid_y_size = 0.9
free_flow_velocity = 1
dx = 0.1
dy = 0.01


PreProcess = Preprocessor(dx=dx, dy=dy, grid_y_size=grid_y_size, free_flow_velocity=free_flow_velocity)
PreProcess.create_grid()
PreProcess.initialize_velocity_field()
Result = Solver_explicit(PreProcess.grid_nodes_x,
                        PreProcess.grid_nodes_y,
                        PreProcess.grid_u_velocity,
                        PreProcess.grid_v_velocity,
                        PreProcess.dx,
                        PreProcess.dy)
Result.solve()

Result.Blasius_delta()
Result.Blasius_delta1_and_delta2()
Result.Blasius_t_wall_and_Cf()

Blasius = blasius_solution()

plt.figure(1)
plt.plot(Result.x_delta_position,Result.y_delta_position)


plt.figure(2)
plt.plot(Result.x_delta1_position,Result.y_delta1_position)

plt.figure(3)
plt.plot(Result.x_delta2_position,Result.y_delta2_position)

plt.figure(4)
plt.plot(Result.Cf_x_cordinate,Result.Cf_friction_coeff)

plt.show()

