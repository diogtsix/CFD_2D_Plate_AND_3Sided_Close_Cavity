import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog
from solvers.Cavity_solver import cavity_solver
from preprocessor.Cavity_preprocessor import cavity_preprocessor
from postprocessor.Cavity_postprocess import cavity_postprocess
import time 


# Function to run the simulation with parameters from the GUI
def run_simulation():
    wall_length = float(entry_wall_length.get())
    dx = float(entry_dx.get())
    dy = float(entry_dy.get())
    dt = float(entry_dt.get())
    relax_factor = float(entry_relax_factor.get())
    wall_vel = float(entry_wall_vel.get())
    reynolds = int(entry_reynolds.get())
    error = float(entry_error.get())
    
    pre = cavity_preprocessor(wall_length, dx, dy, dt, relax_factor, wall_vel, reynolds, error)
    pre.initialize_velocity_field()

    solve = cavity_solver(pre)
    solve.ADI_method()

    postprocessObj = cavity_postprocess(solve.preprocessor, solve)
    postprocessObj.errorConvergence()
    postprocessObj.vorticyContour()
    postprocessObj.streamFunctionContour()

    plt.show()

# GUI setup
root = tk.Tk()
root.title("CFD_Closed_Cavity Simulation Parameters")

# Default values for each parameter
default_values = ["1", "0.01", "0.01", "0.0001", "1", "1", "1", "0.0001"]

# Create and place entry widgets for each parameter
parameters = ["Wall Length :", "dx :", "dy :", "dt :", "Relaxation Factor :", "Wall Velocity :", "Reynolds :", "Error :"]
entries = []
for i, parameter in enumerate(parameters):
    tk.Label(root, text=parameter).grid(row=i, column=0)
    entry = tk.Entry(root)
    entry.grid(row=i, column=1)
    entries.append(entry)

for entry, default in zip(entries, default_values):
    entry.insert(0, default)
    
entry_wall_length, entry_dx, entry_dy, entry_dt, entry_relax_factor, entry_wall_vel, entry_reynolds, entry_error = entries

# Button to run the simulation
run_button = tk.Button(root, text="Run Simulation", command=run_simulation)
run_button.grid(row=len(parameters), column=0, columnspan=2)

root.mainloop()