import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog
from solver_explicit import Solver_explicit
from solver_implicit import Solver_implicit
from preprocessor import Preprocessor


class CFDApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("CFD Solver")
        self.root.geometry('300x300')

        # Create buttons
        self.explicit_button = tk.Button(self.root, text="Solver_explicit", command=self.run_explicit_solver, padx=10, pady=5)
        self.explicit_button.pack()

        self.implicit_button = tk.Button(self.root, text="Solver_implicit", command=self.run_implicit_solver, padx=10, pady=5)
        self.implicit_button.pack()

        self.compare_button = tk.Button(self.root, text="Compare", command=self.compare_solvers, padx=10, pady=5)
        self.compare_button.pack()

        self.root.mainloop()

    def run_explicit_solver(self):
        solver_type = 'explicit'
        dx, dy = self.get_step_sizes(solver_type)
        if dx and dy:
            self.run_solver(Solver_explicit, dx, dy)

    def run_implicit_solver(self):
        solver_type = 'implicit'
        dx, dy = self.get_step_sizes(solver_type)
        if dx and dy:
            self.run_solver(Solver_implicit, dx, dy)

    def run_solver(self, solver_class, dx, dy):
        grid_y_size = 0.1
        free_flow_velocity = 1

        PreProcess = Preprocessor(dx=dx, dy=dy, grid_y_size=grid_y_size, free_flow_velocity=free_flow_velocity)
        PreProcess.create_grid()
        PreProcess.initialize_velocity_field()

        Result = solver_class(PreProcess.grid_nodes_x,
                              PreProcess.grid_nodes_y,
                              PreProcess.grid_u_velocity,
                              PreProcess.grid_v_velocity,
                              PreProcess.dx,
                              PreProcess.dy)

        Result.solve()

    def compare_solvers(self):
        solver_type = 'compare'
        dx, dy = self.get_step_sizes(solver_type)
        if dx and dy:
            grid_y_size = 0.1
            free_flow_velocity = 1

        # Create Preprocessor instance for explicit solver
            explicit_PreProcess = Preprocessor(dx=dx, dy=dy, grid_y_size=grid_y_size, free_flow_velocity=free_flow_velocity)
            explicit_PreProcess.create_grid()
            explicit_PreProcess.initialize_velocity_field()

            explicit_solver = Solver_explicit(explicit_PreProcess.grid_nodes_x,
                                          explicit_PreProcess.grid_nodes_y,
                                          explicit_PreProcess.grid_u_velocity,
                                          explicit_PreProcess.grid_v_velocity,
                                          explicit_PreProcess.dx,
                                          explicit_PreProcess.dy)

            explicit_solver.solve()

        # Create Preprocessor instance for implicit solver
            implicit_PreProcess = Preprocessor(dx=dx, dy=dy, grid_y_size=grid_y_size, free_flow_velocity=free_flow_velocity)
            implicit_PreProcess.create_grid()
            implicit_PreProcess.initialize_velocity_field()

            implicit_solver = Solver_implicit(implicit_PreProcess.grid_nodes_x,
                                          implicit_PreProcess.grid_nodes_y,
                                          implicit_PreProcess.grid_u_velocity,
                                          implicit_PreProcess.grid_v_velocity,
                                          implicit_PreProcess.dx,
                                          implicit_PreProcess.dy)

            implicit_solver.solve()

    def get_step_sizes(self, solver_type):
        if solver_type == 'implicit':
            dx = simpledialog.askfloat("Step Size", "Enter step size in x direction:",initialvalue="0.1")
            dy = simpledialog.askfloat("Step Size", "Enter step size in y direction:",initialvalue="0.01")
        elif solver_type == 'explixit':  
            dx = simpledialog.askfloat("Step Size", "Enter step size in x direction:",initialvalue="0.001")
            dy = simpledialog.askfloat("Step Size", "Enter step size in y direction:",initialvalue="0.01")        
        else:      
            dx = simpledialog.askfloat("Step Size", "Enter step size in x direction:",initialvalue="0.1")
            dy = simpledialog.askfloat("Step Size", "Enter step size in y direction:",initialvalue="0.01")        
        
            
        return dx, dy


if __name__ == "__main__":
    CFDApp()