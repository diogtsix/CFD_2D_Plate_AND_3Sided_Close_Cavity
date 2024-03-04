# CFD_2D_Plate_AND_3Sided_Close_Cavity
This is repository contains 2 projects that deal with CFD methods on real problems. The first one is the CFD_2D_Plate and the second the CFD_CLOSE_CAVITY
For both projects extensive theoritical and result analysis has been conducted and can be found in the documentation folder (3 pdf documents).

# CFD 2D Plate Project

  # Description:
  
This project is dedicated to solving the fluid dynamics around a flat plate using both explicit and implicit computational fluid dynamics (CFD) methods. It aims to analyze the boundary layer formation and characteristics of fluid flow in two dimensions, providing insights into friction coefficients, velocity profiles, and boundary layer thicknesses.

  # Components Used:

Python: For all computational modeling and solution algorithms.
NumPy: For efficient numerical calculations and matrix operations.
Matplotlib: For visualizing the results, including velocity fields and boundary layer profiles.
Tkinter: For creating a simple GUI to run simulations and visualize results interactively.

  # How to Use:

Clone the repository to your local machine.
Ensure Python and the required libraries (NumPy, Matplotlib, Tkinter) are installed.
Navigate to the project directory and run Plate_main.py to launch the application.
Use the GUI to select the solver type (Explicit/Implicit), set grid resolution (dx, dy), and start the simulation.
Visualize the results directly through the GUI or explore the output data files for further analysis.
Note: Detailed documentation and example usage are included in the repository to help get you started.


#  3Sided  - Sided  Lid Driven Close Cavity Project 
  # Project Description
This repository hosts a comprehensive Computational Fluid Dynamics (CFD) project focusing on simulating laminar flow within a two-dimensional rectangular cavity. It is developed using Python and embraces Object-Oriented Programming (OOP) principles to model the flow of incompressible fluids. The project applies the Alternating Direction Implicit (ADI) method and the Thomas algorithm for numerically solving partial differential equations, offering insights into fluid behavior, vorticity patterns, and the convergence of solutions.

  # Components
- Preprocessor: Initializes the simulation grid, setting up boundary conditions and mesh details.
- Solver: Implements the ADI method and the Thomas algorithm for solving the flow and vorticity equations iteratively.
- Postprocessor: Provides visualization of the simulation, including stream function contours, vorticity contours, and error convergence plots.
- 
  # How to Use
Setup: Clone the repository and ensure Python 3.x is installed along with necessary libraries (numpy, matplotlib, and tkinter).
Configuration: Adjust simulation parameters (e.g., Reynolds number, grid size) in the provided configuration file or directly within the simulation scripts.
Run Simulation: Execute the Cavity_main.py simulation script. A GUI will prompt you to enter or confirm simulation parameters.
View Results: Observe the real-time generation of flow and vorticity patterns. Analyze the error convergence to assess the simulation's stability.

  # Contributing
Contributions to enhance the simulation's accuracy, efficiency, or visualization capabilities are welcome. Please submit pull requests or open issues to discuss proposed changes.