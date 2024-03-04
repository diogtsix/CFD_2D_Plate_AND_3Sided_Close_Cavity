# CFD_2D_Plate
This is a project where we implement different CFD algorithms to calculate the flow field above a flat plate. 


start=time.time()
pre = cavity_preprocessor(1, 0.01, 0.01, 0.0001, 1, 1, 500, 0.0001)
pre.initialize_velocity_field()

solve = cavity_solver(pre)
solve.ADI_method()

end = time.time()
print( end - start)

postprocessObj = cavity_postprocess(solve.preprocessor, solve)

fig1 = postprocessObj.errorConvergence()
fig2 = postprocessObj.vorticyContour()
fig3 = postprocessObj.streamFunctionContour()

plt.show()