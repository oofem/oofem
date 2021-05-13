import oofempy
import numpy as np
import pyvista as pv
 
dr=oofempy.OOFEMTXTDataReader("concrete_3point.in")
problem=oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)
problem.init()
problem.solveYourself()
vtkxmlPy = oofempy.vtkmemory(1, problem, domain_all=True, tstep_all=True, dofman_all=True, element_all=True, vars=(1,4), primvars=(1,), cellvars = (47,), stype=1, pythonExport=1)
vtkxmlPy.initialize()
vtkxmlPy.doOutput(problem.giveCurrentStep(), False)

for p in vtkxmlPy.getVTKPieces():
#p = vtkxmlPy.getVTKPieces()[0]
    print ("Piece:", p)
    print(p.getVertices())
    print(p.getCellConnectivity())
    print(p.getCellTypes(vtkxmlPy))
    disp = p.getPrimaryVertexValues(oofempy.UnknownType.DisplacementVector)
    sig = p.getInternalVertexValues(oofempy.InternalStateType.IST_StressTensor)
    sigx = sig[:, 0]
    
    grid = pv.UnstructuredGrid(p.getCellConnectivity(), p.getCellTypes(vtkxmlPy), p.getVertices())
    grid.point_arrays['Sigma_xx'] = sigx
    grid['Disp'] = disp
    print(grid.active_vectors)
    warped = grid.warp_by_vector('Disp', factor=1000.)
    p = pv.Plotter()
    p.add_mesh(warped, scalars='Sigma_xx')
    p.add_mesh(grid, style='wireframe', color='black')
    p.set_viewup((0,1,0))
    p.show()
problem.terminateAnalysis()

 
