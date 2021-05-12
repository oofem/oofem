
Postprocessing
##############

This section covers the basic postprocessing when using oofem python bindings.

One can always directly access results using oofem API after solving each solution step. 
Another option consist of using existing or developing custom export module. 
OOFEM comes with many built-in export modules to produce output in desired format.

In the following sections, we will cover some techniques in more detail.

VTK postprocessing
==================
This section illustrates how to query the data from oofem in python and use vtk library to visualize the results.

Prerequisites
--------------
Here we will use pyvista python module (https://pyvista-doc.readthedocs.io/en/latest/) to interface with Visualization toolkit (VTK).
The recommended procedure to install pyvista is to use python package installer:

.. code-block:: bash

    pip install pyvista
    
The approach uses built-in vtkmemory export module. This module allows in memory pythonic access to oofem grid data and variables.
First, we use existing model. The model is read from input file and solved:

.. code-block:: python3

    import oofempy
    import numpy as np
    import pyvista as pv
    
    dr=oofempy.OOFEMTXTDataReader("concrete_3point.in")
    problem=oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)
    problem.init()
    problem.solveYourself()

Next, we create instance of oofem vtkmemory export module and associate it to our problem. 
The constructor allows to filter the output to specific solution steps 
(here we use ``tstep_all=True`` to match all solution steps) and allows to select variables to export 
(we request dipslacement vector as primary variables, stress and strain tensors as internal 
variables and element number as cell variables, see oofem Input manual for details).

Further, the module is initialized and output is prepared. 

.. code-block:: python3

    vtkxmlPy = oofempy.vtkmemory(1, problem, domain_all=True, tstep_all=True, dofman_all=True, element_all=True, vars=(1,4), primvars=(1,), cellvars = (47,), stype=1, pythonExport=1)
    vtkxmlPy.initialize()
    vtkxmlPy.doOutput(problem.giveCurrentStep(), False)

The vtkmemory module allows in memory access to so called VTKPieces, objects representing piece (subset of mesh) to visualize.
VTKPieces contain all needed information and data about piece geometry, connectivity and about exported variables.
We start looping over the pieces, request data needed to instanciate pyvista UnstructuredGrid instance and 
attach exported variables to it.

.. code-block:: python3

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

You can follow https://github.com/oofem/oofem/blob/master/bindings/python/examples/vtkdemo.py for a full example.