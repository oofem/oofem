if (USE_CEMHYD)
    set (tm_cemhyd
        Materials/cemhyd/cemhydmat.C
        )
endif ()

set (tm_emodel
    EngineeringModels/nltransienttransportproblem.C
    EngineeringModels/nonstationarytransportproblem.C
    EngineeringModels/stationarytransportproblem.C
    EngineeringModels/transienttransportproblem.C
    )

set (tm_element
    Elements/transportelement.C
    Elements/quad1_ht.C
    Elements/qquad1_ht.C
    Elements/tr1_ht.C
    Elements/quadaxisym1_ht.C
    Elements/traxisym1_ht.C
    Elements/brick1_ht.C
    Elements/qbrick1_ht.C
    Elements/wedge_ht.C
    Elements/qwedge_ht.C
    Elements/tetrah1_ht.C
    Elements/tr1darcy.C
    Elements/LatticeElements/latticetransportelement.C
    Elements/LatticeElements/lattice2d_mt.C
    Elements/LatticeElements/lattice3d_mt.C
    )

set (tm_boundary
    BoundaryCondition/gravitypressure.C
    BoundaryCondition/transportgradientneumann.C
    BoundaryCondition/transportgradientdirichlet.C
    BoundaryCondition/transportgradientperiodic.C
    )

set (tm_cs
    simpletransportcrosssection.C
    )

set (tm_material
    Materials/anisomassmat.C
    Materials/bazantnajjarmat.C
    Materials/hemobaznajmat.C
    Materials/hemokunzelmat.C
    Materials/hemotkmat.C
    Materials/hydram.C
    Materials/hydratingconcretemat.C
    Materials/hydratinghemomat.C
    Materials/hydratingisoheatmat.C
    Materials/isoheatmat.C
    Materials/isolinmoisturemat.C
    Materials/isomoisturemat.C
    Materials/LatticeMaterials/latticetransmat.C
    Materials/nlisomoisturemat.C
    Materials/nonlinearmassmat.C
    Materials/rvestokesflow.C
    Materials/transportmaterial.C
    Materials/twophasemat.C
    )

set (tm
    ${tm_emodel}
    ${tm_element}
    ${tm_material}
    ${tm_cs}
    ${tm_new}
    ${tm_cemhyd}
    ${tm_hydram}
    ${tm_boundary}
    )

add_library (tm
    ${LIB_TYPE}
    ${tm}
    )
