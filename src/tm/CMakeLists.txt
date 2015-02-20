if (USE_CEMHYD)
    set (tm_cemhyd
        cemhyd/cemhydmat.C
        )
endif ()

set (tm_emodel
    stationarytransportproblem.C
    nonstationarytransportproblem.C
    nltransienttransportproblem.C
    darcyflow.C
    transienttransportproblem.C
    )

set (tm_element
    transportelement.C
    quad1_ht.C
    tr1_ht.C
    quadaxisym1_ht.C
    traxisym1_ht.C
    brick1_ht.C
    qbrick1_ht.C
    tetrah1_ht.C
    tr1darcy.C
    latticetransportelement.C
    lattice2d_mt.C
    )

set (tm_boundary
    gravitypressure.C
    )

set (tm_cs
    simpletransportcrosssection.C
    )

set (tm_status
    transportmaterial.C
    latticetransmat.C
    )

set (tm_material
    hydram.C
    hydratingisoheatmat.C
    hydratinghemomat.C
    hydratingconcretemat.C
    hemotkmat.C
    hemokunzelmat.C
    hemobaznajmat.C
    isoheatmat.C
    isomoisturemat.C
    isolinmoisturemat.C
    bazantnajjarmat.C
    nlisomoisturemat.C
    nonlinearmassmat.C
    anisomassmat.C
    rvestokesflow.C
    )

set (tm
    ${tm_emodel}
    ${tm_element}
    ${tm_material}
    ${tm_cs}
    ${tm_status}
    ${tm_new}
    ${tm_cemhyd}
    ${tm_hydram}
    ${tm_boundary}
    )

add_library (tm
    ${LIB_TYPE}
    ${tm}
    )
