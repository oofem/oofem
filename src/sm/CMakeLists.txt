set (sm_boundary_conditions
    Deprecated/prescribedgenstrainshell7.C
        prescribeddispsliphomogenization.C
        prescribeddispslipbcdirichletrc.C
        prescribeddispslipbcneumannrc.C
        prescribeddispslipmultiple.C
        transversereinfconstraint.C
    )

set (sm_crosssection
    CrossSections/structuralcrosssection.C
    CrossSections/simplecrosssection.C
    CrossSections/layeredcrosssection.C
    CrossSections/fiberedcs.C
    CrossSections/structuralinterfacecrosssection.C
    CrossSections/variablecrosssection.C
    CrossSections/warpingcrosssection.C
    CrossSections/latticecrosssection.C
    )

set (sm_emodel
    EngineeringModels/structengngmodel.C
    EngineeringModels/staticstructural.C
    EngineeringModels/linearstatic.C
    EngineeringModels/nlinearstatic.C
    EngineeringModels/nlineardynamic.C
    EngineeringModels/eigenvaluedynamic.C
    EngineeringModels/linearstability.C
    EngineeringModels/adaptlinearstatic.C
    EngineeringModels/adaptnlinearstatic.C
    EngineeringModels/deidynamic.C
    EngineeringModels/diidynamic.C
    EngineeringModels/nldeidynamic.C
    EngineeringModels/incrementallinearstatic.C
    EngineeringModels/freewarping.C 
    EngineeringModels/structuralmaterialevaluator.C
    EngineeringModels/qclinearstatic.C
    EngineeringModels/pdelta.C
    )

set (sm_loads
    Loads/tf1.C
    Loads/usrdeftempfield.C
    Loads/foreigntempfieldload.C
    Loads/rotatingboundary.C
    Loads/constantpressureload.C
    Loads/structtemperatureload.C
    Loads/structeigenstrainload.C
    )

set (sm_error
    ErrorEstimators/huertaerrorestimator.C
    ErrorEstimators/zzerrorestimator.C
    ErrorEstimators/combinedzzsiee.C
    ErrorEstimators/scalarerrorindicator.C
    ErrorEstimators/directerrorindicatorrc.C
    )

set (sm_export
    export/dmexportmodule.C
    export/poiexportmodule.C
    export/gnuplotexportmodule.C
    )

set (sm_barrier
    polylinenonlocalbarrier.C
    symmetrybarrier.C
    )

set (sm_xfem
    xfem/xfemstructuralelementinterface.C
    xfem/enrichmentitems/crack.C
    xfem/enrichmentitems/delamination.C
    xfem/enrichmentitems/shellcrack.C
    xfem/propagationlaws/plhoopstresscirc.C
    xfem/propagationlaws/plprincipalstrain.C
    xfem/propagationlaws/plmaterialforce.C
    xfem/propagationlaws/plczdamageradius.C
    xfem/xfemstructuremanager.C
    xfem/matforceevaluator.C
    xfem/nucleationcriteria/ncprincipalstress.C
    xfem/nucleationcriteria/ncprincipalstrain.C
    EngineeringModels/xfemsolverinterface.C
    )

set (sm_unsorted
    particle.C
    refinedelement.C
    refinedmesh.C
    gpinitmodule.C
    stressstrainbasevector.C stressvector.C strainvector.C
    mappers/primvarmapper.C
    Contact/cdefnode2node.C
    Contact/celnode2node.C
    Contact/ActiveBc/node2nodepenaltycontact.C
    Contact/ActiveBc/node2nodelagrangianmultipliercontact.C
    )

set (sm_parallel
    FETISolver/feticommunicator.C
    FETISolver/fetiboundarydofman.C
    FETISolver/fetisolver.C
    )

set (sm_quasicontinuum
    Quasicontinuum/quasicontinuum.C
    Quasicontinuum/quasicontinuumnumberingscheme.C
    Quasicontinuum/quasicontinuumvtkxmlexportmodule.C    
    Quasicontinuum/fullsolveddomain.C
    Quasicontinuum/geometrygenerator.C
    )

set (sm
    ${sm_crosssection}
    ${sm_loads}
    ${sm_emodel}
    ${sm_export}
    ${sm_barrier}
    ${sm_error}
    ${sm_xfem}
    ${sm_unsorted}
    ${sm_new}
    ${sm_boundary_conditions}
    ${sm_quasicontinuum} 
    )

if (USE_PARALLEL)
    list (APPEND sm ${sm_parallel})
endif ()

add_library(sm
    ${LIB_TYPE}
    ${sm}
    )
add_dependencies(sm version)

# Not included:
#pnldeidynamiccomm
