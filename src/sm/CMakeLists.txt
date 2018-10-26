set (sm_element
    Elements/structuralelementevaluator.C
    Elements/igaelements.C
    Elements/structuralelement.C
    Elements/nlstructuralelement.C
    Elements/structural2delement.C
    Elements/structural3delement.C
    Elements/3D/space3delementevaluator.C
    Elements/3D/lspace.C 
    Elements/3D/lspacebb.C
    Elements/3D/qspace.C 
    Elements/3D/qspacegrad.C
    Elements/3D/q27space.C
    Elements/3D/ltrspace.C 
    Elements/3D/qtrspace.C 
    Elements/3D/qtrspacegrad.C
    Elements/3D/lwedge.C 
    Elements/3D/qwedge.C 
    Elements/3D/macrolspace.C
    Elements/3D/qspacegrad.C
    Elements/3D/qwedgegrad.C
    Elements/Bars/truss1d.C 
    Elements/Bars/truss2d.C 
    Elements/Bars/truss3d.C
    Elements/Bars/truss3dnl.C
    Elements/Bars/qtruss1d.C 
    Elements/Bars/qtruss1dgrad.C
    Elements/Beams/libeam2d.C 
    Elements/Beams/libeam2dnl.C
    Elements/Beams/libeam3d.C 
    Elements/Beams/libeam3dnl.C 
    Elements/Beams/libeam3dnl2.C 
    Elements/Beams/libeam3d2.C
    Elements/Beams/beambaseelement.C
    Elements/Beams/beam2d.C 
    Elements/Beams/beam3d.C
    Elements/Plates/cct.C
    Elements/Plates/dkt.C 
    Elements/Plates/qdkt.C
    Elements/Plates/quad1mindlin.C
    Elements/Plates/dkt3d.C
    Elements/Shells/rershell.C 
    Elements/Shells/tr_shell01.C 
    Elements/Shells/tr_shell02.C 
    Elements/Shells/cct3d.C
    Elements/Shells/quad1mindlinshell3d.C
    Elements/Shells/shell7base.C
    Elements/Shells/tr2shell7.C
    Elements/Shells/shell7basexfem.C
    Elements/Shells/tr2shell7xfem.C
    Elements/Shells/mitc4.C
    Elements/Axisymmetry/q4axisymm.C
    Elements/Axisymmetry/axisymm3d.C 
    Elements/Axisymmetry/l4axisymm.C
    Elements/PlaneStrain/quad1planestrain.C 
    Elements/PlaneStrain/qplanestrain.C 
    Elements/PlaneStrain/qplanestraingrad.C
    Elements/PlaneStrain/qtrplanestrain.C 
    Elements/PlaneStrain/qtrplanestraingrad.C
    Elements/PlaneStrain/trplanestrain.C 
    Elements/PlaneStress/planestresselementevaluator.C
    Elements/PlaneStress/trplanrot.C
    Elements/PlaneStress/trplanstrss.C 
    Elements/PlaneStress/qplanstrss.C 
    Elements/PlaneStress/q9planstrss.C 
    Elements/PlaneStress/qtrplstr.C 
    Elements/PlaneStress/planstrss.C
    Elements/PlaneStress/planstrssxfem.C
    Elements/PlaneStress/trplanstrssxfem.C
    Elements/PlaneStress/qtrplanstrssxfem.C
    Elements/PlaneStress/qtrplstrgrad.C
    Elements/PlaneStress/trplanrot3d.C
    Elements/PlaneStress/qplanestressgrad.C
    Elements/PlaneStress/trplanestressrotallman.C
    Elements/PlaneStress/basicelement.C
    Elements/PlaneStress/trplanestressrotallman3d.C
    Elements/PlaneStress/linquad3d_planestress.C
    Elements/AbaqusUserElement.C
    Elements/htselement.C
    Elements/latticestructuralelement.C
    Elements/graddpelement.C 
    Elements/lumpedmasselement.C
    Elements/springelement.C
    Elements/lattice2d.C
    Elements/tet21ghostsolid.C
    Elements/quad1platesubsoil.C 
    Elements/quad2platesubsoil.C
    Elements/tria1platesubsoil.C
    Elements/linedistributedspring.C
    Elements/trwarp.C
    Elements/PlaneStress/basicelement.C
    Elements/nodalspringelement.C
    Elements/structural2delement.C
    Elements/structural3delement.C
    #Elements/phasefieldelement.C
    #qplanstrssphf.C
    )

set (sm_interface_elements
    Elements/Interfaces/cohsur3d.C
    Elements/Interfaces/structuralinterfaceelement.C
    Elements/Interfaces/intelpoint.C
    Elements/Interfaces/intelline1.C
    Elements/Interfaces/intelline1intpen.C
    Elements/Interfaces/intelline2.C
    Elements/Interfaces/intelline2intpen.C
    Elements/Interfaces/intelsurftr1.C
    Elements/Interfaces/Deprecated/interfaceelement1d.C
    Elements/Interfaces/Deprecated/interfaceelem2dlin.C
    Elements/Interfaces/Deprecated/interfaceelem2dquad.C
    Elements/Interfaces/Deprecated/interfaceelem3dtrlin.C
    )

set (sm_boundary_conditions
    Deprecated/prescribedgenstrainshell7.C
    )

set (sm_crosssection
    CrossSections/structuralcrosssection.C
    CrossSections/simplecrosssection.C
    CrossSections/layeredcrosssection.C
    CrossSections/fiberedcs.C
    CrossSections/structuralinterfacecrosssection.C
    CrossSections/variablecrosssection.C
    CrossSections/warpingcrosssection.C
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

set (sm_material
    Materials/abaqususermaterial.C
    Materials/structuralmaterial.C
    Materials/linearelasticmaterial.C
    Materials/isolinearelasticmaterial.C
    Materials/anisolinearelasticmaterial.C
    Materials/ortholinearelasticmaterial.C
    Materials/BoneMaterials/trabbonematerial.C 
    Materials/BoneMaterials/trabbonenl3d.C 
    Materials/BoneMaterials/trabbonenl.C 
    Materials/BoneMaterials/trabbone3d.C 
    Materials/BoneMaterials/trabboneembed.C 
    Materials/BoneMaterials/trabbonenlembed.C 
    Materials/BoneMaterials/trabbonegrad3d.C
    Materials/ConcreteMaterials/concrete2.C 
    Materials/ConcreteMaterials/concrete3.C 
    Materials/ConcreteMaterials/concretedpm.C 
    Materials/ConcreteMaterials/concretedpm2.C
    Materials/ConcreteMaterials/idm1.C 
    Materials/ConcreteMaterials/idmnl1.C 
    Materials/ConcreteMaterials/idmgrad1.C
    Materials/ConcreteMaterials/mazarsmodel.C 
    Materials/ConcreteMaterials/mazarsmodelnl.C
    Materials/ConcreteMaterials/mplasticmaterial2.C
    Materials/ConcreteMaterials/compodamagemat.C
    Materials/ConcreteMaterials/concretefcm.C
    Materials/ConcreteMaterials/concretefcmviscoelastic.C
    Materials/ConcreteMaterials/frcfcm.C
    Materials/ConcreteMaterials/frcfcmnl.C
    Materials/RheoChainMaterials/rheoChM.C 
    Materials/RheoChainMaterials/maxwellChM.C 
    Materials/RheoChainMaterials/kelvinChM.C
    Materials/RheoChainMaterials/kelvinChSolM.C
    Materials/RheoChainMaterials/cebfip78.C
    Materials/RheoChainMaterials/doublepowerlaw.C
    Materials/RheoChainMaterials/b3mat.C
    Materials/RheoChainMaterials/b3solidmat.C
    Materials/RheoChainMaterials/mps.C
    Materials/RheoChainMaterials/mpsdammat.C
    Materials/RheoChainMaterials/eurocode2creep.C
    Materials/structuralfe2material.C
    Materials/rcsd.C 
    Materials/rcsde.C 
    Materials/rcsdnl.C
    Materials/mooneyrivlin.C
    Materials/m4.C
    Materials/microplane_m1.C
    Materials/druckerPragerPlasticitySM.C
    Materials/druckerpragercutmat.C
    Materials/j2mat.C
    Materials/mdm.C
    Materials/micromaterial.C
    Materials/hyperelasticmaterial.C
    Materials/misesmat.C 
    Materials/misesmatgrad.C 
    Materials/misesmatnl.C
    Materials/rankinemat.C
    Materials/rankinematgrad.C
    Materials/rankinematnl.C
    Materials/latticematstatus.C
    Materials/latticedamage2d.C
    Materials/dustmat.C
    Materials/lsmastermat.C
    Materials/lsmastermatgrad.C
    Materials/simplevitrificationmaterial.C
    Materials/anisodamagemodel.C
    Materials/winklerpasternak.C
    Materials/winklermodel.C
    Materials/structmatsettable.C
    Materials/tutorialmaterial.C
    Materials/graddpmaterialextensioninterface.C
    Materials/rcm2.C
    Materials/isodamagemodel.C
    Materials/steelrelaxmat.C
    Materials/fcm.C
    )

if (USE_PYTHON)
    list (APPEND sm_material Materials/structuralpythonmaterial.C)
endif ()

set (sm_interface_materials
    Materials/InterfaceMaterials/structuralinterfacematerial.C
    Materials/InterfaceMaterials/bondceb.C
    Materials/InterfaceMaterials/cohint.C
    Materials/InterfaceMaterials/intmatelastic.C
    Materials/InterfaceMaterials/intmatbilinczfagerstrom.C
    Materials/InterfaceMaterials/intmatbilinczjansson.C
    Materials/InterfaceMaterials/intmatbilinearcz.C
    Materials/InterfaceMaterials/intmatbilinczelastic.C
    Materials/InterfaceMaterials/intmatbilinczfagerstromrate.C
    Materials/InterfaceMaterials/intmatisodamage.C
    Materials/InterfaceMaterials/intmatisodamagetable.C
    Materials/InterfaceMaterials/intmatcoulombcontact.C
    Materials/InterfaceMaterials/intmatdummycz.C    
    Materials/InterfaceMaterials/Deprecated/isointerfacedamage01.C
    Materials/InterfaceMaterials/Deprecated/isointerfacedamage02.C
    Materials/InterfaceMaterials/Deprecated/simpleinterfacemat.C
    Materials/InterfaceMaterials/masonry02.C
    Materials/InterfaceMaterials/mat_cebfip90.C
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

set (sm_status
    Materials/structuralms.C
    Materials/InterfaceMaterials/structuralinterfacematerialstatus.C
    Materials/deformationtheorymaterial.C
    Materials/microplanematerial.C 
    Materials/microplanematerial_bazant.C
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
    )

set (sm_parallel
    FETISolver/feticommunicator.C
    FETISolver/fetiboundarydofman.C
    FETISolver/fetisolver.C
    )

set (sm_obsolete
    Materials/Obsolete/steel1.C
    Materials/Obsolete/plasticmaterial.C
    Materials/Obsolete/perfectlyplasticmaterial.C
    Materials/Obsolete/j2plasticmaterial.C
    Materials/Obsolete/mplasticmaterial.C
    Materials/Obsolete/j2mplasticmaterial.C
    Materials/Obsolete/rankinepm.C
    )

set (sm_quasicontinuum
    Quasicontinuum/quasicontinuum.C
    Quasicontinuum/quasicontinuumnumberingscheme.C
    Quasicontinuum/quasicontinuumvtkxmlexportmodule.C    
    Quasicontinuum/fullsolveddomain.C
    Quasicontinuum/geometrygenerator.C
    )

set (sm
    ${sm_element}
    ${sm_interface_elements}
    ${sm_crosssection}
    ${sm_loads}
    ${sm_emodel}
    ${sm_material}
    ${sm_interface_materials}
    ${sm_status}
    ${sm_export}
    ${sm_barrier}
    ${sm_error}
    ${sm_xfem}
    ${sm_unsorted}
    ${sm_obsolete}
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

# Not included:
#pnldeidynamiccomm
