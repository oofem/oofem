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
    Elements/3D/q27space.C
    Elements/3D/ltrspace.C 
    Elements/3D/ltrspaceboundary.C
    Elements/3D/ltrspaceboundarytruss.C
    Elements/3D/ltrspaceboundarymembrane.C
    Elements/3D/ltrspaceboundarybeam.C
    Elements/3D/ltrspaceboundaryplate.C
    Elements/3D/ltrspaceboundaryvoigt.C
    Elements/3D/qtrspace.C 
    Elements/3D/lwedge.C 
    Elements/3D/qwedge.C 
    Elements/3D/macrolspace.C
    Elements/Bars/truss1d.C 
    Elements/Bars/truss2d.C 
    Elements/Bars/truss3d.C
    Elements/Bars/truss3dnl.C
    Elements/Bars/truss3dnl2.C
    Elements/Bars/qtruss1d.C 
    Elements/Beams/libeam2d.C 
    Elements/Beams/libeam2dnl.C
    Elements/Beams/libeam3d.C 
    Elements/Beams/libeam3dnl.C 
    Elements/Beams/libeam3dnl2.C 
    Elements/Beams/libeam3d2.C
    Elements/Beams/beambaseelement.C
    Elements/Beams/beam2d.C 
    Elements/Beams/beam3d.C
    Elements/Beams/libeam3dboundary.C
    Elements/Beams/libeam3dboundarytruss.C
    Elements/Beams/libeam3dboundarymembrane.C
    Elements/Beams/libeam3dboundarybeam.C
    Elements/Beams/libeam3dboundaryplate.C
    Elements/Beams/libeam3dboundaryvoigt.C
    Elements/Plates/cct.C
    Elements/Plates/dkt.C 
    Elements/Plates/qdkt.C
    Elements/Plates/quad1mindlin.C
    Elements/Plates/dkt3d.C
    Elements/Shells/rershell.C 
    Elements/Shells/tr_shell01.C 
    Elements/Shells/tr_shell02.C
    Elements/Shells/tr_shell11.C
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
    Elements/PlaneStrain/qtrplanestrain.C 
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
    Elements/PlaneStress/trplanrot3d.C
    Elements/PlaneStress/trplanestressrotallman.C
    Elements/PlaneStress/basicelement.C
    Elements/PlaneStress/trplanestressrotallman3d.C
    Elements/PlaneStress/linquad3d_planestress.C
    Elements/PlaneStress/qtrplstrslip.C
    Elements/PlaneStress/qplanstrssslip.C
    Elements/AbaqusUserElement.C
    Elements/htselement.C
    Elements/lumpedmasselement.C
    Elements/springelement.C
    Elements/LatticeElements/latticestructuralelement.C
    Elements/LatticeElements/lattice2d.C
    Elements/LatticeElements/lattice2dboundary.C
    Elements/LatticeElements/lattice3d.C
    Elements/LatticeElements/lattice3dboundary.C
    Elements/LatticeElements/lattice3dboundarytruss.C
    Elements/LatticeElements/latticelink3d.C
    Elements/LatticeElements/latticelink3dboundary.C
    Elements/LatticeElements/latticebeam3d.C
    Elements/LatticeElements/latticebeam3dboundary.C
    Elements/tet21ghostsolid.C
    Elements/quad1platesubsoil.C 
    Elements/quad2platesubsoil.C
    Elements/tria1platesubsoil.C
    Elements/tria2platesubsoil.C
    Elements/linedistributedspring.C
    Elements/trwarp.C
    Elements/PlaneStress/basicelement.C
    Elements/nodalspringelement.C
    #Gradient damage elements
    Elements/GradientDamage/graddamageelement.C 
    Elements/GradientDamage/Bars/truss1dgraddamage.C
    Elements/GradientDamage/Bars/qtruss1dgraddamage.C
    Elements/GradientDamage/PlaneStrain/quad1planestraingraddamage.C 
    Elements/GradientDamage/PlaneStrain/qplanestraingraddamage.C
    Elements/GradientDamage/PlaneStrain/qtrplanestraingraddamage.C
    Elements/GradientDamage/PlaneStress/trplstrgraddamage.C 
    Elements/GradientDamage/PlaneStress/qtrplstrgraddamage.C
    Elements/GradientDamage/PlaneStress/planestressgraddamage.C
    Elements/GradientDamage/PlaneStress/qplanestressgraddamage.C
    Elements/GradientDamage/3D/qspacegraddamage.C
    Elements/GradientDamage/3D/qtrspacegraddamage.C
    Elements/GradientDamage/3D/qspacegraddamage.C
    Elements/GradientDamage/3D/qwedgegraddamage.C
    #Mixed pressure elements
    Elements/MixedPressure/basemixedpressureelement.C
    Elements/MixedPressure/PlaneStrain/qtrplanestrainp1.C
    Elements/MixedPressure/PlaneStrain/quad1planestrainp1.C

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
    Elements/Interfaces/intelsurfquad1.C
    Elements/Interfaces/bondlink3d.C
    Elements/Interfaces/bondlink3dboundary.C
        #    Elements/bondlink3dboundarybeam.C	
#   	 Elements/bondlink3dboundaryplate.C		
    Elements/Interfaces/Deprecated/interfaceelement1d.C
    Elements/Interfaces/Deprecated/interfaceelem2dlin.C
    Elements/Interfaces/Deprecated/interfaceelem2dquad.C
    Elements/Interfaces/Deprecated/interfaceelem3dtrlin.C
)

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
    Materials/ConcreteMaterials/idmgrad.C
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
    Materials/structuralslipfe2material.C
    Materials/rcsd.C 
    Materials/rcsde.C 
    Materials/rcsdnl.C
    Materials/m4.C
    Materials/microplane_m1.C
    Materials/druckerPragerPlasticitySM.C
    Materials/druckerpragercutmat.C
    Materials/j2mat.C
    Materials/mdm.C
    Materials/micromaterial.C
    Materials/misesmat.C 
    Materials/misesmatgrad.C 
    Materials/misesmatnl.C
    Materials/rankinemat.C
    Materials/rankinematgrad.C
    Materials/rankinematnl.C
    Materials/dustmat.C
    Materials/lsmastermat.C
    Materials/lsmastermatgrad.C
    Materials/simplevitrificationmaterial.C
    Materials/anisodamagemodel.C
    Materials/winklerpasternak.C
    Materials/winklermodel.C
    Materials/structmatsettable.C
    Materials/tutorialmaterial.C
    Materials/graddamagematerialextensioninterface.C    
    Materials/rcm2.C
    Materials/isodamagemodel.C
    Materials/steelrelaxmat.C
    Materials/fcm.C
    #Mixed pressure material interface
    Materials/MixedPressure/mixedpressurematerialextensioninterface.C
    #Hyperelastic materials
    Materials/HyperelasticMaterials/basehyperelasticmaterial.C    
    Materials/HyperelasticMaterials/mooneyrivlincompressiblematerial.C
    Materials/HyperelasticMaterials/blatzkomaterial.C
    Materials/HyperelasticMaterials/simopistermaterial.C
    Materials/HyperelasticMaterials/ogdencompressiblematerial.C
    Materials/HyperelasticMaterials/hyperelasticmaterial1d.C
    #Lattice models
    Materials/LatticeMaterials/latticematstatus.C
    Materials/LatticeMaterials/latticestructuralmaterial.C	
    Materials/LatticeMaterials/latticelinearelastic.C
    Materials/LatticeMaterials/latticedamage.C
    Materials/LatticeMaterials/latticeviscoelastic.C
    Materials/LatticeMaterials/latticedamageviscoelastic.C
    Materials/LatticeMaterials/latticebondplasticity.C
    Materials/LatticeMaterials/latticeplasticitydamage.C
    Materials/LatticeMaterials/latticeslip.C
    Materials/LatticeMaterials/latticeplasticitydamageviscoelastic.C
    Materials/1D/isoasymm1d.C
)

if (USE_MFRONT)
    list (APPEND sm_material Materials/mfrontusermaterial.C)
endif ()

if (USE_BOOSTPYTHON_BINDINGS AND USE_PYTHON_EXTENSION)
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
    Materials/InterfaceMaterials/linkslip.C
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
    Contact/ActiveBc/node2nodepenaltycontact.C
    Contact/ActiveBc/node2nodelagrangianmultipliercontact.C
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
