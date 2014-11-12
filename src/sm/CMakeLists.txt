set (sm_element
    structuralelementevaluator.C
    space3delementevaluator.C
    planestresselementevaluator.C
    igaelements.C
    structuralelement.C
    nlstructuralelement.C
    lspace.C lspacebb.C
    qspace.C qspacegrad.C
    q27space.C
    ltrspace.C qtrspace.C qtrspacegrad.C
    lwedge.C qwedge.C qwedgegrad.C
    htselement.C
    truss1d.C truss2d.C truss3d.C
    trplanstrss.C trplanestrain.C trplanrot.C
    libeam2d.C libeam2dnl.C
    libeam3d.C libeam3dnl.C libeam3dnl2.C libeam3d2.C
    beam2d.C beam3d.C
    #interfaceelement1d.C interfaceelem2dquad.C interfaceelem3dtrlin.C
    quad1planestrain.C rershell.C tr_shell01.C qplanstrss.C q9planstrss.C qtrplstr.C q4axisymm.C
    planstrss.C
    planstrssxfem.C
    trplanstrssxfem.C
    axisymm3d.C l4axisymm.C
    macrolspace.C
    cct.C
    cohsur3d.C
    cct3d.C
    quad1mindlin.C
    quad1mindlinshell3d.C
    trplanrot3d.C
    lumpedmasselement.C
    springelement.C
    qtrplstrgrad.C
    qplanestressgrad.C
    latticestructuralelement.C
    lattice2d.C
    shell7base.C
    tr2shell7.C
    shell7basexfem.C
    tr2shell7xfem.C
    trplanestressrotallman.C
    structuralinterfaceelement.C
    tet21ghostsolid.C
    graddpelement.C graddpmaterialextensioninterface.C
    qtruss1d.C qtruss1dgrad.C
    qplanestrain.C qplanestraingrad.C
    qtrplanestrain.C qtrplanestraingrad.C
    qspacegrad.C
    dkt.C qdkt.C
    quad1platesubsoil.C tria1platesubsoil.C
    
    trwarp.C
    #qplanstrssphf.C
    )

set (sm_interface_elements
    InterfaceElements/intelpoint.C
    InterfaceElements/intelline1.C
    InterfaceElements/intelline2.C
    InterfaceElements/intelsurftr1.C
    InterfaceElements/Deprecated/interfaceelement1d.C
    InterfaceElements/Deprecated/interfaceelem2dlin.C
    InterfaceElements/Deprecated/interfaceelem2dquad.C
    InterfaceElements/Deprecated/interfaceelem3dtrlin.C
    )

set (sm_boundary_conditions
    prescribedgenstrainshell7.C
    )

set (sm_crosssection
    structuralcrosssection.C
    simplecrosssection.C
    layeredcrosssection.C
    fiberedcs.C
    structuralinterfacecrosssection.C
    variablecrosssection.C
    )

set (sm_emodel
    structengngmodel.C
    staticstructural.C
    linearstatic.C
    nlinearstatic.C
    xfemstatic.C
    nlineardynamic.C
    eigenvaluedynamic.C
    linearstability.C
    adaptlinearstatic.C
    adaptnlinearstatic.C
    deidynamic.C
    diidynamic.C
    nldeidynamic.C
    incrementallinearstatic.C
    freewarping.C
    )

set (sm_loads
    tf1.C
    usrdeftempfield.C
    rotatingboundary.C
    constantpressureload.C
    structtemperatureload.C
    structeigenstrainload.C
    )

set (sm_material
    abaqususermaterial.C
    structuralmaterial.C
    structuralinterfacematerial.C
    linearelasticmaterial.C
    isolinearelasticmaterial.C
    ortholinearelasticmaterial.C
    anisolinearelasticmaterial.C
    concrete2.C concrete3.C concretedpm.C concretedpm2.C
    cebfip78.C
    doublepowerlaw.C
    b3mat.C
    b3solidmat.C
    rcsd.C rcsde.C rcsdnl.C
    m4.C
    idm1.C idmnl1.C idmgrad1.C
    mazarsmodel.C mazarsmodelnl.C
    druckerPragerPlasticitySM.C
    druckerpragercutmat.C
    masonry02.C
    j2mat.C
    mat_cebfip90.C
    mdm.C
    compodamagemat.C
    micromaterial.C
    hyperelasticmaterial.C
    misesmat.C misesmatgrad.C misesmatnl.C
    trabbonematerial.C trabbonenl3d.C trabbonenl.C trabbone3d.C trabboneembed.C trabbonenlembed.C trabbonegrad3d.C
    cohint.C
    mps.C
    rankinemat.C
    rankinematgrad.C
    rankinematnl.C
    kelvinChSolM.C
    latticematstatus.C
    latticedamage2d.C
    dustmat.C
    lsmastermat.C
    lsmastermatgrad.C
    InterfaceMaterials/intmatdummycz.C
    simplevitrificationmaterial.C
    anisodamagemodel.C
    winklerpasternak.C
    structmatsettable.C
    bondceb.C
    mooneyrivlin.C
	)

set (sm_interface_materials
    InterfaceMaterials/intmatbilinczfagerstrom.C
    InterfaceMaterials/intmatbilinczjansson.C
    InterfaceMaterials/intmatbilinearcz.C
    InterfaceMaterials/intmatbilinczelastic.C
    InterfaceMaterials/intmatbilinczfagerstromrate.C
    InterfaceMaterials/intmatisodamage.C
    InterfaceMaterials/intmatisodamagetable.C
    InterfaceMaterials/intmatcoulombcontact.C
    InterfaceMaterials/Deprecated/isointerfacedamage01.C
    InterfaceMaterials/Deprecated/isointerfacedamage02.C
    InterfaceMaterials/Deprecated/simpleinterfacemat.C
    )

set (sm_error
    huertaerrorestimator.C
    zzerrorestimator.C
    combinedzzsiee.C
    scalarerrorindicator.C
    directerrorindicatorrc.C
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
    structuralms.C
	structuralinterfacematerialstatus.C
    rheoChM.C maxwellChM.C kelvinChM.C
    rcm2.C
    isodamagemodel.C
    mplasticmaterial2.C
    deformationtheorymaterial.C
    microplanematerial.C microplanematerial_bazant.C
    )

set (sm_xfem
	xfem/xfemstructuralelementinterface.C
	xfem/enrichmentitems/crack.C
	xfem/enrichmentitems/delamination.C
	xfem/enrichmentitems/shellcrack.C
	xfem/propagationlaws/plhoopstresscirc.C
	xfem/xfemstructuremanager.C
	xfemsolverinterface.C
	)

set (sm_unsorted
    particle.C
    refinedelement.C
    refinedmesh.C
    microplane.C
    gpinitmodule.C
    stressstrainbasevector.C stressvector.C strainvector.C
    mappers/primvarmapper.C
    )

set (sm_parallel
    feticommunicator.C
    fetiboundarydofman.C
    fetisolver.C
    )

set (sm_obsolete
    steel1.C
    plasticmaterial.C
    perfectlyplasticmaterial.C
    j2plasticmaterial.C
    mplasticmaterial.C
    j2mplasticmaterial.C
    rankinepm.C
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
