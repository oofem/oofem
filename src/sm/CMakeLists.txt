set (sm_element
    lspace.C lspacebb.C qspace.C ltrspace.C
    truss1d.C truss2d.C truss3d.C
    trplanstrss.C trplanestrain.C trplanrot.C
    libeam2d.C libeam2dnl.C
    libeam3d.C libeam3dnl.C libeam3dnl2.C libeam3d2.C
    beam2d.C beam3d.C
    interfaceelement1d.C interfaceelem2dquad.C interfaceelem3dtrlin.C
    quad1planestrain.C rershell.C qplanstrss.C qtrplstr.C q4axisymm.C
    planstrss.C
    planstrssxfem.C
    axisymm3d.C l4axisymm.C
    macrolspace.C
    cct.C
    cohsur3d.C
    cct3d.C
    trplanrot3d.C
    lumpedmasselement.C
    tr_shell01.C
    springelement.C
    qtrplstrgrad.C
    qplanestressgrad.C
    latticestructuralelement.C
    lattice2d.C
    )

set (sm_crosssection
    layeredcrosssection.C
    fiberedcs.C
    )

set (sm_emodel
    linearstatic.C
    nlinearstatic.C
    nlineardynamic.C
    eigenvaluedynamic.C
    linearstability.C
    adaptlinearstatic.C
    adaptnlinearstatic.C
    deidynamic.C
    diidynamic.C
    nldeidynamic.C
    incrementallinearstatic.C
    )

set (sm_loads
    tf1.C
    usrdeftempfield.C
    rotatingboundary.C
    )

set (sm_material
    ortholinearelasticmaterial.C
    concrete2.C concrete3.C concretedpm.C concretedpm2.C
    cebfip78.C
    doublepowerlaw.C
    b3mat.C
    b3solidmat.C
    rcsd.C rcsde.C rcsdnl.C
    m4.C
    idm1.C idmnl1.C
    mazarsmodel.C mazarsmodelnl.C
    druckerPragerPlasticitySM.C
    druckerpragercutmat.C
    masonry02.C
    isointerfacedamage01.C
    j2mat.C
    mat_cebfip90.C
    hellmat.C
    mdm.C
    compodamagemat.C
    micromaterial.C
    hyperelasticmaterial.C
    misesmat.C misesmatgrad.C misesmatnl.C
    graddpelement.C qtruss1d.C qtruss1dgrad.C qplanestrain.C qplanestraingrad.C qtrplanestrain.C qtrplanestraingrad.C qspacegrad.C
    trabbonematerial.C trabbonenl3d.C trabbonenl.C trabbone3d.C trabboneembed.C trabbonenlembed.C
    cohint.C
    mps.C
    rankinemat.C
    rankinematgrad.C
    rankinematnl.C
    kelvinChSolM.C
    simpleinterfacemat.C
    latticematstatus.C
    latticedamage2d.C
    dustmat.C
    )

set (sm_error
    huertaerrorestimator.C
    zzerrorestimator.C
    combinedzzsiee.C
    scalarerrorindicator.C
    directerrorindicatorrc.C
    )

set (sm_export
    dmexportmodule.C
    poiexportmodule.C
    gpexportmodule.C
    )

set (sm_barrier
    polylinenonlocalbarrier.C
    symmetrybarrier.C
    )

set (sm_mesher
    t3dinterface.C
    targe2interface.C
    freeminterface.C
    subdivision.C
    )

set (sm_status
    rheoChM.C maxwellChM.C kelvinChM.C
    rcm2.C
    isodamagemodel.C
    hydram.C
    mplasticmaterial2.C
    deformationtheorymaterial.C
    microplanematerial.C microplanematerial_bazant.C
    )

set (sm_unsorted
    particle.C
    refinedelement.C
    refinedmesh.C
    microplane.C
    gpinitmodule.C
    )

set (sm_mapping
    mmaclosestiptransfer.C
    mmashapefunctprojection.C
    mmaleastsquareprojection.C
    mmacontainingelementprojection.C
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

add_library(sm
    ${sm_element}
    ${sm_crosssection}
    ${sm_loads}
    ${sm_emodel}
    ${sm_material}
    ${sm_status}
    ${sm_export}
    ${sm_barrier}
    ${sm_error}
    ${sm_mesher}
    ${sm_mapping}
    ${sm_unsorted}
    ${sm_parallel}
    ${sm_obsolete}
    ${sm_new}
    )

# Not included:
#pnldeidynamiccomm
