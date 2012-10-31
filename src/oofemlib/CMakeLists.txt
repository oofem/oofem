include_directories(${CMAKE_CURRENT_SOURCE_DIR}/xfem)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/iga)

# Much to organize
set (oofemlib_unsorted
    femcmpnn.C domain.C timestep.C metastep.C gausspnt.C
    cltypes.C clock.C dictionr.C freestor.C
    conTable.C compiler.C error.C mathfem.C logger.C util.C
    initmodulemanager.C initmodule.C initial.C
    homogenize.C
    localgaussianrandomgenerator.C
    nonlocalbarrier.C
    trianglemesherinterface.C
    geotoolbox.C geometry.C
    externalfieldgenerator.C
    # Semi sorted:
    errorestimator.C meshqualityerrorestimator.C remeshingcrit.C
    oofemtxtinputrecord.C buffereddatareader.C oofemtxtdatareader.C
    sloangraph.C sloangraphnode.C sloanlevelstruct.C
    eleminterpunknownmapper.C primaryunknownmapper.C materialmappingalgorithm.C
    structuralms.C
    nonlocalmaterialext.C randommaterialext.C
    inputrecord.C tokenizer.C parser.C
    spatiallocalizer.C dummylocalizer.C octreelocalizer.C
    integrationrule.C gaussintegrationrule.C lobattoir.C
    smoothednodalintvarfield.C dofmanvalfield.C
    # Deprecated?
    strreader.C
    spoolessparsemtrx.C
    rowcol.C skyline.C skylineu.C
    ldltfact.C
    dss.C dsssolver.C
    inverseit.C subspaceit.C gjacobi.C
    )

#oofegutils
#structeigenstrainload

set (oofemlib_iml
    iml/dyncomprow.C iml/dyncompcol.C iml/symcompcol.C iml/compcol.C
    iml/precond.C iml/voidprecond.C iml/icprecond.C iml/iluprecond.C iml/ilucomprowprecond.C iml/diagpre.C
    iml/imlsolver.C
    )

set (oofemlib_oofeg oofeggraphiccontext.C)

set (oofemlib_nodalrecovery
    nodalrecoverymodel.C
    zznodalrecoverymodel.C
    nodalaveragingrecoverymodel.C
    sprnodalrecoverymodel.C
    )

set (oofemlib_dofs
    dof.C
    masterdof.C
    simpleslavedof.C
    slavedof.C
    activedof.C
    )

set (oofemlib_nm
    sparselinsystemnm.C
    nrsolver.C
    nrsolver2.C
    linesearch.C
    spoolessolver.C
    calmls.C
    slepcsolver.C
    petscsolver.C
    )

set (oofemlib_fields
    fieldmanager.C
    field.C
    primaryfield.C
    intvarfield.C
    maskedprimaryfield.C
    dofdistributedprimaryfield.C
    )

set (oofemlib_ltf
    loadtime.C
    constant.C
    heavisideltf.C
    peak.C
    piecewis.C
    piecewisper.C
    piecewisblock.C
    usrdeftimefunct.C
    )

set (oofemlib_cs
    crosssection.C
    emptycs.C
    simplecrosssection.C
    structuralcrosssection.C
    )

set (oofemlib_td
    #particletopologydescription.C
    )

set (oofemlib_dofman
    dofmanager.C
    node.C
    slavenode.C
    hangingnode.C
    rigidarmnode.C
    elementinternaldofman.C
    elementside.C
    )

set (oofemlib_material
    material.C
    structuralmaterial.C
    isolinearelasticmaterial.C
    linearelasticmaterial.C
    )

set (oofemlib_export
    outputmanager.C
    exportmodule.C
    exportmodulemanager.C
    vtkexportmodule.C
    vtkxmlexportmodule.C
    homexportmodule.C
    )

set (oofemlib_iga
    iga/iga.C
    iga/feibspline.C
    iga/feinurbs.C
    iga/feitspline.C
    iga/structuralelementevaluator.C
    iga/planestresselementevaluator.C
    iga/space3delementevaluator.C
    iga/igaelements.C
    )

set (oofemlib_arrays
    intarray.C
    flotarry.C
    stressstrainbasevector.C
    stressvector.C
    strainvector.C
    matrix.C
    flotmtrx.C
    petscsparsemtrx.C
    )

set (oofemlib_engng
    engngm.C
    structengngmodel.C
    )

set (oofemlib_element
    element.C
    structuralelement.C
    nlstructuralelement.C
    )

set (oofemlib_fei
    feinterpol.C
    fei1dlin.C
    fei1dquad.C
    fei1dhermite.C
    fei2dtrlin.C
    fei2dtrquad.C
    fei2dlinehermite.C
    fei2dlinelin.C
    fei2dlinequad.C
    fei2dquadconst.C
    fei2dquadlin.C
    fei2dquadquad.C
    fei3dlinelin.C
    fei3dtrlin.C
    fei3dtetquad.C
    fei3dhexalin.C
    fei3dhexaquad.C
    )

set (oofemlib_xfem
    xfem/delaunay.C
    xfem/xfemelementinterface.C
    xfem/enrichmentfunction.C
    xfem/patch.C
    xfem/xfemmanager.C
    xfem/enrichmentitem.C
    xfem/patchintegrationrule.C
    )

set (oofemlib_parallel
    communicator.C
    processcomm.C
    problemcomm.C
    combuff.C
    domaintransactionmanager.C
    petscordering.C
    petsccontext.C
    dyncombuff.C
    loadbalancer.C
    parmetisloadbalancer.C
    nonlocalmatwtp.C
    )

set (oofemlib_bc
    generalbc.C
    boundary.C
    prescribedgradient.C
    mixedgradientpressurebc.C
    mixedgradientpressuredirichlet.C
    mixedgradientpressureneumann.C
    load.C
    boundaryload.C
    pointload.C
    linearedgeload.C
    constantedgeload.C
    constantsurfaceload.C
    nodload.C
    reinforcement.C
    surfacetensionbc.C
    )

add_library(oofemlib
    ${oofemlib_iml} # TODO: Make this optional
    ${oofemlib_unsorted}
    ${oofemlib_oofeg}
    ${oofemlib_nodalrecovery}
    ${oofemlib_dofs}
    ${oofemlib_nm}
    ${oofemlib_td}
    ${oofemlib_fields}
    ${oofemlib_dofman}
    ${oofemlib_material}
    ${oofemlib_export}
    ${oofemlib_ltf}
    ${oofemlib_cs}
    ${oofemlib_arrays}
    ${oofemlib_engng}
    ${oofemlib_element}
    ${oofemlib_fei}
    ${oofemlib_xfem}
    ${oofemlib_parallel}
    ${oofemlib_bc}
    ${oofemlib_iga}
    )
