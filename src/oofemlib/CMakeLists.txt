# Much to organize
set (core_unsorted
    femcmpnn.C domain.C timestep.C metastep.C gausspnt.C
    cltypes.C timer.C dictionary.C freestor.C
    connectivitytable.C compiler.C error.C mathfem.C logger.C util.C
    initmodulemanager.C initmodule.C initialcondition.C
    homogenize.C
    localgaussianrandomgenerator.C
    nonlocalbarrier.C
    trianglemesherinterface.C
    geotoolbox.C geometry.C
    externalfieldgenerator.C
    datastream.C
    rvematerial.C
    weakperiodicbc.C
    # Semi sorted:
    errorestimator.C meshqualityerrorestimator.C remeshingcrit.C
    sloangraph.C sloangraphnode.C sloanlevelstruct.C
    eleminterpunknownmapper.C primaryunknownmapper.C materialmappingalgorithm.C
    nonlocalmaterialext.C randommaterialext.C
    inputrecord.C oofemtxtinputrecord.C dynamicinputrecord.C
    dynamicdatareader.C buffereddatareader.C oofemtxtdatareader.C tokenizer.C parser.C
    spatiallocalizer.C dummylocalizer.C octreelocalizer.C
    integrationrule.C gaussintegrationrule.C lobattoir.C
    smoothednodalintvarfield.C dofmanvalfield.C
    # Deprecated?
    spoolessparsemtrx.C
    rowcol.C skyline.C skylineu.C
    ldltfact.C
    dss.C dsssolver.C
    inverseit.C subspaceit.C gjacobi.C
    )

#oofegutils

set (core_iml
    iml/dyncomprow.C iml/dyncompcol.C iml/symcompcol.C iml/compcol.C
    iml/precond.C iml/voidprecond.C iml/icprecond.C iml/iluprecond.C iml/ilucomprowprecond.C iml/diagpre.C
    iml/imlsolver.C
    )

set (core_oofeg oofeggraphiccontext.C)

set (core_nodalrecovery
    nodalrecoverymodel.C
    zznodalrecoverymodel.C
    nodalaveragingrecoverymodel.C
    sprnodalrecoverymodel.C
    )

set (core_dofs
    dof.C
    masterdof.C
    simpleslavedof.C
    slavedof.C
    activedof.C
    )

set (core_nm
    sparselinsystemnm.C
    nrsolver.C
    linesearch.C
    spoolessolver.C
    calmls.C
    slepcsolver.C
    petscsolver.C
    )

set (core_fields
    fieldmanager.C
    field.C
    primaryfield.C
    intvarfield.C
    maskedprimaryfield.C
    dofdistributedprimaryfield.C
    )

set (core_ltf
    loadtimefunction.C
    constantfunction.C
    heavisideltf.C
    peakfunction.C
    piecewiselinfunction.C
    periodicpiecewiselinfunction.C
    userdefinedloadtimefunction.C
    )

set (core_cs
    crosssection.C
    emptycs.C
    simplecrosssection.C
    )

set (core_td
    #particletopologydescription.C
    )

set (core_dofman
    dofmanager.C
    node.C
    slavenode.C
    hangingnode.C
    rigidarmnode.C
    elementinternaldofman.C
    elementside.C
    )

set (core_material
    material.C
    )

set (core_export
    outputmanager.C
    exportmodule.C
    exportmodulemanager.C
    vtkexportmodule.C
    vtkxmlexportmodule.C
    homexportmodule.C
    matlabexportmodule.C
    gpexportmodule.C
    )

set (core_iga
    iga/iga.C
    iga/feibspline.C
    iga/feinurbs.C
    iga/feitspline.C
    )

set (core_arrays
    intarray.C
    floatarray.C
    stressstrainbasevector.C
    stressvector.C
    strainvector.C
    matrix.C
    floatmatrix.C
    petscsparsemtrx.C
    )

set (core_engng
    engngm.C
    staggeredproblem.C
    )

set (core_element
    element.C
    )

set (core_fei
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
    fei2dquadbiquad.C
    fei3dlinelin.C
    fei3dtrquad.C
    fei3dtetlin.C
    fei3dtetquad.C
    fei3dhexalin.C
    fei3dhexaquad.C
    fei3dwedgelin.C 
    fei3dwedgequad.C
    )

set (core_xfem
    xfem/delaunay.C
    xfem/xfemelementinterface.C
    xfem/enrichmentfunction.C
    xfem/patch.C
    xfem/xfemmanager.C
    xfem/enrichmentitem.C
    xfem/patchintegrationrule.C
    xfem/enrichmentdomain.C
    )

set (core_parallel
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

set (core_bc
    generalboundarycondition.C
    boundarycondition.C
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
    nodalload.C
    reinforcement.C
    surfacetensionbc.C
    )

set (core
    ${core_iml} # TODO: Make this optional
    ${core_unsorted}
    ${core_oofeg}
    ${core_nodalrecovery}
    ${core_dofs}
    ${core_nm}
    ${core_td}
    ${core_fields}
    ${core_dofman}
    ${core_material}
    ${core_export}
    ${core_ltf}
    ${core_cs}
    ${core_arrays}
    ${core_engng}
    ${core_element}
    ${core_fei}
    ${core_xfem}
    ${core_parallel}
    ${core_bc}
    ${core_iga}
    )

add_library (core
    ${LIB_TYPE}
    ${core}
    )
