# Much to organize
set (core_unsorted
    classfactory.C
    femcmpnn.C domain.C timestep.C metastep.C gausspoint.C
    cltypes.C timer.C dictionary.C
    connectivitytable.C error.C mathfem.C logger.C util.C
    initmodulemanager.C initmodule.C initialcondition.C
    assemblercallback.C
    homogenize.C
    nonlocalbarrier.C
    geotoolbox.C geometry.C
    datastream.C
    set.C
    weakperiodicbc.C
    solutionbasedshapefunction.C
    fracturemanager.C
    neumannmomentload.C
    parallelcontext.C
    contact/contactmanager.C
    contact/contactdefinition.C
    contact/contactelement.C
    prescribedmean.C
    # Semi sorted:
    errorestimator.C meshqualityerrorestimator.C remeshingcrit.C
    sloangraph.C sloangraphnode.C sloanlevelstruct.C
    eleminterpunknownmapper.C primaryunknownmapper.C materialmappingalgorithm.C
    nonlocalmaterialext.C randommaterialext.C
    inputrecord.C oofemtxtinputrecord.C dynamicinputrecord.C
    dynamicdatareader.C oofemtxtdatareader.C tokenizer.C parser.C
    spatiallocalizer.C dummylocalizer.C octreelocalizer.C
    integrationrule.C gaussintegrationrule.C lobattoir.C
    smoothednodalintvarfield.C dofmanvalfield.C
    # Deprecated?
    rowcol.C skyline.C skylineu.C
    ldltfact.C
    inverseit.C subspaceit.C gjacobi.C
    #
    symcompcol.C compcol.C
    )

#oofegutils

if (USE_SPOOLES)
    list (APPEND core_unsorted spoolessparsemtrx.C spoolessolver.C)
endif ()

if (USE_IML)
    list (APPEND core_unsorted
        iml/dyncomprow.C iml/dyncompcol.C
        iml/precond.C iml/voidprecond.C iml/icprecond.C iml/iluprecond.C iml/ilucomprowprecond.C iml/diagpre.C
        iml/imlsolver.C
        )
endif ()

if (USE_SLEPC)
    list (APPEND core_unsorted slepcsolver.C)
endif ()


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
    sparsenonlinsystemnm.C
    nrsolver.C
    dynamicrelaxationsolver.C
    linesearch.C
    calmls.C
    staggeredsolver.C
    )

if (USE_PETSC)
    list (APPEND core_nm petscsolver.C petscsparsemtrx.C)
endif ()

if (USE_MKLPARDISO)
    list (APPEND core_nm mklpardisosolver.C)
endif ()

set (core_fields
    fieldmanager.C
    field.C
    primaryfield.C
    intvarfield.C
    maskedprimaryfield.C
    dofdistributedprimaryfield.C
    )

set (core_ltf
    function.C
    scalarfunction.C
    constantfunction.C
    heavisidetimefunction.C
    peakfunction.C
    piecewiselinfunction.C
    periodicpiecewiselinfunction.C
    calculatorfunction.C
    localgaussianrandomfunction.C
    interpolatingfunction.C
    )

set (core_cs
    crosssection.C
    emptycs.C
    )

set (core_td
    particletopologydescription.C
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
    dummymaterial.C
    )

set (core_export
    outputmanager.C
    exportmodule.C
    exportmodulemanager.C
    outputexportmodule.C
    errorcheckingexportmodule.C
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
    floatmatrix.C
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
    feinterpol1d.C
    feinterpol2d.C
    feinterpol3d.C
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
    fei3dtrlin.C
    fei3dtrquad.C
    fei3dtetlin.C
    fei3dtetquad.C
    fei3dhexalin.C
    fei3dhexaquad.C
    fei3dhexatriquad.C
    fei3dwedgelin.C 
    fei3dwedgequad.C
    )

set (core_xfem
    xfem/delaunay.C
    xfem/xfemelementinterface.C
    xfem/enrichmentfunction.C
    xfem/xfemmanager.C
    xfem/enrichmentitem.C
    xfem/geometrybasedei.C
    xfem/hybridei.C
    xfem/directei.C
    xfem/listbasedei.C
    xfem/inclusion.C
    xfem/patchintegrationrule.C
    xfem/propagationlaw.C
    xfem/XFEMDebugTools.C
    xfem/enrichmentfronts/enrichmentfront.C
    xfem/enrichmentfronts/enrichmentfrontdonothing.C
    xfem/enrichmentfronts/enrichmentfrontextend.C
    xfem/enrichmentfronts/enrichmentfrontreducefront.C
    xfem/enrichmentfronts/enrichmentfrontlinbranchfuncradius.C
    xfem/enrichmentfronts/enrichmentfrontlinbranchfunconeel.C
    xfem/enrichmentfronts/enrichmentfrontcohesivebranchfunconeel.C
    xfem/enrichmentfronts/enrichmentfrontintersection.C
    xfem/integrationrules/discsegintegrationrule.C
    xfem/xfemtolerances.C
    )

set (core_parallel
    communicator.C
    processcomm.C
    problemcomm.C
    combuff.C
    domaintransactionmanager.C
    parallelordering.C
    dyncombuff.C
    loadbalancer.C
    wallclockloadbalancermonitor.C
    nonlocalmatwtp.C
    )

if (USE_PARMETIS)
    list (APPEND core_parallel parmetisloadbalancer.C)
endif ()

set (core_bc
    generalboundarycondition.C
    boundarycondition.C
    prescribedgradienthomogenization.C
    prescribedgradient.C
    prescribedgradientbcneumann.C
    prescribedgradientbcweak.C
    prescribedgradientbcweakperiodic.C
    prescribedgradientbcweakdirichlet.C
    prescribedgradientbcperiodic.C
    mixedgradientpressurebc.C
    mixedgradientpressuredirichlet.C
    mixedgradientpressureneumann.C
    mixedgradientpressureweakperiodic.C
    load.C
    boundaryload.C
    pointload.C
    linearedgeload.C
    constantedgeload.C
    constantsurfaceload.C
    nodalload.C
    deadweight.C
    reinforcement.C
    surfacetensionbc.C
    linearconstraintbc.C
    )

if (USE_PYTHON)
    list (APPEND core_bc userdefdirichletbc.C)
    list (APPEND core_unsorted pythonexpression.C)
endif ()

set (core_mapping
    mmaclosestiptransfer.C
    mmashapefunctprojection.C
    mmaleastsquareprojection.C
    mmacontainingelementprojection.C
    matstatmapperint.C
    )

set (core_mesher
    t3dinterface.C
    targe2interface.C
    freeminterface.C
    subdivision.C
    trianglemesherinterface.C
    )

set (core
    ${core_unsorted}
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
    ${core_bc}
    ${core_iga}
    ${core_mesher}
    ${core_mapping}
    )

if (USE_PARALLEL)
    list (APPEND core ${core_parallel})
endif ()

if (USE_OOFEG)
    list (APPEND core ${core_oofeg})
endif ()


add_library (core
    ${LIB_TYPE}
    ${core}
    )
