set (fm_emodel
    fluidmodel.C
    cbs.C
    supg.C
    stokesflow.C
    stokesflowvelocityhomogenization.C
    )

set (fm_element
    fmelement.C
    supgelement.C
    supgelement2.C
    cbselement.C
    tr1_2d_cbs.C
    tr1_2d_supg.C
    tr1_2d_supg2.C
    tr1_2d_supg_axi.C
    tr1_2d_supg2_axi.C
    tet1_3d_supg.C
    tr21_2d_supg.C
    quad10_2d_supg.C
    tr21stokes.C
    tet21stokes.C
    hexa21stokes.C
    tr1bubblestokes.C
    tet1bubblestokes.C
    line2boundaryelement.C
    fluidmaterialevaluator.C
    )

set (fm_material
    fluidcrosssection.C
    newtonianfluid.C
    fe2fluidmaterial.C
    nonlinearfluidmaterial.C
    )

set (fm_loads
    tractionpressurebc.C
    )

set (fm_status
    fluiddynamicmaterial.C
    twofluidmaterial.C
    binghamfluid2.C
    )

set (fm_interface
    leplic.C
    )

set (fm_levelset
    levelsetpcs.C
    fastmarchingmethod.C
    )

if (USE_PFEM)
    set (fm_pfem
        pfem/delaunaytriangulator.C
	pfem/edge2d.C
	pfem/delaunaytriangle.C
	pfem/pfem.C
	pfem/pfemelement.C
	pfem/pfemelement2d.C
	pfem/pfemparticle.C
	pfem/pfemnumberingschemes.C
	pfem/delaunaytriangle.C
	pfem/tr1_2d_pfem.C
	pfem/fluidstructureproblem.C
	pfem/interactionload.C
	pfem/interactionpfemparticle.C
	pfem/interactionboundarycondition.C
	)
endif ()

set (fm
    ${fm_emodel}
    ${fm_element}
    ${fm_material}
    ${fm_status}
    ${fm_loads}
    ${fm_interface}
    ${fm_levelset}
    ${fm_new}
    ${fm_pfem}
    )

add_library (fm
    ${LIB_TYPE}
    ${fm}
    )
