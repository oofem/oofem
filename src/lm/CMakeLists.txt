set (lm_elements
    Elements/latticestructuralelement.C
    Elements/lattice2d.C
    Elements/lattice2dboundary.C
    Elements/lattice3d.C
    Elements/lattice3dboundary.C
    Elements/lattice3dboundarytruss.C
    Elements/latticelink3d.C
    Elements/latticelink3dboundary.C
    Elements/lattice3dnl.C
    Elements/latticelink3dnl.C
    Elements/latticetransportelement.C
    Elements/lattice2d_mt.C
    Elements/lattice3d_mt.C
    Elements/lattice3dboundary_mt.C
    )

set (lm_materials
    Materials/latticematstatus.C
    Materials/latticestructuralmaterial.C
    Materials/latticelinearelastic.C
    Materials/latticedamage.C
    Materials/latticeviscoelastic.C
    Materials/latticedamageviscoelastic.C
    Materials/latticebondplasticity.C
    Materials/latticeplasticitydamage.C
    Materials/latticeslip.C
    Materials/latticeplasticitydamageviscoelastic.C
    Materials/latticedamageshell.C
    Materials/latticeframeelastic.C
    Materials/latticeframeconcreteplastic.C
    Materials/latticeframesteelplastic.C
    Materials/latticetransmat.C
    )

set (lm_crosssections
    CrossSections/latticecrosssection.C
    )

set (lm_bc
    latticedirichletcoupling.C
    latticeneumanncoupling.C
    )

set (lm_export
    vtkxmllatticeexportmodule.C
    )

set (lm
    ${lm_elements}
    ${lm_materials}
    ${lm_crosssections}
    ${lm_bc}
    ${lm_export}
    )

add_library (lm
    ${LIB_TYPE}
    ${lm}
    )
