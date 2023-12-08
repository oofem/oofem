set (sm_material
    abaqususermaterial.C
    structuralmaterial.C
    linearelasticmaterial.C
    isolinearelasticmaterial.C
    anisolinearelasticmaterial.C
    ortholinearelasticmaterial.C
    BoneMaterials/trabbonematerial.C 
    BoneMaterials/trabbonenl3d.C 
    BoneMaterials/trabbonenl.C 
    BoneMaterials/trabbone3d.C 
    BoneMaterials/trabboneembed.C 
    BoneMaterials/trabbonenlembed.C 
    BoneMaterials/trabbonegrad3d.C
    ConcreteMaterials/concrete2.C 
    ConcreteMaterials/concrete3.C 
    ConcreteMaterials/concretedpm.C 
    ConcreteMaterials/concretedpm2.C
    ConcreteMaterials/idm1.C 
    ConcreteMaterials/idmnl1.C 
    ConcreteMaterials/idmgrad.C
    ConcreteMaterials/mazarsmodel.C 
    ConcreteMaterials/mazarsmodelnl.C
    ConcreteMaterials/mplasticmaterial2.C
    ConcreteMaterials/compodamagemat.C
    ConcreteMaterials/concretefcm.C
    ConcreteMaterials/concretefcmviscoelastic.C
    ConcreteMaterials/frcfcm.C
    ConcreteMaterials/frcfcmnl.C
    RheoChainMaterials/rheoChM.C 
    RheoChainMaterials/maxwellChM.C 
    RheoChainMaterials/kelvinChM.C
    RheoChainMaterials/kelvinChSolM.C
    RheoChainMaterials/cebfip78.C
    RheoChainMaterials/doublepowerlaw.C
    RheoChainMaterials/b3mat.C
    RheoChainMaterials/b3solidmat.C
    RheoChainMaterials/mps.C
    RheoChainMaterials/mpsdammat.C
    RheoChainMaterials/eurocode2creep.C
    structuralfe2material.C
    structuralslipfe2material.C
    rcsd.C 
    rcsde.C 
    rcsdnl.C
    m4.C
    microplane_m1.C
    druckerPragerPlasticitySM.C
    druckerpragercutmat.C
    j2mat.C
    mdm.C
    micromaterial.C
    misesmat.C 
    misesmatgrad.C 
    misesmatnl.C
    rankinemat.C
    rankinematgrad.C
    rankinematnl.C
    dustmat.C
    lsmastermat.C
    lsmastermatgrad.C
    simplevitrificationmaterial.C
    anisodamagemodel.C
    winklerpasternak.C
    winklermodel.C
    structmatsettable.C
    tutorialmaterial.C
    graddamagematerialextensioninterface.C    
    rcm2.C
    isodamagemodel.C
    steelrelaxmat.C
    fcm.C
    #Mixed pressure material interface
    MixedPressure/mixedpressurematerialextensioninterface.C
    #Hyperelastic materials
    HyperelasticMaterials/basehyperelasticmaterial.C    
    HyperelasticMaterials/mooneyrivlincompressiblematerial.C
    HyperelasticMaterials/blatzkomaterial.C
    HyperelasticMaterials/simopistermaterial.C
    HyperelasticMaterials/ogdencompressiblematerial.C
    HyperelasticMaterials/hyperelasticmaterial1d.C
    #Lattice models
    LatticeMaterials/latticematstatus.C
    LatticeMaterials/latticestructuralmaterial.C	
    LatticeMaterials/latticelinearelastic.C
    LatticeMaterials/latticedamage.C
    LatticeMaterials/latticeviscoelastic.C
    LatticeMaterials/latticedamageviscoelastic.C
    LatticeMaterials/latticebondplasticity.C
    LatticeMaterials/latticeplasticitydamage.C
    LatticeMaterials/latticeslip.C
    LatticeMaterials/latticeplasticitydamageviscoelastic.C
    1D/isoasymm1d.C
)

if (USE_MFRONT)
    list (APPEND sm_material mfrontusermaterial.C)
endif ()
    
if (USE_BOOSTPYTHON_BINDINGS AND USE_PYTHON_EXTENSION)
    list (APPEND sm_material structuralpythonmaterial.C)
endif ()

set (sm_status
    structuralms.C
    InterfaceMaterials/structuralinterfacematerialstatus.C
    deformationtheorymaterial.C
    microplanematerial.C 
    microplanematerial_bazant.C
    )

set (sm_interface_materials
    InterfaceMaterials/structuralinterfacematerial.C
    InterfaceMaterials/bondceb.C
    InterfaceMaterials/cohint.C
    InterfaceMaterials/intmatelastic.C
    InterfaceMaterials/intmatbilinczfagerstrom.C
    InterfaceMaterials/intmatbilinczjansson.C
    InterfaceMaterials/intmatbilinearcz.C
    InterfaceMaterials/intmatbilinczelastic.C
    InterfaceMaterials/intmatbilinczfagerstromrate.C
    InterfaceMaterials/intmatisodamage.C
    InterfaceMaterials/intmatisodamagetable.C
    InterfaceMaterials/intmatcoulombcontact.C
    InterfaceMaterials/intmatdummycz.C    
    InterfaceMaterials/Deprecated/isointerfacedamage01.C
    InterfaceMaterials/Deprecated/isointerfacedamage02.C
    InterfaceMaterials/Deprecated/simpleinterfacemat.C
    InterfaceMaterials/masonry02.C
    InterfaceMaterials/mat_cebfip90.C
    InterfaceMaterials/linkslip.C
    )

set (sm_obsolete
    Obsolete/steel1.C
    Obsolete/plasticmaterial.C
    Obsolete/perfectlyplasticmaterial.C
    Obsolete/j2plasticmaterial.C
    Obsolete/mplasticmaterial.C
    Obsolete/j2mplasticmaterial.C
    Obsolete/rankinepm.C
    )

set (sm_materials
    ${sm_material}
    ${sm_interface_materials}
    ${sm_status}
    ${sm_obsolete}
    )

add_library(sm_materials
    ${LIB_TYPE}
    ${sm_materials}
    )