/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "dummymaterial.h"
#include "isolinearelasticmaterial.h"

#ifdef __SM_MODULE
 #include "ortholinearelasticmaterial.h"
 #include "perfectlyplasticmaterial.h"
 #include "steel1.h"
 #include "concrete2.h"
 #include "concrete3.h"
 #include "cebfip78.h"
 #include "doublepowerlaw.h"
 #include "b3mat.h"
 #include "b3solidmat.h"
 #include "mps.h"
 #include "j2plasticmaterial.h"
 #include "rcsd.h"
 #include "rcsde.h"
 #include "rcsdnl.h"
 #include "m4.h"
 #include "idm1.h"
 #include "idmnl1.h"
 #include "mazarsmodel.h"
 #include "mazarsmodelnl.h"
 #include "druckerPragerPlasticitySM.h"
 #include "druckerpragermat.h"
 #include "j2mplasticmaterial.h"
 #include "rankinepm.h"
 #include "masonry02.h"
 #include "isointerfacedamage01.h"
 #include "j2mat.h"
 #include "mat_cebfip90.h"
 #include "hellmat.h"
 #include "mdm.h"
 #include "compodamagemat.h"
 #include "micromaterial.h"
 #include "hyperelasticmaterial.h"
 #include "misesmat.h"
 #include "rankinematnl.h"
 #include "rankinematgrad.h"
 #include "trabbonematerial.h"
 #include "trabbonenl.h"
 #include "trabbone3d.h"
 #include "trabboneembed.h"
 #include "trabbonenlembed.h"
 #include "trabbonenl3d.h"
 #include "concretedpm.h"
 #include "concretedpm2.h"
 #include "cohint.h"
 #include "latticedamage2d.h"
 #include "simpleinterfacemat.h"
// #include "abaqususermaterial.h"
#endif // __SM_MODULE

#ifdef __TM_MODULE
 #include "isoheatmat.h"
 #include "hemotkmat.h"
 #include "hydratingisoheatmat.h"
 #include "hydratinghemomat.h"
 #include "hydratingconcretemat.h"
 #include "anisomassmat.h"
 #include "nonlinearmassmat.h"
 #ifdef __CEMHYD_MODULE
  #include "cemhydmat.h"
 #endif
#endif //__TM_MODULE

#ifdef __FM_MODULE
 #include "newtonianfluid.h"
 #include "fe2fluidmaterial.h"
 #include "surfacetensionmaterial.h"
 #include "twofluidmaterial.h"
 #include "binghamfluid2.h"
#endif // __FM_Module

#ifdef __SM_MODULE
REGISTER_CLASS(DummyMaterial, "dummymat", DummyMaterialClass)
REGISTER_CLASS(IsotropicLinearElasticMaterial, "isole", IsotropicLinearElasticMaterialClass)
REGISTER_CLASS(OrthotropicLinearElasticMaterial, "orthole", OrthotropicLinearElasticMaterialClass)
REGISTER_CLASS(Steel1, "steel1", Steel1MaterialClass)
REGISTER_CLASS(Concrete2, "concrete2", Concrete2Class)
REGISTER_CLASS(Concrete3, "concrete3", Concrete3Class)
REGISTER_CLASS(CebFip78Material, "cebfip78", CebFip78MaterialClass)
REGISTER_CLASS(DoublePowerLawMaterial, "doublepowerlaw", DoublePowerLawMaterialClass)
REGISTER_CLASS(B3Material, "b3mat", B3MaterialClass)
REGISTER_CLASS(B3SolidMaterial, "b3solidmat", B3SolidMaterialClass)
REGISTER_CLASS(MPSMaterial, "mps", MPSMaterialClass)
REGISTER_CLASS(J2plasticMaterial, "j2mat", J2plasticMaterialClass)
REGISTER_CLASS(RCSDNLMaterial, "rcsdnl", RCSDNLMaterialClass)
REGISTER_CLASS(RCSDEMaterial, "rcsde", RCSDEMaterialClass)
REGISTER_CLASS(RCSDMaterial, "rcsd", RCSDMaterialClass)
REGISTER_CLASS(M4Material, "microplane_m4", M4MaterialClass)
REGISTER_CLASS(IsotropicDamageMaterial1, "idm1", IsotropicDamageMaterial1Class)
REGISTER_CLASS(IDNLMaterial, "idmnl1", IDNLMaterialClass)
REGISTER_CLASS(MazarsNLMaterial, "mazarsmodelnl", MazarsNLMaterialClass)
REGISTER_CLASS(MazarsMaterial, "mazarsmodel", MazarsMaterialClass)
REGISTER_CLASS(DruckerPragerPlasticitySM, "druckerprager", DruckerPragerPlasticitySMClass)
REGISTER_CLASS(DruckerPragerMat, "druckerpragermat", DruckerPragerMatClass)
REGISTER_CLASS(J2MPlasticMaterial, "j2mmat", J2MPlasticMaterialClass)
REGISTER_CLASS(RankinePlasticMaterial, "rankine", RankinePlasticMaterialClass)
REGISTER_CLASS(Masonry02, "masonry02", Masonry02Class)
REGISTER_CLASS(IsoInterfaceDamageMaterial, "isointrfdm01", IsoInterfaceDamageMaterialClass)
REGISTER_CLASS(J2Mat, "j22mat", J2MatClass)
REGISTER_CLASS(CebFipSlip90Material, "cebfipslip90", CebFipSlip90MaterialClass)
REGISTER_CLASS(HellmichMaterial, "hellmat", HellmichMaterialClass)
REGISTER_CLASS(MDM, "mdm", MDMClass)
REGISTER_CLASS(CompoDamageMat, "compdammat", CompoDamageMatClass)
REGISTER_CLASS(MicroMaterial, "micromat", MicroMaterialClass)
REGISTER_CLASS(HyperElasticMaterial, "hyperelmat", HyperElasticMaterialClass)
REGISTER_CLASS(MisesMatGrad, "misesmatgrad", MisesMatGradClass)
REGISTER_CLASS(MisesMatNl, "misesmatnl", MisesMatNlClass)
REGISTER_CLASS(MisesMat, "misesmat", MisesMatClass)
REGISTER_CLASS(RankineMatGrad, "rankmatgrad", RankineMatGradClass)
REGISTER_CLASS(RankineMatNl, "rankmatnl", RankineMatNlClass)
REGISTER_CLASS(RankineMat, "rankmat", RankineMatClass)
REGISTER_CLASS(TrabBoneNL3D, "trabbonenl3d", TrabBoneNL3DClass)
REGISTER_CLASS(TrabBoneEmbed, "trabboneembed", TrabBoneEmbedClass)
REGISTER_CLASS(TrabBoneNLEmbed, "trabbonenlembed", TrabBoneNLEmbedClass)
REGISTER_CLASS(TrabBoneNL, "trabbonenl", TrabBoneNLClass)
REGISTER_CLASS(TrabBone3D, "trabbone3d", TrabBone3DClass)
REGISTER_CLASS(TrabBoneMaterial, "trabbone", TrabBoneMaterialClass)
REGISTER_CLASS(ConcreteDPM, "concretedpm", ConcreteDPMClass)
REGISTER_CLASS(ConcreteDPM, "concreteidm", ConcreteDPMClass) // for compatibility with old inputfiles
REGISTER_CLASS(CohesiveInterfaceMaterial, "cohint", CohesiveInterfaceMaterialClass)
REGISTER_CLASS(SimpleInterfaceMaterial, "simpleintermat", SimpleInterfaceMaterialClass)
REGISTER_CLASS(ConcreteDPM2, "con2dpm", ConcreteDPM2Class)
REGISTER_CLASS(LatticeDamage2d, "latticedamage2d", LatticeDamage2dClass)
//REGISTER_CLASS(AbaqusUserMaterial, "abaqususermaterial", AbaqusUserMaterialClass)
#endif //__SM_MODULE
#ifdef __TM_MODULE
REGISTER_CLASS(IsotropicHeatTransferMaterial, "isoheat", IsotropicHeatTransferMaterialClass)
REGISTER_CLASS(HeMoTKMaterial, "hemotk", HeMoTKMaterialClass)
REGISTER_CLASS(HydratingConcreteMat, "hydratingconcretemat", HydratingConcreteMatClass)
REGISTER_CLASS(AnisotropicMassTransferMaterial, "anisomass", AnisotropicMassTransferMaterialClass)
REGISTER_CLASS(NonlinearMassTransferMaterial, "nonlinmass", NonlinearMassTransferMaterialClass)
 #ifdef __CEMHYD_MODULE
REGISTER_CLASS(CemhydMat, "cemhydmat", CemhydMatClass)
 #endif //__CEMHYD_MODULE
#endif //__TM_MODULE
#if defined ( __SM_MODULE ) && defined ( __TM_MODULE )
REGISTER_CLASS(HydratingIsoHeatMaterial, "hisoheat", HydratingIsoHeatMaterialClass)
REGISTER_CLASS(HydratingHeMoMaterial, "hhemotk", HydratingHeMoMaterialClass)
#endif //defined(__SM_MODULE) && defined (__TM_MODULE)
#ifdef __FM_MODULE
REGISTER_CLASS(NewtonianFluidMaterial, "newtonianfluid", NewtonianFluidMaterialClass)
REGISTER_CLASS(TwoFluidMaterial, "twofluidmat", TwoFluidMaterialClass)
REGISTER_CLASS(BinghamFluidMaterial2, "binghamfluid2", BinghamFluidMaterial2Class)
REGISTER_CLASS(BinghamFluidMaterial2, "binghamfluid", BinghamFluidMaterial2Class)
REGISTER_CLASS(FE2FluidMaterial, "fe2fluidmaterial", FE2FluidMaterialClass)
REGISTER_CLASS(SurfaceTensionMaterial, "surfacetension", SurfaceTensionMaterialClass)
#endif // __FM_MODULE

