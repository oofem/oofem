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

#ifdef __SM_MODULE

// Elements of SM module
 #include "truss2d.h"
 #include "trplanstrss.h"
 #include "trplanrot.h"
 #include "trplanrot3d.h"
 #include "libeam2d.h"
 #include "planstrss.h"
 #include "quad1planestrain.h"
 #include "qplanestrain.h"
 #include "qplanestraingrad.h"
 #include "qtrplanestraingrad.h"
 #include "qtrplanestrain.h"
 #include "qtruss1d.h"
 #include "misesmatgrad.h"
 #include "misesmatnl.h"
 #include "qtruss1dgrad.h"
 #include "qspacegrad.h"
 #include "qplanstrss.h"
 #include "qplanestressgrad.h"
 #include "qtrplstr.h"
 #include "qtrplstrgrad.h"
 #include "lspace.h"
 #include "lspacebb.h"
 #include "qspace.h"
 #include "axisymm3d.h"
 #include "q4axisymm.h"
 #include "l4axisymm.h"
 #include "ltrspace.h"
 #include "beam2d.h"
 #include "beam3d.h"
 #include "libeam2dnl.h"
 #include "libeam3dnl.h"
 #include "truss3d.h"
 #include "trplanestrain.h"
 #include "libeam3dnl2.h"
 #include "libeam3d.h"
 #include "libeam3d2.h"
 #include "truss1d.h"
 #include "cct.h"
 #include "cct3d.h"
 #include "tr_shell01.h"
 #include "rershell.h"
 #include "interfaceelem2dquad.h"
 #include "interfaceelement1d.h"
 #include "interfaceelem3dtrlin.h"
 #include "macrolspace.h"
 #include "planstrssxfem.h"
 #include "cohsur3d.h"
 #include "lumpedmasselement.h"
 #include "springelement.h"
 #include "particle.h"
 #include "lattice2d.h"

// iga elements
 #include "igaelements.h"
#endif

#ifdef __TM_MODULE
// Elements of TM module
 #include "quad1_ht.h"
 #include "tr1_ht.h"
 #include "quadaxisym1_ht.h"
 #include "traxisym1_ht.h"
 #include "brick1_ht.h"
 #include "tetrah1_ht.h"
#endif

#ifdef __FM_MODULE
// Elements
 #include "tr1_2d_cbs.h"
 #include "tr21stokes.h"
 #include "tet21stokes.h"
 #include "linesurfacetension.h"
 #include "line2surfacetension.h"
 #include "line2boundaryelement.h"
 #include "tr1_2d_supg_axi.h"
 #include "tr1_2d_supg2_axi.h"
 #include "tr1_2d_supg2.h"
 #include "tr1_2d_supg.h"
 #include "quad10_2d_supg.h"
 #include "tet1_3d_supg.h"
 #include "tr21_2d_supg.h"
#endif

#ifdef __SM_MODULE
REGISTER_CLASS(PlaneStress2dXfem, "planestress2dxfem", PlaneStress2dXfemClass)
REGISTER_CLASS(PlaneStress2d, "planestress2d", PlaneStress2dClass)
REGISTER_CLASS(Quad1PlaneStrain, "quad1planestrain", Quad1PlaneStrainClass)
REGISTER_CLASS(QTrPlaneStrainGrad, "qtrplanestraingrad", QTrPlaneStrainGradClass)
REGISTER_CLASS(QTrPlaneStrain, "qtrplanestrain", QTrPlaneStrainClass)
REGISTER_CLASS(QPlaneStrainGrad, "qplanestraingrad", QPlaneStrainGradClass)
REGISTER_CLASS(QPlaneStrain, "qplanestrain", QPlaneStrainClass)
REGISTER_CLASS(TrPlaneStress2d, "trplanestress2d", TrPlaneStress2dClass)
REGISTER_CLASS(TrPlaneStrRot3d, "trplanestrrot3d", TrPlaneStrRot3dClass)
REGISTER_CLASS(TrPlaneStrRot, "trplanestrrot", TrPlaneStrRotClass)
REGISTER_CLASS(QPlaneStressGrad, "qplanestressgrad", QPlaneStressGradClass)
REGISTER_CLASS(QPlaneStress2d, "qplanestress2d", QPlaneStress2dClass)
REGISTER_CLASS(QTrPlaneStressGrad, "qtrplstrgrad", QTrPlaneStressGradClass)
REGISTER_CLASS(QTrPlaneStress2d, "qtrplstr", QTrPlaneStress2dClass)
REGISTER_CLASS(Axisymm3d, "axisymm3d", Axisymm3dClass)
REGISTER_CLASS(Q4Axisymm, "q4axisymm", Q4AxisymmClass)
REGISTER_CLASS(L4Axisymm, "l4axisymm", L4AxisymmClass)
REGISTER_CLASS(LSpaceBB, "lspacebb", LSpaceBBClass)
REGISTER_CLASS(LSpace, "lspace", LSpaceClass)
REGISTER_CLASS(QSpaceGrad, "qspacegrad", QSpaceGradClass)
REGISTER_CLASS(QSpace, "qspace", QSpaceClass)
REGISTER_CLASS(CCTPlate3d, "cctplate3d", CCTPlate3dClass)
REGISTER_CLASS(CCTPlate, "cctplate", CCTPlateClass)
//REGISTER_CLASS(LTRSpaceWithEmbeddedCrack, "ltrspaceec", )
REGISTER_CLASS(LTRSpace, "ltrspace", LTRSpaceClass)
REGISTER_CLASS(Truss2d, "truss2d", Truss2dClass)
REGISTER_CLASS(RerShell, "rershell", RerShellClass)
REGISTER_CLASS(TR_SHELL01, "tr_shell01", TR_SHELL01Class)
REGISTER_CLASS(Beam2d, "beam2d", Beam2dClass)
REGISTER_CLASS(Beam3d, "beam3d", Beam3dClass)
REGISTER_CLASS(LIBeam2dNL, "libeam2dNL", LIBeam2dNLClass)
REGISTER_CLASS(LIBeam2d, "libeam2d", LIBeam2dClass)
REGISTER_CLASS(LIBeam3dNL2, "libeam3dnl2", LIBeam3dNL2Class)
REGISTER_CLASS(LIBeam3dNL, "libeam3dNL", LIBeam3dNLClass)
REGISTER_CLASS(Truss3d, "truss3d", Truss3dClass)
REGISTER_CLASS(TrPlaneStrain, "trplanestrain", TrPlaneStrainClass)
REGISTER_CLASS(LIBeam3d2, "libeam3d2", LIBeam3d2Class)
REGISTER_CLASS(LIBeam3d, "libeam3d", LIBeam3dClass)
REGISTER_CLASS(QTruss1dGrad, "qtruss1dgrad", QTruss1dGradClass)
REGISTER_CLASS(QTruss1d, "qtruss1d", QTruss1dClass)
REGISTER_CLASS(Truss1d, "truss1d", Truss1dClass)
REGISTER_CLASS(InterfaceElem2dQuad, "interface2dquad", InterfaceElem2dQuadClass)
REGISTER_CLASS(InterfaceElement3dTrLin, "interface3dtrlin", InterfaceElement3dTrLinClass)
REGISTER_CLASS(InterfaceElem1d, "interface1d", InterfaceElem1dClass)
REGISTER_CLASS(MacroLSpace, "macrolspace", MacroLSpaceClass)
REGISTER_CLASS(LumpedMassElement, "lumpedmass", LumpedMassElementClass)
REGISTER_CLASS(SpringElement, "spring", SpringElementClass)
REGISTER_CLASS(CohesiveSurface3d, "cohsur3d", CohesiveSurface3dClass)
REGISTER_CLASS(BsplinePlaneStressElement, "bsplineplanestresselement", BsplinePlaneStressElementClass)
REGISTER_CLASS(NURBSPlaneStressElement, "nurbsplanestresselement", NURBSPlaneStressElementClass)
REGISTER_CLASS(TSplinePlaneStressElement, "tsplineplanestresselement", TSplinePlaneStressElementClass)
REGISTER_CLASS(NURBSSpace3dElement, "nurbs3delement", NURBSSpace3dElementClass)
REGISTER_CLASS(Lattice2d, "lattice2d", Lattice2dClass)
#endif

#ifdef __FM_MODULE
REGISTER_CLASS(TR1_2D_CBS, "tr1cbs", TR1_2D_CBSClass)
REGISTER_CLASS(TR1_2D_SUPG_AXI, "tr1supgaxi", TR1_2D_SUPG_AXIClass)
REGISTER_CLASS(TR1_2D_SUPG2_AXI, "tr1supg2axi", TR1_2D_SUPG2_AXIClass)
REGISTER_CLASS(TR1_2D_SUPG2, "tr1supg2", TR1_2D_SUPG2Class)
REGISTER_CLASS(TR1_2D_SUPG, "tr1supg", TR1_2D_SUPGClass)
REGISTER_CLASS(Tet1_3D_SUPG, "tet1supg", Tet1_3D_SUPGClass)
REGISTER_CLASS(TR21_2D_SUPG, "tr21supg", TR21_2D_SUPGClass)
REGISTER_CLASS(Quad10_2D_SUPG, "quad1supg", Quad10_2D_SUPGClass)
REGISTER_CLASS(Tr21Stokes, "tr21stokes", Tr21StokesElementClass)
REGISTER_CLASS(Tet21Stokes, "tet21stokes", Tet21StokesElementClass)
REGISTER_CLASS(Line2BoundaryElement, "line2boundary", Line2BoundaryElementClass)
REGISTER_CLASS(LineSurfaceTension, "linesurfacetension", LineSurfaceTensionElementClass)
REGISTER_CLASS(Line2SurfaceTension, "line2surfacetension", Line2SurfaceTensionElementClass)
#endif

#ifdef __TM_MODULE
REGISTER_CLASS(Quad1_ht, "quad1ht", Quad1_htClass)
REGISTER_CLASS(Quad1_hmt, "quad1hmt", Quad1_hmtClass)
REGISTER_CLASS(QuadAxisym1_ht, "quadaxisym1ht", QuadAxisym1_htClass)
REGISTER_CLASS(QuadAxisym1_hmt, "quadaxisym1hmt", QuadAxisym1_hmtClass)
REGISTER_CLASS(Brick1_ht, "brick1ht", Brick1_htClass)
REGISTER_CLASS(Brick1_hmt, "brick1hmt", Brick1_hmtClass)
REGISTER_CLASS(Tetrah1_ht, "tetrah1ht", Tetrah1_htClass)
REGISTER_CLASS(Tetrah1_hmt, "tetrah1hmt", Tetrah1_hmtClass)
REGISTER_CLASS(Tr1_ht, "tr1ht", Tr1_htClass)
REGISTER_CLASS(Tr1_hmt, "tr1hmt", Tr1_hmtClass)
REGISTER_CLASS(TrAxisym1_ht, "traxisym1ht", TrAxisym1_htClass)
#endif
