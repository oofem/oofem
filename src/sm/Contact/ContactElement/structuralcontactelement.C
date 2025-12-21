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
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "structuralcontactelement.h"
#include "feinterpol.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "integrationrule.h"
#include "crosssection.h"
#include "gaussintegrationrule.h"
#include "lobattoir.h"
//@todo: add lobatto ir
namespace oofem {

REGISTER_Element(StructuralContactElement_LineLin);
REGISTER_Element(StructuralContactElement_TrLin);
REGISTER_Element(StructuralContactElement_QuadLin);


StructuralContactElement :: StructuralContactElement(int n, Domain *aDomain) : ContactElement(n, aDomain)
{

}


void
StructuralContactElement :: computeNmatrixAt(const FloatArray &lcoord, FloatMatrix &answer)
{
    FloatArray n;
    this->giveInterpolation()->evalN( n, lcoord,  FEIElementGeometryWrapper(this));
    answer.beNMatrixOf(n, nsd);
}     



void StructuralContactElement :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.XF
{
  if ( integrationRulesArray.size() == 0 ) {
    integrationRulesArray.resize(1);
    //integrationRulesArray [ 0 ] = std::make_unique< GaussIntegrationRule >(1, this);
    integrationRulesArray [ 0 ] = std::make_unique<LobattoIntegrationRule>(1, this);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
  }
}


  

FEI3dQuadLin StructuralContactElement_QuadLin :: interpolation;

StructuralContactElement_QuadLin :: StructuralContactElement_QuadLin(int n, Domain *aDomain) : StructuralContactElement(n, aDomain)
{
  this->numberOfDofMans = 4;
  this->nsd = 3;
  this->numberOfGaussPoints = 4;
  
}
  
FEInterpolation *
StructuralContactElement_QuadLin :: giveInterpolation() const { return & interpolation; }


FloatArray
StructuralContactElement_QuadLin :: computeNormalVectorAt(const FloatArray &lCoords)
{
  FloatArray normal;
  auto norm = this->giveInterpolation()->boundarySurfaceEvalNormal(normal, 0, lCoords, FEIElementGeometryWrapper(this));
  return normal*norm;
  }
  

void
StructuralContactElement_QuadLin :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {
    D_u, D_v, D_w
  };
}




FEI3dTrLin StructuralContactElement_TrLin :: interpolation;

StructuralContactElement_TrLin :: StructuralContactElement_TrLin(int n, Domain *aDomain) : StructuralContactElement(n, aDomain)
{
  this->numberOfDofMans = 3;
  this->nsd = 3;
  this->numberOfGaussPoints = 3;
  
}
  
FEInterpolation *
StructuralContactElement_TrLin :: giveInterpolation() const { return & interpolation; }


FloatArray
StructuralContactElement_TrLin :: computeNormalVectorAt(const FloatArray &lCoords)
{
  FloatArray normal;
  auto norm = this->giveInterpolation()->boundarySurfaceEvalNormal(normal, 0, lCoords, FEIElementGeometryWrapper(this));
  return normal*norm;
  }
  

void
StructuralContactElement_TrLin :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {
    D_u, D_v, D_w
  };
}

  





FEI2dLineLin StructuralContactElement_LineLin :: interpolation(1, 2);

StructuralContactElement_LineLin :: StructuralContactElement_LineLin(int n, Domain *aDomain) : StructuralContactElement(n, aDomain)
{
  this->numberOfDofMans = 2;
  this->nsd = 2;
  this->numberOfGaussPoints = 2;
}
  
FEInterpolation *
StructuralContactElement_LineLin :: giveInterpolation() const { return & interpolation; }


FloatArray
StructuralContactElement_LineLin :: computeNormalVectorAt(const FloatArray &lCoords)
{
  FloatArray normal;
  auto norm = this->giveInterpolation()->boundaryEvalNormal(normal, 0, lCoords, FEIElementGeometryWrapper(this));
  return normal*norm;
}
  

void
StructuralContactElement_LineLin :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {
    D_u, D_v
  };
}


  

  

} // end namespace oofem
