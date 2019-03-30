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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "Elements/GradientDamage/PlaneStrain/quad1planestraingraddamage.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "Materials/rcm2.h"
#endif

namespace oofem {
REGISTER_Element(Quad1PlaneStrainGradDamage);
  IntArray Quad1PlaneStrainGradDamage :: locationArray_u = {1, 2, 4, 5, 7, 8, 10, 11};
  IntArray Quad1PlaneStrainGradDamage :: locationArray_d = {3, 6, 9, 12};

Quad1PlaneStrainGradDamage :: Quad1PlaneStrainGradDamage(int n, Domain *aDomain) : Quad1PlaneStrain(n, aDomain), GradientDamageElement()
{
  nPrimNodes = 4;
    nPrimVars = 2;
    nSecNodes = 4;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
  
}





  

void
Quad1PlaneStrainGradDamage :: giveDofManDofIDMask(int inode, IntArray &answer) const

{
    answer = {D_u, D_v, G_0};  
}


void
Quad1PlaneStrainGradDamage :: giveDofManDofIDMask_u(IntArray &answer) const
{
  answer = {D_u, D_v};
}


void
Quad1PlaneStrainGradDamage :: giveDofManDofIDMask_d(IntArray &answer) const
{
  answer = {G_0};   
}

  
  

void
Quad1PlaneStrainGradDamage :: computeNdMatrixAt(GaussPoint *gp, FloatArray &answer)
{
    this->interp.evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}

void
Quad1PlaneStrainGradDamage :: computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;
    this->interp.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dnx);
}

  
void
Quad1PlaneStrainGradDamage :: giveLocationArray_u(IntArray &answer)
{
  answer = locationArray_u;
}

  
void
Quad1PlaneStrainGradDamage :: giveLocationArray_d(IntArray &answer)
{
  answer = locationArray_d;
}
  
void
Quad1PlaneStrainGradDamage :: postInitialize()
{
  GradientDamageElement:: postInitialize();
  NLStructuralElement :: postInitialize();
}


  

} // end namespace oofem


