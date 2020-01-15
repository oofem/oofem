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

#include "sm/Elements/nodalspringelement.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(NodalSpringElement);

NodalSpringElement :: NodalSpringElement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 1;
}

void
NodalSpringElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    answer.beDiagonal(this->springConstants);
}


void
NodalSpringElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  int ndofs = computeNumberOfDofs();
  FloatArray u;
  this->computeVectorOf(VM_Total, tStep, u);

  answer.resize(ndofs);
  
  for (int i = 1; i<=ndofs; i++) {
    answer.at(i) = this->springConstants.at(i)*u.at(i);
  }
}


bool
NodalSpringElement :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
  return false;
}


void
NodalSpringElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer= this->dofMask;
}

void 
NodalSpringElement::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
  if (this->masses.isNotEmpty()) {
    answer.beDiagonal(masses);
  } else {
    int ndofs = this->computeNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
  }
}

int
NodalSpringElement :: computeNumberOfGlobalDofs()
{
  return this->computeNumberOfDofs();
}


void
NodalSpringElement :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, dofMask, _IFT_NodalSpringElement_dofmask);
    IR_GIVE_FIELD(ir, springConstants, _IFT_NodalSpringElement_springConstants);

    this->masses.resize(dofMask.giveSize());
    this->masses.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, masses, _IFT_NodalSpringElement_masses);
    
    int ndofs = dofMask.giveSize();
    if (ndofs != springConstants.giveSize()) {
      OOFEM_ERROR ("Spring constants size not equal to number of DOFs");
    }

    if (ndofs != masses.giveSize()) {
      OOFEM_ERROR ("Masses array size not equal to number of DOFs");
    }

    // from element
    activityTimeFunction = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, activityTimeFunction, _IFT_Element_activityTimeFunction);
    IR_GIVE_FIELD(ir, dofManArray, _IFT_Element_nodes);
}

void NodalSpringElement :: printOutputAt(FILE *File, TimeStep *tStep)
{
  FloatArray F;
  this->giveInternalForcesVector(F, tStep);
  
  fprintf(File, "NodalSpring element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
  fprintf(File, "  Force/Moment ");
  for (int i = 1; i<= F.giveSize(); i++) {
    fprintf(File, "%.4e ", F.at(i));
  }
  fprintf(File, "\n");
}
} // end namespace oofem
