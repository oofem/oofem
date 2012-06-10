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

#include "springelement.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

SpringElement :: SpringElement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    springConstant  = 0.0;
}

void
SpringElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    /* spring stiffness matrix in local coordinate system (along orientation axis) */
    answer.resize(2,2);
    answer.at(1,1)=answer.at(2,2)= this->springConstant;
    answer.at(1,2)=answer.at(2,1)=-this->springConstant;
}


void
SpringElement :: giveInternalForcesVector(FloatArray &answer, TimeStep * atTime, int useUpdatedGpRecord)
{
    double f = this->computeSpringInternalForce(atTime);
    answer.resize(2);
    answer.at(1) = -f;
    answer.at(2) = f;
}


bool
SpringElement :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    /*
     * Spring is defined as 1D element along orientation axis (or around orientation axis for torsional springs)
     * The transformation from local (1d) to global system typically expand dimensions
     */
    if ((this->mode == SE_1D_SPRING) || (this->mode == SE_2D_TORSIONALSPRING_XZ)) {
      answer.resize(2,2);
      answer.at(1,1)=answer.at(2,2) = 1.0;
    } else if (this->mode == SE_2D_SPRING_XY) {
      answer.resize(2,4);
      answer.at(1,1)=this->dir.at(1);
      answer.at(1,2)=this->dir.at(2);
      answer.at(2,3)=this->dir.at(1);
      answer.at(2,4)=this->dir.at(2);
    } else if ((this->mode == SE_3D_SPRING) || (this->mode == SE_3D_TORSIONALSPRING)) {
      answer.resize(2,6);
      answer.at(1,1)=this->dir.at(1);
      answer.at(1,2)=this->dir.at(2);
      answer.at(1,3)=this->dir.at(3);
      answer.at(2,4)=this->dir.at(1);
      answer.at(2,5)=this->dir.at(2);
      answer.at(2,6)=this->dir.at(3);
    }
    return 1;
}


void
SpringElement :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    if (this->mode == SE_1D_SPRING) {
        answer.setValues(1, D_u);
    } else if (this->mode == SE_2D_SPRING_XY) {
        answer.setValues(2, D_u, D_v);
    } else if (this->mode == SE_2D_TORSIONALSPRING_XZ) {
        answer.setValues(1, R_v);
    } else if (this->mode == SE_3D_SPRING) {
        answer.setValues(3, D_u, D_v, D_w);
    } else if (this->mode == SE_3D_TORSIONALSPRING) {
        answer.setValues(3, R_u, R_v, R_w);
    }
}

double
SpringElement :: computeSpringInternalForce(TimeStep *stepN)
{
    FloatArray u;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    return (this->springConstant * (u.at(2)-u.at(1)));
}


int
SpringElement :: computeNumberOfGlobalDofs(EquationID eid)
{
    if ((this->mode == SE_1D_SPRING) || (this->mode == SE_2D_TORSIONALSPRING_XZ))  {
        return 2;
    } else if (this->mode == SE_2D_SPRING_XY) {
        return 4;
    } else if ((this->mode == SE_3D_SPRING) || (this->mode == SE_3D_TORSIONALSPRING)) {
        return 6;
    }
    return 0;
}


IRResultType
SpringElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // first call parent
    StructuralElement :: initializeFrom(ir);

    int _mode;
    IR_GIVE_FIELD(ir, _mode, IFT_SpringElement_mode, "mode"); // Macro
    IR_GIVE_FIELD(ir, springConstant, IFT_SpringElement_springConstant, "k"); // Macro

    this->mode = (SpringElementType) _mode;
    if (mode != SE_1D_SPRING) {
      IR_GIVE_OPTIONAL_FIELD(ir, this->dir, IFT_SpringElement_orientation, "orientation"); // Macro
      this->dir.normalize();
    }
    return IRRT_OK;
}

void SpringElement :: printOutputAt(FILE *File, TimeStep *stepN)
{
    fprintf(File, "spring element %d :\n", number);
    fprintf(File, "  spring force or moment % .4e", this->computeSpringInternalForce(stepN));
    fprintf(File, "\n");
}


} // end namespace oofem
