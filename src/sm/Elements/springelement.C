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

#include "../sm/Elements/springelement.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(SpringElement);

SpringElement :: SpringElement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    springConstant  = 0.0;
}

void
SpringElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    /* spring stiffness matrix in local coordinate system (along orientation axis) */
    answer.resize(2, 2);
    answer.at(1, 1) = answer.at(2, 2) = this->springConstant;
    answer.at(1, 2) = answer.at(2, 1) = -this->springConstant;
}


void
SpringElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    double f = this->computeSpringInternalForce(tStep);
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
    if ( ( this->mode == SE_1D_SPRING ) || ( this->mode == SE_2D_TORSIONALSPRING_XZ ) ) {
        answer.resize(2, 2);
        answer.at(1, 1) = answer.at(2, 2) = 1.0;
    } else if ( this->mode == SE_2D_SPRING_XY ) {
        answer.resize(2, 4);
        answer.at(1, 1) = this->dir.at(1);
        answer.at(1, 2) = this->dir.at(2);
        answer.at(2, 3) = this->dir.at(1);
        answer.at(2, 4) = this->dir.at(2);
    } else if ( this->mode == SE_2D_SPRING_XZ ) {
        answer.resize(2, 4);
        answer.at(1, 1) = this->dir.at(1);
        answer.at(1, 2) = this->dir.at(3);
        answer.at(2, 3) = this->dir.at(1);
        answer.at(2, 4) = this->dir.at(3);
    } else if ( ( this->mode == SE_3D_SPRING ) || ( this->mode == SE_3D_TORSIONALSPRING ) ) {
        answer.resize(2, 6);
        answer.at(1, 1) = this->dir.at(1);
        answer.at(1, 2) = this->dir.at(2);
        answer.at(1, 3) = this->dir.at(3);
        answer.at(2, 4) = this->dir.at(1);
        answer.at(2, 5) = this->dir.at(2);
        answer.at(2, 6) = this->dir.at(3);
    }
    return 1;
}


void
SpringElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( this->mode == SE_1D_SPRING ) {
        answer = {D_u};
    } else if ( this->mode == SE_2D_SPRING_XY ) {
       answer = {D_u, D_v};
    } else if ( this->mode == SE_2D_SPRING_XZ ) {
       answer = {D_u, D_w};
    } else if ( this->mode == SE_2D_TORSIONALSPRING_XZ ) {
        answer = {R_v};
    } else if ( this->mode == SE_3D_SPRING ) {
        answer = {D_u, D_v, D_w};
    } else if ( this->mode == SE_3D_TORSIONALSPRING ) {
        answer = {R_u, R_v, R_w};
    }
}

double
SpringElement :: computeSpringInternalForce(TimeStep *tStep)
{
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
    return ( this->springConstant * ( u.at(2) - u.at(1) ) );
}

void 
SpringElement::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
  answer.resize(2, 2);
  answer.at(1,1)=answer.at(2,2) = this->mass/2.0;
  answer.at(1,2)=answer.at(2,1) = 0.0;
}

int
SpringElement :: computeNumberOfGlobalDofs()
{
    if ( ( this->mode == SE_1D_SPRING ) || ( this->mode == SE_2D_TORSIONALSPRING_XZ ) ) {
        return 2;
    } else if ( ( this->mode == SE_2D_SPRING_XY ) || ( this->mode == SE_2D_SPRING_XZ ) ) {
        return 4;
    } else if ( ( this->mode == SE_3D_SPRING ) || ( this->mode == SE_3D_TORSIONALSPRING ) ) {
        return 6;
    }
    return 0;
}


IRResultType
SpringElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int _mode;
    IR_GIVE_FIELD(ir, _mode, _IFT_SpringElement_mode);
    IR_GIVE_FIELD(ir, springConstant, _IFT_SpringElement_springConstant);
    this->mass = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->mass, _IFT_SpringElement_mass);

    this->mode = ( SpringElementType ) _mode;
    if ( mode != SE_1D_SPRING ) {
        IR_GIVE_OPTIONAL_FIELD(ir, this->dir, _IFT_SpringElement_orientation);
        this->dir.normalize();
    }
    return StructuralElement :: initializeFrom(ir);
}

void SpringElement :: printOutputAt(FILE *File, TimeStep *tStep)
{
    fprintf(File, "spring element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
    fprintf(File, "  spring force or moment %.4e", this->computeSpringInternalForce(tStep) );
    fprintf(File, "\n");
}
} // end namespace oofem
