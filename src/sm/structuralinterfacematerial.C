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

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
namespace oofem {

int
StructuralInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_InterfaceJump ) {
        answer = status->giveJump();
        return 1;

    } else if ( type == IST_InterfaceTraction ) {
        answer = status->giveTraction();
        return 1;

    } else if ( type == IST_InterfaceFirstPKTraction ) {
        answer = status->giveFirstPKTraction();
        return 1;

    } else if ( type == IST_DeformationGradientTensor ) {
        answer.beVectorForm( status->giveF() );
        return 1;

    } else {
        return Material :: giveIPValue(answer, gp, type, atTime);
    }
}


IRResultType
StructuralInterfaceMaterial :: initializeFrom(InputRecord *ir)
{

    const char *__proc = "initializeFrom";  // Required by IR_GIVE_FIELD macro
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, this->useNumericalTangent, _IFT_StructuralInterfaceMaterial_useNumericalTangent);

    return IRRT_OK;
}


void
StructuralInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
}


#if 1
void
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_Num(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{

    // Numerical tangent
    // Computes the material stiffness using a central difference method
	
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );
    if ( status ) {
        FloatMatrix F;
        F = status->giveTempF();
        int dim = F.giveNumberOfRows();
        answer.resize(dim,dim);
        answer.zero();
        const double eps = 1.0e-9;
        FloatArray T, TPlus, TMinus;
        
        
        FloatArray jump, jumpPlus, jumpMinus, Kcolumn;
        jump = status->giveTempJump();

        for ( int i = 1; i <= dim; i++ ) {
            jumpPlus = jumpMinus = jump;
            jumpPlus.at(i)  += eps; 
            jumpMinus.at(i) -= eps;
            this->giveFirstPKTraction_3d(TPlus, gp, jumpPlus, F, tStep);
            this->giveFirstPKTraction_3d(TMinus, gp, jumpMinus, F, tStep);

            Kcolumn = (TPlus - TMinus);
            answer.setColumn(Kcolumn, i);
        }
        answer.times( 1.0/(2*eps) );

        this->giveFirstPKTraction_3d(T, gp, jump, F, tStep); // reset temp values by recomputing the stress

    }

}

#endif

} // end namespace oofem
