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

#include "mat_cebfip90.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(CebFipSlip90Material);

CebFipSlip90Material :: CebFipSlip90Material(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{ }


double
CebFipSlip90Material :: giveEngTraction_1d(double jump, GaussPoint *gp, TimeStep *tStep) const
{
    CebFipSlip90MaterialStatus *status = static_cast< CebFipSlip90MaterialStatus * >( this->giveStatus(gp) );
    double slip = jump; ///@todo This isn't slip, this is the normal displacement. This code should be rewritten to deal with a 3D jump.
    // compute value of loading function if strainLevel crit apply
    double f = fabs(slip) - status->giveKappa();

    double tempKappa;
    if ( f <= 0.0 ) {
        // kappa do not grow
        tempKappa = status->giveKappa();
    } else {
        // kappa grow
        tempKappa = fabs(slip);
        // evaluate damage parameter
    }

    double answer = 0.0;
    if ( tempKappa > 1.e-12 ) {
        answer = ( this->computeBondForce(tempKappa) / tempKappa ) * slip;
    }

    // update gp
    status->letTempJumpBe({jump, 0., 0.});
    status->letTempTractionBe({answer, 0., 0.});
    status->setTempKappa(tempKappa);

    return answer;
}


FloatMatrixF<1,1>
CebFipSlip90Material :: give1dStiffnessMatrix_Eng(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    CebFipSlip90MaterialStatus *status = static_cast< CebFipSlip90MaterialStatus * >( this->giveStatus(gp) );
    FloatMatrixF<1,1> answer;

    if ( mode == ElasticStiffness || mode == SecantStiffness ) {
        double kappa = status->giveKappa();

        if ( kappa > 0.0 ) {
            answer.at(1, 1) = ( this->computeBondForce(kappa) / kappa );
        } else {
            answer.at(1, 1) = computeBondForceStiffness(0.0);
        }
    } else if ( mode == TangentStiffness ) {
        answer.at(1, 1) = computeBondForceStiffness( status->giveTempKappa() );
    }  else {
        OOFEM_ERROR("unknown MatResponseMode (%s)", __MatResponseModeToString(mode) );
    }
    return answer;
}


int
CebFipSlip90Material :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    CebFipSlip90MaterialStatus *status = static_cast< CebFipSlip90MaterialStatus * >( this->giveStatus(gp) );

    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


void
CebFipSlip90Material :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, tmax, _IFT_CebFipSlip90Material_tmax);
    IR_GIVE_FIELD(ir, tres, _IFT_CebFipSlip90Material_tres);

    IR_GIVE_FIELD(ir, s1, _IFT_CebFipSlip90Material_s1);
    IR_GIVE_FIELD(ir, s2, _IFT_CebFipSlip90Material_s2);
    IR_GIVE_FIELD(ir, s3, _IFT_CebFipSlip90Material_s3);

    alpha = 0.4;
}


void
CebFipSlip90Material :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);
    input.setField(this->tmax, _IFT_CebFipSlip90Material_tmax);
    input.setField(this->tres, _IFT_CebFipSlip90Material_tres);

    input.setField(this->s1, _IFT_CebFipSlip90Material_s1);
    input.setField(this->s2, _IFT_CebFipSlip90Material_s2);
    input.setField(this->s3, _IFT_CebFipSlip90Material_s3);
}



double
CebFipSlip90Material :: computeBondForce(double s) const
{
    if ( s <= s1 ) {
        return tmax *pow( ( s / s1 ), alpha );
    } else if ( ( s >= s1 ) && ( s <= s2 ) ) {
        return tmax;
    } else if ( ( s >= s2 ) && ( s <= s3 ) ) {
        return tmax - ( tmax - tres ) * ( s - s1 ) / ( s3 - s2 );
    } else {
        return tres;
    }
}


double
CebFipSlip90Material :: computeBondForceStiffness(double s) const
{
    if ( s <= s1 / 1000. ) {
        s = s1 / 1000.;
        return alpha *tmax *pow( ( s / s1 ), alpha - 1 ) / s1;
    } else if ( s <= s1 ) {
        return alpha *tmax *pow( ( s / s1 ), alpha - 1 ) / s1;
    } else if ( ( s >= s1 ) && ( s <= s2 ) ) {
        return 1.e-6; // should be zero
    } else if ( ( s >= s2 ) && ( s <= s3 ) ) {
        return -( tmax - tres ) / ( s3 - s2 );
    } else {
        return 1.e-6; // should be zero
    }
}


CebFipSlip90MaterialStatus :: CebFipSlip90MaterialStatus(GaussPoint *g) : StructuralInterfaceMaterialStatus(g)
{}


void
CebFipSlip90MaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "kappa %f", this->kappa);

    fprintf(file, "}\n");
}


void
CebFipSlip90MaterialStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
}


void
CebFipSlip90MaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
}


void
CebFipSlip90MaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralInterfaceMaterialStatus :: saveContext(stream, mode);

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
CebFipSlip90MaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralInterfaceMaterialStatus :: restoreContext(stream, mode);

    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
} // end namespace oofem
