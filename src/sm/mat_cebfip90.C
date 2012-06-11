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

#include "mat_cebfip90.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
CebFipSlip90Material :: CebFipSlip90Material(int n, Domain *d) : StructuralMaterial(n, d)
//
// constructor
//
{ }


CebFipSlip90Material :: ~CebFipSlip90Material()
//
// destructor
//
{ }

int
CebFipSlip90Material :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( mode == _1dInterface ) {
        return 1;
    }

    return 0;
}


void
CebFipSlip90Material :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseForm form,
                                                      MatResponseMode mode,
                                                      GaussPoint *gp,
                                                      TimeStep *atTime)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
    _error("give3dMaterialStiffnessMatrix: not implemented");
}


void
CebFipSlip90Material :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                             const FloatArray &totalStrain,
                                             TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    CebFipSlip90MaterialStatus *status = ( CebFipSlip90MaterialStatus * ) this->giveStatus(gp);
    FloatArray reducedTotalStrainVector;
    double f, slip, tempKappa;

    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);

    //crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);
    slip = reducedTotalStrainVector.at(1);
    // compute value of loading function if strainLevel crit apply
    f = fabs(slip) - status->giveKappa();

    if ( f <= 0.0 ) {
        // kappa do not grow
        tempKappa = status->giveKappa();
    } else {
        // kappa grow
        tempKappa = fabs(slip);
        // evaluate damage parameter
    }

    answer.resize(1);
    if ( tempKappa <= 1.e-12 ) {
        answer.at(1) = 0.0;
    } else {
        answer.at(1) = ( this->computeBondForce(tempKappa) / tempKappa ) * slip;
    }

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
}


void
CebFipSlip90Material :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                 MatResponseForm form, MatResponseMode rMode,
                                                 GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dInterface:
        give1dInterfaceMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
        break;
    default:
        StructuralMaterial :: giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
    }
}


int
CebFipSlip90Material :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
//
// returns the size of reduced stress-strain vector
// according to mode given by gp.
//
{
    switch ( mode ) {
    case _1dInterface:
        return 1;

    default:
        return StructuralMaterial :: giveSizeOfReducedStressStrainVector(mode);
    }
}


int
CebFipSlip90Material :: giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind)
//
// this function returns index of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain component in Full or reduced
// stressStrainVector acording to stressStrain mode of given gp.
//
{
    //MaterialMode mode  = gp -> giveMaterialMode ();

    if ( mmode == _1dInterface ) {
        return ind;
    } else {
        return StructuralMaterial :: giveStressStrainComponentIndOf(form, mmode, ind);
    }
}


void
CebFipSlip90Material :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                             MaterialMode mmode) const
//
// this function returns mask of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain vector in full or
// reduced StressStrainVector
// acording to stressStrain mode of given gp.
//
//
// mask has size of reduced or full StressStrain Vector and  i-th component
// is index to full or reduced StressStrainVector where corresponding
// stressStrain resides.
//
// Reduced form is sub-vector (of stress or strain components),
// where components corresponding to imposed zero stress (plane stress,...)
// are not included. On the other hand, if zero strain component is imposed
// (Plane strain, ..) this condition must be taken into account in geometrical
// relations, and corresponding component is included in reduced vector.
//
{
    if ( mmode == _1dInterface ) {
        answer.resize(1);
        answer.at(1) = 1;
    } else {
        StructuralMaterial :: giveStressStrainMask(answer, form, mmode);
    }
}


void
CebFipSlip90Material :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                        const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode == _1dInterface ) {
        answer = charVector3d;
        return;
    } else {
        StructuralMaterial :: giveReducedCharacteristicVector(answer, gp, charVector3d);
    }
}


void
CebFipSlip90Material :: giveFullCharacteristicVector(FloatArray &answer,
                                                     GaussPoint *gp,
                                                     const FloatArray &strainVector)
//
// returns full 3d general strain vector from strainVector in reducedMode
// based on StressStrainMode in gp. Included are strains which
// perform nonzero work.
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
// you must assigng your stress strain mode to one of the folloving modes (or add new)
// FullForm of MaterialStiffnessMatrix must have the same form.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _1dInterface ) {
        answer = strainVector;
        return;
    } else {
        StructuralMaterial :: giveFullCharacteristicVector(answer, gp, strainVector);
    }
}


void
CebFipSlip90Material :: give1dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                               GaussPoint *gp, TimeStep *atTime)
{
    double kappa;
    CebFipSlip90MaterialStatus *status = ( CebFipSlip90MaterialStatus * ) this->giveStatus(gp);
    answer.resize(1, 1);

    if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) ) {
        kappa = status->giveKappa();

        if ( kappa > 0.0 ) {
            answer.at(1, 1) = ( this->computeBondForce(kappa) / kappa );
        } else {
            answer.at(1, 1) = computeBondForceStiffness(0.0);
        }
    } else if ( rMode == TangentStiffness ) {
        answer.at(1, 1) = computeBondForceStiffness( status->giveTempKappa() );
    }  else {
        _error2( "give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
    }
}


int
CebFipSlip90Material :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
}


InternalStateValueType
CebFipSlip90Material :: giveIPValueType(InternalStateType type)
{
    return StructuralMaterial :: giveIPValueType(type);
}


int
CebFipSlip90Material :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
}


int
CebFipSlip90Material :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
}


void
CebFipSlip90Material :: giveThermalDilatationVector(FloatArray &answer,
                                                    GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    answer.resize(1);
    answer.at(1) = 0.0;
}


IRResultType
CebFipSlip90Material :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, tmax, IFT_CebFipSlip90Material_tmax, "tmax"); // Macro
    IR_GIVE_FIELD(ir, tres, IFT_CebFipSlip90Material_tres, "tres"); // Macro

    IR_GIVE_FIELD(ir, s1, IFT_CebFipSlip90Material_s1, "s1"); // Macro
    IR_GIVE_FIELD(ir, s2, IFT_CebFipSlip90Material_s2, "s2"); // Macro
    IR_GIVE_FIELD(ir, s3, IFT_CebFipSlip90Material_s3, "s3"); // Macro

    alpha = 0.4;
    return StructuralMaterial :: initializeFrom(ir);
}


int
CebFipSlip90Material :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    StructuralMaterial :: giveInputRecordString(str, keyword);

    sprintf(buff, " tmax %e tres %e s1 %e s2 %e s3 %e", tmax, tres, s1, s2, s3);
    str += buff;

    return 1;
}


double
CebFipSlip90Material :: computeBondForce(double s)
{
    if ( s <= s1 ) {
        return tmax * pow( ( s / s1 ), alpha );
    } else if ( ( s >= s1 ) && ( s <= s2 ) ) {
        return tmax;
    } else if ( ( s >= s2 ) && ( s <= s3 ) ) {
        return tmax - ( tmax - tres ) * ( s - s1 ) / ( s3 - s2 );
    } else {
        return tres;
    }
}


double
CebFipSlip90Material :: computeBondForceStiffness(double s)
{
    if ( s <= s1 / 1000. ) {
        s = s1 / 1000.;
        return alpha * tmax * pow( ( s / s1 ), alpha - 1 ) / s1;
    } else if ( s <= s1 ) {
        return alpha * tmax * pow( ( s / s1 ), alpha - 1 ) / s1;
    } else if ( ( s >= s1 ) && ( s <= s2 ) ) {
        return 1.e-6; // should be zero
    } else if ( ( s >= s2 ) && ( s <= s3 ) ) {
        return -( tmax - tres ) / ( s3 - s2 );
    } else {
        return 1.e-6; // should be zero
    }
}


CebFipSlip90MaterialStatus :: CebFipSlip90MaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
}


CebFipSlip90MaterialStatus :: ~CebFipSlip90MaterialStatus()
{ }


void
CebFipSlip90MaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "kappa %f", this->kappa);

    fprintf(file, "}\n");
}


void
CebFipSlip90MaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
}


void
CebFipSlip90MaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->kappa = this->tempKappa;
}


contextIOResultType
CebFipSlip90MaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
CebFipSlip90MaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
} // end namespace oofem
