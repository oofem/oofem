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

#include "rvestokesflow.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "stdio.h"
#include "rveengngmodel.h"
#include "engngm.h"
#include "nummet.h"
#include "unistd.h"
#include "contextioerr.h"
#include "gausspnt.h"
#include <string.h>
#include <sstream>
#include <math.h>

namespace oofem {

rvestokesflowMaterialStatus :: rvestokesflowMaterialStatus(int n, Domain *d, GaussPoint *g) :
    TransportMaterialStatus(n, d, g)
{
    temp_gradPVector.resize(2);
    temp_gradPVector.zero();

    temp_velocityVector.resize(2);
    temp_velocityVector.zero();

    temp_TangentMatrix.resize(2, 2);
    temp_TangentMatrix.zero();

    solutionVector = new FloatArray;
}

rvestokesflowMaterialStatus :: ~rvestokesflowMaterialStatus()
{}

void
rvestokesflowMaterialStatus :: exportFilter(EngngModel *E, GaussPoint *gp, TimeStep *tStep)
{
    // Fix name for RVE output file and print output

    std :: string filename, basefilename;
    std :: ostringstream filenameStringStream;


    filenameStringStream << this->giveDomain()->giveEngngModel()->giveOutputBaseFileName().c_str() << ".rve/Rve_" << gp->giveElement()->giveGlobalNumber() << "_" << gp->giveNumber();

    filename = filenameStringStream.str();

    // Update fields on subscale
    FloatArray grapP = this->giveTempGradP(), seepageVelocity;

    rveEngngModel *rveE;
    rveE = dynamic_cast< rveEngngModel * >(E);

    rveE->rveSetBoundaryConditions(10, grapP);
    rveE->rveGiveCharacteristicData(1, &grapP, &seepageVelocity, tStep);

    basefilename = E->giveOutputBaseFileName();

    E->letOutputBaseFileNameBe(filename);
    E->doStepOutput(tStep);
    E->letOutputBaseFileNameBe(basefilename);
}

void
rvestokesflowMaterialStatus :: initTempStatus()
{
    TransportMaterialStatus :: initTempStatus();

    temp_gradPVector = gradPVector;
    temp_velocityVector = velocityVector;
}

void
rvestokesflowMaterialStatus :: updateYourself(TimeStep *tStep)
{
    TransportMaterialStatus :: updateYourself(tStep);
    gradPVector = temp_gradPVector;
    velocityVector = temp_velocityVector;
    tangentMatrix = temp_TangentMatrix;

    EngngModel *E;

    rvestokesflow *Mtrl;
    Mtrl = dynamic_cast< rvestokesflow * >( gp->giveMaterial() );
    E = dynamic_cast< EngngModel * >(Mtrl->rve);

    Mtrl->suppressStdout();
    this->exportFilter(E, gp, tStep);
    Mtrl->enableStdout();
}


contextIOResultType
rvestokesflowMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = TransportMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
rvestokesflowMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = TransportMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

rvestokesflow :: rvestokesflow(int n, Domain *d) : RVEMaterial(n, d), TransportMaterial(n, d)
{

}

IRResultType rvestokesflow :: initializeFrom(InputRecord *ir)
{
    // this->TransportMaterial :: initializeFrom(ir);
    this->RVEMaterial :: initializeFrom(ir);

    return IRRT_OK;
}

int
rvestokesflow :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    rvestokesflowMaterialStatus *thisMaterialStatus;
    thisMaterialStatus = ( ( rvestokesflowMaterialStatus * ) this->giveStatus(aGaussPoint) );
    FloatMatrix temp;
    answer.resize(3);
    answer.zero();

    switch ( type ) {
    case IST_Velocity:
        answer.copySubVector(thisMaterialStatus->giveVelocityVector(), 1);
        break;
    case IST_PressureGradient:
        answer.copySubVector(thisMaterialStatus->giveGradP(), 1);
        break;
    case IST_TangentNorm:
        temp = thisMaterialStatus->giveTangentMatrix();
        answer.resize(1);
        answer.at(1) = temp.computeFrobeniusNorm();
        break;
    case IST_Tangent:
        temp = thisMaterialStatus->giveTangentMatrix();
        answer.resize(4);
        answer.at(1) = temp.at(1, 1);
        answer.at(2) = temp.at(1, 2);
        answer.at(3) = temp.at(2, 1);
        answer.at(4) = temp.at(2, 2);
        break;
    default:
        return TransportMaterial :: giveIPValueType(type);
    }

    return 1;
}

void
rvestokesflow :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    this->suppressStdout();

    OOFEM_LOG_DEBUG("\n****** Enter giveFluxVector ********************** Element number %u, Gauss point %u, rve @ %p \n", gp->giveElement()->giveGlobalNumber(), gp->giveNumber(), this->rve);

    rvestokesflowMaterialStatus *status = ( ( rvestokesflowMaterialStatus * ) this->giveStatus(gp) );
    FloatArray temp = status->giveTempGradP();

    if ( temp.giveSize() >= 3 && temp.at(1) == eps.at(1) && temp.at(2) == eps.at(2) ) {
    	OOFEM_LOG_DEBUG("This pressure gradient has already been evaluated\n");
        answer = status->giveTempVelocityVector();
    } else {
        FloatArray X;
        rveEngngModel *rveE;

        rveE = dynamic_cast< rveEngngModel * >(this->rve);

        X = eps;

        OOFEM_LOG_DEBUG( "Solve RVE problem with boundary conditions (type=%u) for macroscale pressure gradient gradP=[%f, %f]\n ", BCType, eps.at(1), eps.at(2) );

        // Compute seepage velocity
        rveE->rveSetBoundaryConditions(BCType, eps);
        rveE->rveGiveCharacteristicData(1, &X, &answer, tStep);

        OOFEM_LOG_DEBUG( "Pressure gradient gradP=[%f %f] yields velocity vector [%f %f]\n", X.at(1), X.at(2), answer.at(1), answer.at(2) );

        status->letTempGradPVectorBe(X);
        status->letTempVelocityVectorBe(answer);

        // Compute tangent
        FloatMatrix K;
        rveE->rveGiveCharacteristicData(2, &X, &K, tStep);
        status->letTempTangentMatrixBe(K);
    }

    OOFEM_LOG_DEBUG("****** Exit giveFluxVector **************************************** \n");

    this->enableStdout();
}

void
rvestokesflow :: exportFilter(EngngModel *E, GaussPoint *gp, TimeStep *tStep)
{
    // Fix name for RVE output file and print output

    std :: string filename, basefilename;
    std :: ostringstream tempstring;

    basefilename = this->giveDomain()->giveEngngModel()->giveOutputBaseFileName();
    tempstring << basefilename << ".rve/Rve_" << gp->giveElement()->giveGlobalNumber() << "_" << gp->giveNumber();
    filename = tempstring.str();

    basefilename = E->giveOutputBaseFileName();
    E->letOutputBaseFileNameBe(filename);

    this->suppressStdout();
    E->doStepOutput(tStep);
    this->enableStdout();

    E->letOutputBaseFileNameBe(basefilename);
}

void
rvestokesflow :: giveCharacteristicMatrix(FloatMatrix &answer,  MatResponseForm form, MatResponseMode, GaussPoint *gp, TimeStep *atTime)
{
    this->suppressStdout();

    OOFEM_LOG_INFO("\n****** Enter giveDeviatoricStiffnessMatrix ********************** rve @ %p \n", this->rve);

    rvestokesflowMaterialStatus *status = ( ( rvestokesflowMaterialStatus * ) this->giveStatus(gp) );

    answer = status->giveTempTangentMatrix();
    status->letTempTangentMatrixBe(answer);

    OOFEM_LOG_INFO("****** Exit giveDeviatoricStiffnessMatrix **************************************** \n");

    this->enableStdout();
}

MaterialStatus *
rvestokesflow :: CreateStatus(GaussPoint *gp) const
{
    return new rvestokesflowMaterialStatus(1, this->giveDomain(), gp);
}

}
