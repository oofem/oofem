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

#include "rvestokesflow.h"
#include "util.h"
#include "rveengngmodel.h"
#include "engngm.h"
#include "nummet.h"
#include "contextioerr.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

#include <cstdio>
#include <cstring>
#include <sstream>

namespace oofem {
REGISTER_Material(RVEStokesFlow);

RVEStokesFlowMaterialStatus :: RVEStokesFlowMaterialStatus(int n, Domain *d, GaussPoint *g, EngngModel *rve) :
    TransportMaterialStatus(n, d, g)
{
    temp_gradient.resize(2);
    temp_gradient.zero();

    temp_flux.resize(2);
    temp_flux.zero();

    temp_TangentMatrix.resize(2, 2);
    temp_TangentMatrix.zero();

    solutionVector = new FloatArray;

    this->rve = rve;
}

RVEStokesFlowMaterialStatus :: ~RVEStokesFlowMaterialStatus()
{ }

void
RVEStokesFlowMaterialStatus :: exportFilter(GaussPoint *gp, TimeStep *tStep)
{
    // Fix name for RVE output file and print output

    std :: string filename, basefilename;
    std :: ostringstream filenameStringStream;


    filenameStringStream << this->giveDomain()->giveEngngModel()->giveOutputBaseFileName().c_str() << ".rve/Rve_" << gp->giveElement()->giveGlobalNumber() << "_" << gp->giveNumber();

    filename = filenameStringStream.str();

    // Update fields on subscale
    FloatArray grapP = this->giveTempGradient(), seepageVelocity;

    rveEngngModel *rveE;
    rveE = dynamic_cast< rveEngngModel * >(this->rve);

    rveE->rveSetBoundaryConditions(10, grapP);
    rveE->rveGiveCharacteristicData(1, & grapP, & seepageVelocity, tStep);

    basefilename = this->rve->giveOutputBaseFileName();

    this->rve->letOutputBaseFileNameBe(filename);
    this->rve->doStepOutput(tStep);
    this->rve->letOutputBaseFileNameBe(basefilename);
}

void
RVEStokesFlowMaterialStatus :: initTempStatus()
{
    TransportMaterialStatus :: initTempStatus();
}

void
RVEStokesFlowMaterialStatus :: updateYourself(TimeStep *tStep)
{
    TransportMaterialStatus :: updateYourself(tStep);
    tangentMatrix = temp_TangentMatrix;

    //Mtrl->suppressStdout();
    //Mtrl->enableStdout();
}


contextIOResultType
RVEStokesFlowMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        OOFEM_ERROR("saveContex : can't write into NULL stream");
    }

    if ( ( iores = TransportMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
RVEStokesFlowMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        OOFEM_ERROR("saveContex : can't write into NULL stream");
    }

    if ( ( iores = TransportMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

RVEStokesFlow :: RVEStokesFlow(int n, Domain *d) : RVEMaterial(n, d), TransportMaterial(n, d)
{ }

IRResultType RVEStokesFlow :: initializeFrom(InputRecord *ir)
{
    // this->TransportMaterial :: initializeFrom(ir);
    this->RVEMaterial :: initializeFrom(ir);

    return IRRT_OK;
}

int
RVEStokesFlow :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    RVEStokesFlowMaterialStatus *thisMaterialStatus;
    thisMaterialStatus = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );
    FloatMatrix temp;
    answer.resize(3);
    answer.zero();

    switch ( type ) {
    case IST_Velocity:
        answer.copySubVector(thisMaterialStatus->giveFlux(), 1);
        break;
    case IST_PressureGradient:
        answer.copySubVector(thisMaterialStatus->giveGradient(), 1);
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
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 1;
}

void
RVEStokesFlow :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    this->suppressStdout();

    OOFEM_LOG_DEBUG("\n****** Enter giveFluxVector ********************** Element number %u, Gauss point %u, rve @ %p \n", gp->giveElement()->giveGlobalNumber(), gp->giveNumber(), this->rve);

    RVEStokesFlowMaterialStatus *status = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );
    FloatArray temp = status->giveTempGradient();

    if ( temp.giveSize() >= 3 && temp.at(1) == grad.at(1) && temp.at(2) == grad.at(2) ) {
        OOFEM_LOG_DEBUG("This pressure gradient has already been evaluated\n");
        answer = status->giveTempFlux();
    } else {
        FloatArray X;
        rveEngngModel *rveE;

        rveE = dynamic_cast< rveEngngModel * >(this->rve);

        X = grad;

        OOFEM_LOG_DEBUG( "Solve RVE problem with boundary conditions (type=%u) for macroscale pressure gradient gradP=[%f, %f]\n ", BCType, grad.at(1), grad.at(2) );

        // Compute seepage velocity
        rveE->rveSetBoundaryConditions(BCType, grad);
        rveE->rveGiveCharacteristicData(1, & X, & answer, tStep);

        OOFEM_LOG_DEBUG( "Pressure gradient gradP=[%f %f] yields velocity vector [%f %f]\n", X.at(1), X.at(2), answer.at(1), answer.at(2) );


        status->setTempGradient(X);
        status->setTempFlux(answer);

        // Compute tangent
        FloatMatrix K;
        rveE->rveGiveCharacteristicData(2, & X, & K, tStep);
        status->letTempTangentMatrixBe(K);
    }

    OOFEM_LOG_DEBUG("****** Exit giveFluxVector **************************************** \n");

    this->enableStdout();
}

void
RVEStokesFlow :: exportFilter(EngngModel *E, GaussPoint *gp, TimeStep *tStep)
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
RVEStokesFlow :: giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep)
{
    this->suppressStdout();

    OOFEM_LOG_DEBUG("\n****** Enter giveDeviatoricStiffnessMatrix ********************** rve @ %p \n", this->rve);

    RVEStokesFlowMaterialStatus *status = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );

    answer = status->giveTempTangentMatrix();
    status->letTempTangentMatrixBe(answer);

    OOFEM_LOG_DEBUG("****** Exit giveDeviatoricStiffnessMatrix **************************************** \n");

    this->enableStdout();
}

MaterialStatus *
RVEStokesFlow :: CreateStatus(GaussPoint *gp) const
{
    return new RVEStokesFlowMaterialStatus(1, this->giveDomain(), gp, this->rve);
}
}
