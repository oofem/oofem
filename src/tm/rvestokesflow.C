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
#include "oofemtxtdatareader.h"
#include "../fm/stokesflowvelocityhomogenization.h"
#include "engngm.h"
#include "contextioerr.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

#include <cstdio>
#include <cstring>
#include <sstream>

namespace oofem {
REGISTER_Material(RVEStokesFlow);

int RVEStokesFlow :: n = 1;

RVEStokesFlowMaterialStatus :: RVEStokesFlowMaterialStatus(int n, Domain *d, GaussPoint *g, const std :: string &inputfile) :
    TransportMaterialStatus(n, d, g), oldTangent(true)
{
    OOFEM_LOG_INFO( "************************** Instanciating microproblem from file %s\n", inputfile.c_str() );
    OOFEMTXTDataReader dr( inputfile.c_str() );

    EngngModel *e = InstanciateProblem(& dr, _processor, 0);
    dr.finish();
 
    if ( dynamic_cast< StokesFlowVelocityHomogenization* >(e) ) {
        this->rve.reset( dynamic_cast< StokesFlowVelocityHomogenization* >(e) );
    } else {
        delete e;
        OOFEM_ERROR("Unexpected RVE engineering model");
    }

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
    if ( this->domain->giveEngngModel()->isParallel() && this->domain->giveEngngModel()->giveNumberOfProcesses() > 1 ) {
        name << "." << this->domain->giveEngngModel()->giveRank();
    }
    this->rve->letOutputBaseFileNameBe( name.str() );

    this->rve->setProblemScale(microScale);
    this->rve->checkProblemConsistency();
    this->rve->initMetaStepAttributes( this->rve->giveMetaStep(1) );
    this->rve->giveNextStep();
    this->rve->init();

    OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", rve.get());
}

RVEStokesFlowMaterialStatus :: ~RVEStokesFlowMaterialStatus()
{
}


void RVEStokesFlowMaterialStatus :: setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep();
    rveTStep->setNumber( tStep->giveNumber() );
    rveTStep->setTime( tStep->giveTargetTime() );
    rveTStep->setTimeIncrement( tStep->giveTimeIncrement() );
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

    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);
}


contextIOResultType
RVEStokesFlowMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = TransportMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
RVEStokesFlowMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = TransportMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

RVEStokesFlow :: RVEStokesFlow(int n, Domain *d) : TransportMaterial(n, d)
{ }

IRResultType RVEStokesFlow :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, this->rveFilename, _IFT_RVEStokesFlow_fileName);

    SupressRVEoutput = 0;

    IR_GIVE_OPTIONAL_FIELD(ir, SupressRVEoutput, _IFT_RVEStokesFlow_supressoutput);

    return IRRT_OK;
}


int
RVEStokesFlow :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    RVEStokesFlowMaterialStatus *thisMaterialStatus;
    thisMaterialStatus = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );
    FloatMatrix temp;
    answer.clear();

    switch ( type ) {
    case IST_Velocity:
        answer.copySubVector(thisMaterialStatus->giveFlux(), 1);
        return 1;
    case IST_PressureGradient:
        answer.copySubVector(thisMaterialStatus->giveGradient(), 1);
        return 1;
    case IST_TangentNorm:
        temp = thisMaterialStatus->giveTangentMatrix();
        answer.resize(1);
        answer.at(1) = temp.computeFrobeniusNorm();
        return 1;
    case IST_Tangent:
        temp = thisMaterialStatus->giveTangentMatrix();
        answer.resize(4);
        answer.at(1) = temp.at(1, 1);
        answer.at(2) = temp.at(1, 2);
        answer.at(3) = temp.at(2, 1);
        answer.at(4) = temp.at(2, 2);
        return 1;
    default:
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 0;
}

void
RVEStokesFlow :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{

    this->suppressStdout();

    OOFEM_LOG_DEBUG("\n****** Enter giveFluxVector ********************** Element number %u, Gauss point %u\n", 
                    gp->giveElement()->giveGlobalNumber(), gp->giveNumber());

    RVEStokesFlowMaterialStatus *status = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );
    StokesFlowVelocityHomogenization *rveE = status->giveRVE();

    OOFEM_LOG_DEBUG( "Solve RVE problem for macroscale pressure gradient gradP=[%f, %f, %f]\n ",
                     grad.at(1), grad.at(2), grad.giveSize() == 3 ? grad.at(3) : 0. );

    // Compute seepage velocity
    rveE->applyPressureGradient(grad);
    status->setTimeStep(tStep);
    rveE->solveYourselfAt(rveE->giveCurrentStep());
    rveE->computeSeepage(answer, tStep);

    OOFEM_LOG_DEBUG( "Pressure gradient gradP=[%f %f] yields velocity vector [%f %f]\n", grad.at(1), grad.at(2), answer.at(1), answer.at(2) );


    status->setTempGradient(grad);
    status->setTempFlux(answer);
    status->oldTangent = true;

    OOFEM_LOG_DEBUG("****** Exit giveFluxVector **************************************** \n");

    this->enableStdout();
}


void
RVEStokesFlow :: giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep)
{
    this->suppressStdout();

    OOFEM_LOG_DEBUG("\n****** Enter giveDeviatoricStiffnessMatrix **********************\n");

    RVEStokesFlowMaterialStatus *status = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );

    if ( status->oldTangent ) {
        // Compute tangent
        status->giveRVE()->computeTangent(answer, tStep);
        status->letTempTangentMatrixBe(answer);
        status->oldTangent = false;
    } else {
        answer = status->giveTempTangentMatrix();
    }

    OOFEM_LOG_DEBUG("****** Exit giveDeviatoricStiffnessMatrix **************************************** \n");

    this->enableStdout();
}

MaterialStatus *
RVEStokesFlow :: CreateStatus(GaussPoint *gp) const
{
    return new RVEStokesFlowMaterialStatus(n++, this->giveDomain(), gp, this->rveFilename);
}

void
RVEStokesFlow :: suppressStdout()
{
    //    if (SupressRVEoutput) {
    //        fgetpos(stdout, &stdoutPos);
    //        stdoutFID=dup(fileno(stdout));
    //        freopen(this->rveLogFilename.c_str(), "a", stdout);
    //    }
}

void RVEStokesFlow :: enableStdout()
{
    //    if (SupressRVEoutput) {
    //        fflush(stdout);
    //        dup2(stdoutFID, fileno(stdout));
    //        close (stdoutFID);
    //        clearerr(stdout);
    //        fsetpos(stdout, &stdoutPos);        /* for C9X */
    //    }
}

}
