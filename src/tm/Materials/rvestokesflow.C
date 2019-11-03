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

#include "tm/Materials/rvestokesflow.h"
#include "util.h"
#include "oofemtxtdatareader.h"
#include "fm/stokesflowvelocityhomogenization.h"
#include "engngm.h"
#include "contextioerr.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

#include <cstdio>
#include <cstring>
#include <sstream>

#ifdef __FM_MODULE

namespace oofem {
REGISTER_Material(RVEStokesFlow);

int RVEStokesFlow :: n = 1;

RVEStokesFlowMaterialStatus :: RVEStokesFlowMaterialStatus(int n, int rank, GaussPoint *g, const std :: string &inputfile) :
    TransportMaterialStatus(g), oldTangent(true)
{
    OOFEM_LOG_INFO( "************************** Instanciating microproblem from file %s\n", inputfile.c_str() );
    OOFEMTXTDataReader dr( inputfile.c_str() );

    auto e = InstanciateProblem(dr, _processor, 0);
    dr.finish();
 
    if ( dynamic_cast< StokesFlowVelocityHomogenization* >(e.get()) ) {
        this->rve.reset( dynamic_cast< StokesFlowVelocityHomogenization* >(e.release()) );
    } else {
        OOFEM_ERROR("Unexpected RVE engineering model");
    }

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
    if ( rank >= 0 ) {
        name << "." << rank;
    }
    this->rve->letOutputBaseFileNameBe( name.str() );

    this->rve->setProblemScale(microScale);
    this->rve->checkProblemConsistency();
    this->rve->initMetaStepAttributes( this->rve->giveMetaStep(1) );
    this->rve->giveNextStep();
    this->rve->init();

    OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", rve.get());
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


void
RVEStokesFlowMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    TransportMaterialStatus :: saveContext(stream, mode);
}


void
RVEStokesFlowMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    TransportMaterialStatus :: restoreContext(stream, mode);
}


RVEStokesFlow :: RVEStokesFlow(int n, Domain *d) : TransportMaterial(n, d)
{ }


void RVEStokesFlow :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, this->rveFilename, _IFT_RVEStokesFlow_fileName);
}


int
RVEStokesFlow :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    auto status = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );
    answer.clear();

    switch ( type ) {
    case IST_Velocity:
        answer = status->giveFlux();
        return 1;
    case IST_PressureGradient:
        answer = status->giveGradient();
        return 1;
    case IST_TangentNorm:
        answer.resize(1);
        answer.at(1) = frobeniusNorm(status->giveTangentMatrix());
        return 1;
    case IST_Tangent:
    {
        const auto &temp = status->giveTangentMatrix();
        answer = {temp(0,0), temp(0,1), temp(0,2), temp(1,0), temp(1,1), temp(1,2), temp(2,0), temp(2,1), temp(2,2)};
        return 1;
    }
    default:
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 0;
}


FloatArrayF<3>
RVEStokesFlow :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_LOG_DEBUG("\n****** Enter giveFluxVector ********************** Element number %u, Gauss point %u\n", 
                    gp->giveElement()->giveGlobalNumber(), gp->giveNumber());

    auto status = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );
    auto rveE = status->giveRVE();

    OOFEM_LOG_DEBUG( "Solve RVE problem for macroscale pressure gradient gradP=[%f, %f, %f]\n ", grad[0], grad[1], grad[2] );

    // Compute seepage velocity
    rveE->applyPressureGradient(grad);
    status->setTimeStep(tStep);
    rveE->solveYourselfAt(rveE->giveCurrentStep());
    FloatArray answer;
    rveE->computeSeepage(answer, tStep);
    answer.resizeWithValues(3);

    OOFEM_LOG_DEBUG( "Pressure gradient gradP=[%f %f] yields velocity vector [%f %f]\n", grad.at(1), grad.at(2), answer.at(1), answer.at(2) );

    status->setTempGradient(grad);
    status->setTempFlux(answer);
    status->oldTangent = true;

    OOFEM_LOG_DEBUG("****** Exit giveFluxVector **************************************** \n");
    return answer;
}


FloatMatrixF<3,3>
RVEStokesFlow :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_LOG_DEBUG("\n****** Enter giveDeviatoricStiffnessMatrix **********************\n");

    auto status = static_cast< RVEStokesFlowMaterialStatus * >( this->giveStatus(gp) );

    if ( status->oldTangent ) {
        // Compute tangent
        FloatMatrix answer;
        status->giveRVE()->computeTangent(answer, tStep);
        answer.resizeWithData(3, 3);
        status->letTempTangentMatrixBe(answer);
        status->oldTangent = false;
        return answer;
    } else {
        return status->giveTempTangentMatrix();
    }

    OOFEM_LOG_DEBUG("****** Exit giveDeviatoricStiffnessMatrix **************************************** \n");
}


MaterialStatus *
RVEStokesFlow :: CreateStatus(GaussPoint *gp) const
{
    int rank = -1;
    if ( this->domain->giveEngngModel()->isParallel() && this->domain->giveEngngModel()->giveNumberOfProcesses() > 1 ) {
        rank = this->domain->giveEngngModel()->giveRank();
    }
    return new RVEStokesFlowMaterialStatus(n++, rank, gp, this->rveFilename);
}

}
#endif // ifdef __FM_MODULE
