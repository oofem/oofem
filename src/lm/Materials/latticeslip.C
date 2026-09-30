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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#include "latticeslip.h"
#include "latticelinearelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "floatarrayf.h"
#include "Elements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(LatticeSlip);

LatticeSlip::LatticeSlip(int n, Domain *d) : LatticeLinearElastic(n, d)
{}


bool
LatticeSlip::hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeSlip::initializeFrom(const std::shared_ptr < InputRecord > & ir)
{
    LatticeLinearElastic::initializeFrom(ir);


    //Parameter which relates the shear stiffness to the normal stiffness. Reread this again because we use the stiffness matrix from latticelinearelastic material where the default is 1.
    IR_GIVE_FIELD(ir, this->alphaOne, _IFT_LatticeSlip_a1); // Macro

    //Parameter which is used for the definition of bending stiffness. Reread this again because we use the stiffness matrix from latticelinearelastic material where the default is 0. It cannot be zero for latticeslip. Optional; defaults to a1 for backward compatibility.
    this->alphaTwo = this->alphaOne;
    IR_GIVE_OPTIONAL_FIELD(ir, this->alphaTwo, _IFT_LatticeSlip_a2); // Macro

    //Parameter which limits the stress in slip direction.
    IR_GIVE_FIELD(ir, tauMax, _IFT_LatticeSlip_t0); // Macro

    //Three models available
    //0 default. Linear elastic-perfect plastic
    //1 bond slip according first part of CEB model and then constant
    //2 bond slip according to CEB model
    IR_GIVE_OPTIONAL_FIELD(ir, this->type, _IFT_LatticeSlip_type); // Macro

    if ( type == 1 || type == 2 ) {
        IR_GIVE_FIELD(ir, s1, _IFT_LatticeSlip_s1);
        IR_GIVE_FIELD(ir, alpha, _IFT_LatticeSlip_alpha);
    }

    if ( type == 2 ) {
        IR_GIVE_FIELD(ir, s2, _IFT_LatticeSlip_s2);
        IR_GIVE_FIELD(ir, s3, _IFT_LatticeSlip_s3);
        IR_GIVE_FIELD(ir, tauFinal, _IFT_LatticeSlip_tf);
    }
}


std::unique_ptr < MaterialStatus >
LatticeSlip::CreateStatus(GaussPoint * gp) const
{
    return std::make_unique < LatticeSlipStatus > ( gp );
}


FloatArrayF < 6 >
LatticeSlip::giveLatticeStress3d(const FloatArrayF < 6 > & totalStrain, GaussPoint * gp, TimeStep * atTime)
{
    auto status = static_cast < LatticeSlipStatus * > ( this->giveStatus(gp) );
    status->initTempStatus();

    auto oldStress = status->giveLatticeStress();
    auto oldStrain = status->giveLatticeStrain();

    auto stiffnessMatrix = this->give3dLatticeStiffnessMatrix(ElasticStiffness, gp, atTime);

    //evaluate tempKappa (no elastic strain in axial direction)
    double tempKappa = status->giveKappa() + fabs(totalStrain.at(1) - oldStrain.at(1) );

    double deltaSlip = totalStrain.at(1) - oldStrain.at(1);

    /*First component is the slip one for which the stress should be limited using plasiticity (frictional slip between fibre and matrix). The other components are kept elastic. */
    FloatArrayF < 6 > stress;
    stress.at(1) = oldStress.at(1) +
        ( totalStrain.at(1) - oldStrain.at(1) ) * stiffnessMatrix.at(1, 1);

    double bondStress = evaluateBondStress(tempKappa);
    double f = fabs(stress.at(1) ) - bondStress;

    if ( f > 0 ) {//plastic response.
        //Reduced stress by increasing plastic strain.
        stress.at(1) = sgn(stress.at(1) ) * bondStress;
        status->setTempCrackFlag(1);
    }



    //Compute the final stress components
    for ( int i = 2; i <= 6; i++ ) {
        stress.at(i) =  stiffnessMatrix.at(i, i) * totalStrain.at(i);
    }

    status->letTempKappaBe(tempKappa);
    status->letTempLatticeStrainBe(totalStrain);
    status->letTempLatticeStressBe(stress);




    double tempDissipation = status->giveDissipation();
    double tempDeltaDissipation;

    tempDeltaDissipation = computeDeltaDissipation(gp, atTime);

    tempDissipation += tempDeltaDissipation;

    //Set all temp values
    status->setTempDissipation(tempDissipation);
    status->setTempDeltaDissipation(tempDeltaDissipation);

    return stress;
}


double
LatticeSlip::evaluateBondStress(const double kappa) const
{
    if ( this->type == 0 ) { //elastic-perfectly plastic
        return this->tauMax;
    } else if ( this->type == 1 ) { //modified bond model
        if ( kappa <= 0. ) {
            return 0.;
        }
        if ( kappa <= s1 ) {
            return this->tauMax * pow(kappa / s1, alpha);
        }
        return tauMax;
    } else if ( this->type == 2 ) {
        if ( kappa <= 0. ) {
            return 0.;
        }
        if ( kappa <= s1 ) {
            return tauMax * pow(kappa / s1, alpha);
        }
        if ( kappa <= s2 ) {
            return tauMax;
        }
        if ( kappa <= s3 ) {
            return tauMax - ( tauMax - tauFinal ) * ( kappa - s2 ) / ( s3 - s2 );
        }
        return tauFinal;
    } else { //unknown type
        OOFEM_ERROR("Unknown bond model type. Type should be 0, 1 or 2.");
    }

    return 0.; //Should not be here.
}



Interface *
LatticeSlip::giveInterface(InterfaceType type)
{
    return nullptr;
}


LatticeSlipStatus::LatticeSlipStatus(GaussPoint *g) :  LatticeMaterialStatus(g)
{}


void
LatticeSlipStatus::initTempStatus()
{
    LatticeMaterialStatus::initTempStatus();
    this->tempKappa = this->kappa;
}


FloatArrayF < 6 >
LatticeSlip::giveThermalDilatationVector(GaussPoint * gp, TimeStep * tStep) const
//
// returns a FloatArray(6) of initial strain vector
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);

    //Option to add a eigendisplacement instead of strain
    double length = static_cast < LatticeStructuralElement * > ( gp->giveElement() )->giveLength();
    alpha += this->cAlpha / length;

    return {
        alpha, 0., 0., 0., 0., 0.
    };
}

void
LatticeSlipStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus::printOutputAt(file, tStep);
    fprintf(file, "kappa %.8e, dissipation %f, deltaDissipation %f, crackFlag %d \n", this->kappa, this->dissipation, this->deltaDissipation, this->crackFlag);
}


double
LatticeSlip::computeDeltaDissipation(GaussPoint *gp, TimeStep *atTime) const
{
    auto status = static_cast < LatticeSlipStatus * > ( this->giveStatus(gp) );

    const double kappa = status->giveKappa();
    const double tempKappa = status->giveTempKappa();
    const auto &tempStress = status->giveTempLatticeStress();
    const auto &stress = status->giveLatticeStress();

    return ( tempStress.at(1) + stress.at(1) )  * ( tempKappa - kappa ) / 2.;
}


void
LatticeSlipStatus::saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    LatticeMaterialStatus::saveContext(stream, mode);

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
LatticeSlipStatus::restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    LatticeMaterialStatus::restoreContext(stream, mode);
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
LatticeSlipStatus::updateYourself(TimeStep *atTime)
{
    LatticeMaterialStatus::updateYourself(atTime);
    this->kappa = this->tempKappa;
}



int
LatticeSlip::giveIPValue(FloatArray &answer,
                         GaussPoint *gp,
                         InternalStateType type,
                         TimeStep *atTime)
{
    auto status = static_cast < LatticeSlipStatus * > ( this->giveStatus(gp) );

    if ( type == IST_CrackStatuses ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackFlag();
        return 1;
    } else if ( type == IST_DissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDissipation();
        return 1;
    } else if ( type == IST_DeltaDissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDeltaDissipation();
        return 1;
    } else if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = static_cast < LatticeStructuralElement * > ( gp->giveElement() )->giveLength();
        return 1;
    } else {
        return LatticeLinearElastic::giveIPValue(answer, gp, type, atTime);
    }
}
}
