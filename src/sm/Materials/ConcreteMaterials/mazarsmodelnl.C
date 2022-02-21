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

#include "mazarsmodelnl.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "datastream.h"

namespace oofem {
REGISTER_Material(MazarsNLMaterial);

MazarsNLMaterial :: MazarsNLMaterial(int n, Domain *d) : MazarsMaterial(n, d), StructuralNonlocalMaterialExtensionInterface(d)
{}


Interface *
MazarsNLMaterial :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else {
        return nullptr;
    }
}



void
MazarsNLMaterial :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const
{
    /*  Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */
    auto nlstatus = static_cast< MazarsNLMaterialStatus * >( this->giveStatus(gp) );
    FloatArray SDstrainVector;
    double equivStrain;

    this->initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, tStep, VM_Total);

    // compute equivalent strain
    equivStrain = this->computeLocalEquivalentStrain(SDstrainVector, gp, tStep);

    nlstatus->setLocalEquivalentStrainForAverage(equivStrain);
}



double
MazarsNLMaterial :: computeEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< MazarsNLMaterialStatus * >( this->giveStatus(gp) );
    double nonlocalEquivalentStrain = 0.0;

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);

    // compute nonlocal strain increment first
    for ( auto &lir: *this->giveIPIntegrationList(gp) ) {
        auto nonlocStatus = static_cast< MazarsNLMaterialStatus * >( this->giveStatus(lir.nearGp) );
        auto nonlocalContribution = nonlocStatus->giveLocalEquivalentStrainForAverage();
        nonlocalContribution *= lir.weight;

        nonlocalEquivalentStrain += nonlocalContribution;
    }

    nonlocalEquivalentStrain *= 1. / status->giveIntegrationScale();
    this->endIPNonlocalAverage(gp);  // !
    return nonlocalEquivalentStrain;
}

void
MazarsNLMaterial :: initializeFrom(InputRecord &ir)
{
    MazarsMaterial :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, _IFT_MazarsNLMaterial_r);
    if ( R < 0.0 ) {
        R = 0.0;
    }

    this->hReft = this->hRefc = 1.0;
}


double
MazarsNLMaterial :: computeWeightFunction(const double R, const FloatArray &src, const FloatArray &coord) const
{
    // Bell shaped function decaying with the distance.

    double dist = distance(src, coord);

    if ( ( dist >= 0. ) && ( dist <= this->R ) ) {
        double help = ( 1. - dist * dist / ( R * R ) );
        return help * help;
    }

    return 0.0;
}


void
MazarsNLMaterial :: initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) const
{
    /*
     * Perfoms initialization, when damage first appear. The Le characteristic length is
     * set equal to 1.0, it doesnot matter - nonlocal approach is used. The only
     * same value with reference length (which is used in local model, which
     * computeDmaga function is reused).
     */
    MazarsNLMaterialStatus *status = static_cast< MazarsNLMaterialStatus * >( this->giveStatus(gp) );

    status->setLe(1.0);
    status->setLec(1.0);
}


MazarsNLMaterialStatus :: MazarsNLMaterialStatus(GaussPoint *g) :
    MazarsMaterialStatus(g), StructuralNonlocalMaterialStatusExtensionInterface()
{}


void
MazarsNLMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
        fprintf(file, "nonloc-kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
}


void
MazarsNLMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    MazarsMaterialStatus :: initTempStatus();
}


void
MazarsNLMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    MazarsMaterialStatus :: updateYourself(tStep);
}


void
MazarsNLMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    MazarsMaterialStatus :: saveContext(stream, mode);
    //if (!stream.write(&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);
}

void
MazarsNLMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    MazarsMaterialStatus :: restoreContext(stream, mode);
    //if (!stream.read (&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);
}

Interface *
MazarsNLMaterialStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialStatusExtensionInterface * >(this);
    } else {
        return NULL;
    }
}


int
MazarsNLMaterial :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    MazarsNLMaterialStatus *status = static_cast< MazarsNLMaterialStatus * >( this->giveStatus(ip) );

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(tStep);

    return buff.write( status->giveLocalEquivalentStrainForAverage() );
}

int
MazarsNLMaterial :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    int result;
    MazarsNLMaterialStatus *status = static_cast< MazarsNLMaterialStatus * >( this->giveStatus(ip) );
    double localEquivalentStrainForAverage;

    result = buff.read(localEquivalentStrainForAverage);
    status->setLocalEquivalentStrainForAverage(localEquivalentStrainForAverage);
    return result;
}

int
MazarsNLMaterial :: estimatePackSize(DataStream &buff, GaussPoint *ip)
{
    //
    // Note: status localStrainVectorForAverage memeber must be properly sized!
    //
    //MazarsNLMaterialStatus *status = (MazarsNLMaterialStatus*) this -> giveStatus (ip);

    return buff.givePackSizeOfDouble(1);
}

} // end namespace oofem
