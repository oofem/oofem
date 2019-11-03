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

#include "tm/Materials/hydratingisoheatmat.h"
#include "gausspoint.h"
#include "timestep.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HydratingIsoHeatMaterial);

void
HydratingIsoHeatMaterial :: initializeFrom(InputRecord &ir)
{
    int value;
    double dvalue;

    // set k, c - necessary; rc beton Hellmich 2428 kJ/m3
    IsotropicHeatTransferMaterial :: initializeFrom(ir);

    // setup hydration model
    HydrationModelInterface :: initializeFrom(ir);

    dvalue = -2.;
    IR_GIVE_OPTIONAL_FIELD(ir, dvalue, _IFT_HydratingIsoHeatMaterial_hydration);
    hydration = dvalue >= 0.;

    if ( hydration ) {
        // mixture type: 1 - mtLafarge, 2 - mtHuber, 3 - mtC60
        value = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_HydratingIsoHeatMaterial_mix);
        if ( !value ) {
            value = mtLafarge;
        }

        setMixture( ( MixtureType ) value );
        printf("\nHydratingHeatMat %d: using mixture %d.\n", giveNumber(), value);

        if ( ir.hasField(_IFT_HydratingIsoHeatMaterial_noHeat) ) {
            hydrationHeat = false;
            printf( "HydratingHeatMat %d: hydration heat neglected.\n", giveNumber() );
        } else {
            hydrationHeat = true;
        }

        if ( hydrationHeat ) {
            // include hydration internal source in LHS?
            if ( ir.hasField(_IFT_HydratingIsoHeatMaterial_noLHS) ) {
                hydrationLHS = false;
                printf( "HydratingHeatMat %d: hydration heat not included in LHS.\n", giveNumber() );
            } else {
                hydrationLHS = true;
            }
        }
    }
}

void
HydratingIsoHeatMaterial :: setMixture(MixtureType mix)
// creates the hydration model instance if necessary, sets the mixture type
{
    if ( hydrationModel ) {
        hydrationModel->setMixture(mix);
    } else if ( hydration ) {
        OOFEM_ERROR("Can't setup undefined hydrationModel.");
    }
}

bool
HydratingIsoHeatMaterial :: hasInternalSource() const
{
    return hydrationHeat;
}

void
HydratingIsoHeatMaterial :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
// returns in val the hydration heat computed by the hydration model for given hydration degree increment
// current hydration model returns heat in (k)J/m3.
// maybe??? element expects J/kg -> would have to divide by density here
// rate of internal source must be returned, it is multiplied by time increment in element integration.
{
    if ( hydrationHeat ) {
        if ( hydrationModel ) { ///@todo better via HydrationModelInterface
            hydrationModel->computeInternalSourceVector(val, gp, tStep, VM_Incremental); ///@todo mode is VM_Total for nltransientstatic
            val.times( 1. / tStep->giveTimeIncrement() ); // /give('d');
        } else {
            val.zero();
        }

        /*
         * printf("HIsoHeatMat: Ksi %.4f, dksi %.4f, heat %g\n",
         * giveHydrationDegree(gp, tStep, VM_Total), giveHydrationDegree(gp, tStep, VM_Incremental), (val.giveSize())?val.at(1):0);
         */
    } else {
        val.clear();
    }
}

void
HydratingIsoHeatMaterial :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    if ( ms ) {
        ms->setTempField(vec[0]);

        if ( hydration ) {
            /* OBSOLETE
             * FloatArray s = ms->giveStateVector ();
             * if (vec.isEmpty()) OOFEM_ERROR("empty new state vector");
             * aux.resize(2);
             * aux.at(1) = vec.at(1);
             * if (s.isEmpty()||(tStep->giveTime()<=0)) aux.at(2) = initialHydrationDegree; // apply initial conditions
             * else {
             *  aux.at(2) = s.at(2);
             *  if (!castAt || (tStep->giveTime()>=castAt)) aux.at(2) += hydrationModel->dksi (s.at(2), vec.at(1), tStep->giveTimeIncrement()); // compute hydration degree increment
             * }
             */
            HydrationModelInterface :: updateInternalState(vec, gp, tStep);

            // additional file output !!!
            if ( gp->giveNumber() == 1 && giveStatus(gp) ) {
                FILE *vyst = fopen("teplota.out", "a");
                FloatArray aux;
                computeInternalSourceVector(aux, gp, tStep, VM_Incremental);
                if ( aux.isEmpty() ) {
                    aux.resize(1);
                    aux.zero();
                }

                aux.times( 1. / give('d', gp, tStep) );
                fprintf( vyst, "Elem %.3d krok %.2d: t= %.0f, dt=%.0f, %ld. it, ksi= %.12f, T= %.8f, heat=%.8f\n", gp->giveElement()->giveNumber(), tStep->giveNumber(),
                        tStep->giveTargetTime(), tStep->giveTimeIncrement(), tStep->giveSolutionStateCounter(),
                        giveHydrationDegree(gp, tStep, VM_Total), vec.at(1), aux.at(1) * tStep->giveTimeIncrement() );
                fclose(vyst);
            }
        }
    }
}

double
HydratingIsoHeatMaterial :: giveCharacteristicValue(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    if ( rmode == Capacity ) {
        if ( castAt && ( tStep->giveTargetTime() < castAt ) ) {
            return this->give('c', gp, tStep) * this->give('d', gp, tStep) / 1000;                            // Zero capacity before cast
        } else {
            return this->give('c', gp, tStep) * this->give('d', gp, tStep);
        }
    } else if ( !hydrationLHS ) {
        return 0;
    } else if ( hydrationModel ) { //!!! better via HydrationModelInterface
        auto status = static_cast< HeMoTransportMaterialStatus * >( giveStatus(gp) );
        double t = status->giveTempTemperature();
        double h = status->giveTempHumidity(); // TODO CHECK

        return hydrationModel->giveCharacteristicValue(t, h, rmode, gp, tStep) / tStep->giveTimeIncrement();
    } else {
        OOFEM_ERROR("unknown MatResponseMode (%s)", __MatResponseModeToString(rmode) );
        return 0.;
    }
}

void
HydratingIsoHeatMaterial :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
{
    TransportMaterial :: saveIPContext(stream, mode, gp);

    // save hydration model data - maybe should check hydration option?
    HydrationModelInterface :: saveContext(stream, mode);
}

void
HydratingIsoHeatMaterial :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
{
    TransportMaterial :: restoreIPContext(stream, mode, gp);

    // read hydration model data - maybe should check hydration option?
    HydrationModelInterface :: restoreContext(stream, mode);
}

int
HydratingIsoHeatMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    // printf ("IP %d::giveIPValue, IST %d", giveNumber(), type);
    if ( type == IST_HydrationDegree ) {
        //TransportMaterialStatus* status = (TransportMaterialStatus*) this -> giveStatus (gp);
        answer.resize(1);
        //if (hydration)
        answer.at(1) = giveHydrationDegree(gp, tStep, VM_Total);
        //else answer.at(1) = 0;
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

MaterialStatus *
HydratingIsoHeatMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new HydratingTransportMaterialStatus(gp);
}


void
HydratingTransportMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    fprintf(file, " status ");
    HydrationModelStatusInterface :: printOutputAt(file, tStep);
    TransportMaterialStatus :: printOutputAt(file, tStep);
}

// necessary for proper cast to interface, can't be done from outside
Interface *
HydratingTransportMaterialStatus :: giveInterface(InterfaceType type)
{
    if ( type == HydrationModelStatusInterfaceType ) {
        return static_cast< HydrationModelStatusInterface * >(this);
    } else {
        return nullptr;
    }
}
} // end namespace oofem
