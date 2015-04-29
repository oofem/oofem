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

#include "hydratingisoheatmat.h"
#include "gausspoint.h"
#include "timestep.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HydratingIsoHeatMaterial);

IRResultType
HydratingIsoHeatMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    int value;
    double dvalue;

    // set k, c - necessary; rc beton Hellmich 2428 kJ/m3
    result = IsotropicHeatTransferMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    // setup hydration model
    result = HydrationModelInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    dvalue = -2.;
    IR_GIVE_OPTIONAL_FIELD(ir, dvalue, _IFT_HydratingIsoHeatMaterial_hydration);
    if ( dvalue >= 0. ) {
        hydration = 1;
    } else {
        hydration = 0;
    }

    if ( hydration ) {
        // mixture type: 1 - mtLafarge, 2 - mtHuber, 3 - mtC60
        value = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_HydratingIsoHeatMaterial_mix);
        if ( !value ) {
            value = mtLafarge;
        }

        setMixture( ( MixtureType ) value );
        printf("\nHydratingHeatMat %d: using mixture %d.\n", giveNumber(), value);

        if ( ir->hasField(_IFT_HydratingIsoHeatMaterial_noHeat) ) {
            hydrationHeat = 0;
            printf( "HydratingHeatMat %d: hydration heat neglected.\n", giveNumber() );
        } else {
            hydrationHeat = 1;
        }

        if ( hydrationHeat ) {
            // include hydration internal source in LHS?
            if ( ir->hasField(_IFT_HydratingIsoHeatMaterial_noLHS) ) {
                hydrationLHS = 0;
                printf( "HydratingHeatMat %d: hydration heat not included in LHS.\n", giveNumber() );
            } else {
                hydrationLHS = 1;
            }
        }
    }

    return IRRT_OK;
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

int
HydratingIsoHeatMaterial :: hasInternalSource()
// return true if hydration heat source is present
{
    if ( hydrationHeat ) {
        return 1;
    } else {
        return 0;
    }
}

void
HydratingIsoHeatMaterial :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
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
    FloatArray aux;
    if ( ms ) {
        ms->letTempStateVectorBe(vec);
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
            if ( ( gp->giveNumber() == 1 ) && giveStatus(gp) ) {
                FILE *vyst = fopen("teplota.out", "a");
                computeInternalSourceVector(aux, gp, tStep, VM_Incremental);
                if ( aux.isEmpty() ) {
                    aux.resize(1);
                    aux.zero();
                }

                aux.times( 1. / give('d', gp) );
                fprintf( vyst, "Elem %.3d krok %.2d: t= %.0f, dt=%.0f, %ld. it, ksi= %.12f, T= %.8f, heat=%.8f\n", gp->giveElement()->giveNumber(), tStep->giveNumber(),
                        tStep->giveTargetTime(), tStep->giveTimeIncrement(), tStep->giveSolutionStateCounter(),
                        giveHydrationDegree(gp, tStep, VM_Total), vec.at(1), aux.at(1) * tStep->giveTimeIncrement() );
                fclose(vyst);
            }
        }
    }
}

double
HydratingIsoHeatMaterial :: giveCharacteristicValue(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    double answer = 0;
    FloatArray vec;

    if ( rmode == Capacity ) {
        if ( castAt && ( tStep->giveTargetTime() < castAt ) ) {
            answer = capacity * this->give('d', gp) / 1000;                            // Zero capacity before cast
        } else {
            answer = capacity * this->give('d', gp);
        }
    } else if ( !hydrationLHS ) {
        answer = 0;
    } else if ( hydrationModel ) { //!!! better via HydrationModelInterface
        vec = static_cast< TransportMaterialStatus * >( giveStatus(gp) )->giveTempField();
        if ( vec.giveSize() < 2 ) {
            vec.resize(2);
            vec.at(2) = 1.; // saturated if undefined
        }

        answer = hydrationModel->giveCharacteristicValue(vec, rmode, gp, tStep)
        / tStep->giveTimeIncrement();
    } else {
        OOFEM_ERROR("unknown MatResponseMode (%s)", __MatResponseModeToString(rmode) );
    }

    return answer;
}

contextIOResultType
HydratingIsoHeatMaterial :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
// saves full status for this material, also invokes saving
// for sub-objects of this.
{
    contextIOResultType iores;

    // write parent data
    if ( ( iores = TransportMaterial :: saveIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // save hydration model data - maybe should check hydration option?
    if ( ( iores = HydrationModelInterface :: saveContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
HydratingIsoHeatMaterial :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
// restores full status for this material, also invokes restoring for sub-objects of this.
{
    contextIOResultType iores;

    // read parent data
    if ( ( iores = TransportMaterial :: restoreIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read hydration model data - maybe should check hydration option?
    if ( ( iores = HydrationModelInterface :: restoreContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
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
    return new HydratingTransportMaterialStatus(1, domain, gp);
}


void
HydratingTransportMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
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
        return NULL;
    }
}
} // end namespace oofem
