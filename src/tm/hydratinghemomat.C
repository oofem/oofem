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

#include "hydratinghemomat.h"
#include "hydratingisoheatmat.h"
#include "gausspnt.h"
#include "timestep.h"
#include "contextioerr.h"

namespace oofem {
#define PRECAST_CAPACITY_COEFF 1e-2 // coefficient for obtaining capacity before cast of the material : 1e-4, tried 1e-2 for jete (no convergency with 1e-4)

IRResultType
HydratingHeMoMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro
    int value;
    double dvalue;

    // set k, c - necessary; rc beton Hellmich 2428 kJ/m3
    HeMoTKMaterial :: initializeFrom(ir);
    // setup hydration model
    HydrationModelInterface :: initializeFrom(ir);

    dvalue = -2.;
    IR_GIVE_OPTIONAL_FIELD(ir, dvalue, IFT_HydratingHeMoMaterial_hydration, "hydration"); // Macro
    if ( dvalue >= 0. ) {
        hydration = 1;
    } else {
        hydration = 0;
    }

    /* if (ir->hasField("tout")) {
     * teplotaOut = 1;
     * printf("HydratingHeMoMat %d: additional teplota.out output selected.\n", giveNumber());
     * } else */teplotaOut = 0;

    if ( hydration ) {
        // mixture type: 1 - mtLafarge, 2 - mtHuber, 3 - mtC60
        value = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HydratingHeMoMaterial_mix, "mix"); // Macro
        if ( !value ) {
            value = mtLafarge;
        }

        setMixture( ( MixtureType ) value );
        printf("\nHydratingHeMoMat %d: using mixture %d.\n", giveNumber(), value);

        if ( ir->hasField(IFT_HydratingHeMoMaterial_noHeat, "noheat") ) {
            hydrationHeat = 0;
            printf( "HydratingHeMoMat %d: hydration heat neglected.\n", giveNumber() );
        } else {
            hydrationHeat = 1;
        }

        if ( hydrationHeat ) {
            // include hydration internal source in LHS?
            if ( ir->hasField(IFT_HydratingHeMoMaterial_noLHS, "nolhs") ) {
                hydrationLHS = 0;
                printf( "HydratingHeMoMat %d: hydration heat not included in LHS.\n", giveNumber() );
            } else {
                hydrationLHS = 1;
            }
        }
    } else {
        hydrationHeat = false;
        hydrationLHS = false;
    }

    return IRRT_OK;
}

void
HydratingHeMoMaterial :: setMixture(MixtureType mix)
// creates the hydration model instance if necessary, sets the mixture type
{
    if ( hydrationModel ) {
        hydrationModel->setMixture(mix);
    } else if ( hydration ) {
        _error("setMixture: Can't setup undefined hydrationModel.");
    }
}

int
HydratingHeMoMaterial :: hasInternalSource()
// return true if hydration heat source is present
{
    if ( hydrationHeat ) {
        return 1;
    } else {
        return 0;
    }
}

void
HydratingHeMoMaterial :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
// returns in val the hydration heat computed by the hydration model for given hydration degree increment
// current hydration model returns heat in (k)J/m3.
// maybe??? element expects J/kg -> would have to divide by density here
// rate of internal source must be returned, it is multiplied by time increment in element integration.
{
    if ( hydrationHeat ) {
        if ( hydrationModel ) { //!!! better via HydrationModelInterface
            hydrationModel->computeInternalSourceVector(val, gp, atTime, VM_Incremental); //!!! mode is VM_Total for nltransientstatic
            val.times( 1. / atTime->giveTimeIncrement() ); // /give('d');
        } else {
            val.zero();
        }

        /*
         * printf("HIsoHeatMat: Ksi %.4f, dksi %.4f, heat %g\n",
         * giveHydrationDegree(gp, atTime, VM_Total), giveHydrationDegree(gp, atTime, VM_Incremental), (val.giveSize())?val.at(1):0);
         */
    } else {
        val.resize(0);
    }
}

void
HydratingHeMoMaterial :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *atTime)
{
    TransportMaterialStatus *ms = ( TransportMaterialStatus * ) this->giveStatus(gp);
    FloatArray aux;
    if ( ms ) {
        ms->letTempStateVectorBe(vec);
        if ( hydration ) {
            /* OBSOLETE
             * FloatArray s = ms->giveStateVector ();
             * if (vec.isEmpty()) _error("updateInternalState: empty new state vector");
             * aux.resize(2);
             * aux.at(1) = vec.at(1);
             * if (s.isEmpty()||(atTime->giveTime()<=0)) aux.at(2) = initialHydrationDegree; // apply initial conditions
             * else {
             *  aux.at(2) = s.at(2);
             *  if (!castAt || (atTime->giveTime()>=castAt)) aux.at(2) += hydrationModel->dksi (s.at(2), vec.at(1), atTime->giveTimeIncrement()); // compute hydration degree increment
             * }
             */
            // it is necessary to convert the passed state vector to relative humidity expected by the hydration model
            //!!! might be cleaner to choose wc / h in hydration model, but it must be defined which one is passed anyway; so relative humidity was chosen
            //!!! also, the humidity vector might be evaluated by a function (ensure 2 elements and set humidity)
            FloatArray vech = vec;
            if ( vech.giveSize() >= 2 ) {
                vech.at(2) = inverse_sorption_isotherm( vec.at(2) );           // compute relative humidity
            } else {
                vech.resize(2);
                vech.at(2) = 1.; // saturated if undefined
            }

            HydrationModelInterface :: updateInternalState(vech, gp, atTime);

            // additional file output !!!
            if ( teplotaOut && ( gp->giveNumber() == 1 ) && giveStatus(gp) ) {
                FILE *vyst = fopen("teplota.out", "a");
                computeInternalSourceVector(aux, gp, atTime, VM_Incremental);
                if ( aux.isEmpty() ) {
                    aux.resize(1);
                    aux.zero();
                }

                aux.times( 1. / give('d', gp) );
                fprintf( vyst, "Elem %.3d krok %.2d: t= %.0f, dt=%.0f, %ld. it, ksi= %.12f, T= %.8f, heat=%.8f\n", gp->giveElement()->giveNumber(), atTime->giveNumber(),
                        atTime->giveTargetTime(), atTime->giveTimeIncrement(), atTime->giveSolutionStateCounter(),
                        giveHydrationDegree(gp, atTime, VM_Total), vec.at(1), aux.at(1) * atTime->giveTimeIncrement() );
                fclose(vyst);
            }
        }
    }
}

double
HydratingHeMoMaterial :: giveCharacteristicValue(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime)
{
    double answer = 0;
    FloatArray *vec;

    if ( ( rmode >= Capacity_ww ) && ( rmode <= Capacity_wh ) ) { // standard HeMoTK values
        answer = HeMoTKMaterial :: giveCharacteristicValue(rmode, gp, atTime);
        if ( castAt && ( atTime->giveTargetTime() < castAt ) ) {
            answer *= PRECAST_CAPACITY_COEFF;                                  // ~Zero capacity before cast
        }
    } else if ( ( rmode >= IntSource_ww ) && ( rmode <= IntSource_wh ) ) {         // Internal source values
        if ( !hydrationLHS ) {
            answer = 0;
        } else if ( hydrationModel ) {  //!!! better via HydrationModelInterface
            vec = new FloatArray( ( ( TransportMaterialStatus * ) giveStatus(gp) )->giveTempStateVector() );

            if ( vec->giveSize() < 2 ) {
                vec->resize(2);
                vec->at(2) = 1.0; // saturated if undefined
            } else {
                vec->at(2) = inverse_sorption_isotherm( vec->at(2) ); // compute relative humidity
            }

            answer = hydrationModel->giveCharacteristicValue(* vec, rmode, gp, atTime)
                     / atTime->giveTimeIncrement();
            if ( ( rmode == IntSource_ww ) || ( rmode == IntSource_hw ) ) {
                answer *= give_dphi_dw( vec->at(2) );
            }

            delete vec;
        }
    } else {
        _error2( "giveCharacteristicValue: unknown MatResponseMode (%s)", __MatResponseModeToString(rmode) );
    }

    return answer;
}

contextIOResultType
HydratingHeMoMaterial :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
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
HydratingHeMoMaterial :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
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
HydratingHeMoMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    // printf ("IP %d::giveIPValue, IST %d", giveNumber(), type);
    if ( type == IST_HydrationDegree ) {
        //TransportMaterialStatus* status = (TransportMaterialStatus*) this -> giveStatus (aGaussPoint);
        answer.resize(1);
        // zh 24/08/2004 hydration should be selected in HydrationModelInterface->giveHydrationDegree()
        //if (hydration)
        answer.at(1) = giveHydrationDegree(aGaussPoint, atTime, VM_Total);
        //else answer.at(1) = 0;
        return 1;
    } else {
        return HeMoTKMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

InternalStateValueType
HydratingHeMoMaterial :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_HydrationDegree ) {
        return ISVT_SCALAR;
    } else {
        return HeMoTKMaterial :: giveIPValueType(type);
    }
}

int
HydratingHeMoMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return HeMoTKMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

int
HydratingHeMoMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_HydrationDegree ) {
        return 1;
    } else {
        return HeMoTKMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}

MaterialStatus *
HydratingHeMoMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new HydratingTransportMaterialStatus(1, domain, gp);
}
} // end namespace oofem
