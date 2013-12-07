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

#include <algorithm>

#include "isointerfacedamage02.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IsoInterfaceDamageMaterial_2);

IsoInterfaceDamageMaterial_2 :: IsoInterfaceDamageMaterial_2(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    maxOmega = 0.999999;
}


IsoInterfaceDamageMaterial_2 :: ~IsoInterfaceDamageMaterial_2()
//
// destructor
//
{ }

int
IsoInterfaceDamageMaterial_2 :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _2dInterface || mode == _3dInterface;
}


void
IsoInterfaceDamageMaterial_2 :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
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
IsoInterfaceDamageMaterial_2 :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                     const FloatArray &totalStrain,
                                                     TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    IsoInterfaceDamageMaterialStatus_2 *status = static_cast< IsoInterfaceDamageMaterialStatus_2 * >( this->giveStatus(gp) );
    FloatArray strainVector, reducedTotalStrainVector;
    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;

    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);

    //crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, atTime);

    // compute value of loading function if strainLevel crit apply
    f = equivStrain - status->giveKappa();

    if ( f <= 0.0 ) {
        // damage do not grow
        tempKappa = status->giveKappa();
        omega     = status->giveDamage();
    } else {
        // damage grows
        tempKappa = equivStrain;
        // evaluate damage parameter
        this->computeDamageParam(omega, tempKappa, reducedTotalStrainVector, gp);
    }

    this->giveStiffnessMatrix(de, ElasticStiffness, gp, atTime);
    // damage in tension only
    if ( equivStrain >= 0.0 ) {
        de.times(1.0 - omega);
    }

    answer.beProductOf(de, reducedTotalStrainVector);

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
}

void
IsoInterfaceDamageMaterial_2 :: giveStiffnessMatrix(FloatMatrix &answer,
                                                    MatResponseMode rMode,
                                                    GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _2dInterface:
        give2dInterfaceMaterialStiffnessMatrix(answer, rMode, gp, atTime);
        break;
    case _3dInterface:
        give3dInterfaceMaterialStiffnessMatrix(answer, rMode, gp, atTime);
        break;
    default:
        StructuralMaterial :: giveStiffnessMatrix(answer, rMode, gp, atTime);
    }
}


void
IsoInterfaceDamageMaterial_2 :: give2dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                                       GaussPoint *gp, TimeStep *atTime)
{
    double om, un;
    IsoInterfaceDamageMaterialStatus_2 *status = static_cast< IsoInterfaceDamageMaterialStatus_2 * >( this->giveStatus(gp) );


    if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
        // assemble eleastic stiffness
        answer.resize(2, 2);
        answer.at(1, 1) = kn;
        answer.at(2, 2) = ks;
        answer.at(1, 2) = answer.at(2, 1) = 0.0;

        if ( rMode == ElasticStiffness ) {
            return;
        }

        if ( rMode == SecantStiffness ) {
            // Secant stiffness
            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
            }

            return;
        } else {
            // Tangent Stiffness
            FloatArray se(2), e(2);
            e = status->giveTempStrainVector();
            se.beProductOf(answer, e);

            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
                return;

                /* Unreachable code - commented out to supress compiler warnings
                 * double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
                 * if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                 *  answer.at(1, 1) -= se.at(1) * dom;
                 *  answer.at(2, 1) -= se.at(2) * dom;
                 * }
                 */
            }
        }
    }  else {
        _error("give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode");
    }
}


void
IsoInterfaceDamageMaterial_2 :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                                       GaussPoint *gp, TimeStep *atTime)
{
    double om, un;
    IsoInterfaceDamageMaterialStatus_2 *status = static_cast< IsoInterfaceDamageMaterialStatus_2 * >( this->giveStatus(gp) );


    if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
        // assemble eleastic stiffness
        answer.resize(3, 3);
        answer.at(1, 1) = kn;
        answer.at(2, 2) = ks;
        answer.at(3, 3) = ks;
        answer.at(1, 2) = answer.at(2, 1) = answer.at(1, 3) = answer.at(3, 1) = answer.at(2, 3) = answer.at(3, 2) = 0.0;

        if ( rMode == ElasticStiffness ) {
            return;
        }

        if ( rMode == SecantStiffness ) {
            // Secant stiffness
            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
            }

            return;
        } else {
            // Tangent Stiffness
            FloatArray se, e;
            e = status->giveTempStrainVector();
            se.beProductOf(answer, e);

            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
                return;
                /* Unreachable code - commented out to supress compiler warnings
                 * double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
                 * if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                 *  answer.at(1, 1) -= se.at(1) * dom;
                 *  answer.at(2, 1) -= se.at(2) * dom;
                 * }
                 */
            }
        }
    }  else {
        _error("give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode");
    }
}


int
IsoInterfaceDamageMaterial_2 :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    IsoInterfaceDamageMaterialStatus_2 *status = static_cast< IsoInterfaceDamageMaterialStatus_2 * >( this->giveStatus(aGaussPoint) );
    if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveTempDamage();
        return 1;
    } else if ( type == IST_PrincipalDamageTensor ) {
        answer.resize(3);
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_PrincipalDamageTempTensor ) {
        answer.resize(3);
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveTempDamage();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

void
IsoInterfaceDamageMaterial_2 :: giveThermalDilatationVector(FloatArray &answer,
                                                            GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    answer.resize(6);
    answer.zero();
    answer.at(1) = this->tempDillatCoeff;
    answer.at(2) = this->tempDillatCoeff;
    answer.at(3) = this->tempDillatCoeff;
}


IRResultType
IsoInterfaceDamageMaterial_2 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    std :: ifstream is;
    int nbrOfLinesToRead;

    IR_GIVE_FIELD(ir, kn, _IFT_IsoInterfaceDamageMaterial_2_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_IsoInterfaceDamageMaterial_2_ks);
    IR_GIVE_FIELD(ir, ft, _IFT_IsoInterfaceDamageMaterial_2_ft);

    this->e0 = ft / kn;

    //Set limit on the maximum isotropic damage parameter if needed
    IR_GIVE_OPTIONAL_FIELD(ir, maxOmega, _IFT_IsoInterfaceDamageMaterial_2_maxOmega);
    maxOmega = min(maxOmega, 0.999999);
    maxOmega = max(maxOmega, 0.0);

    IR_GIVE_FIELD(ir, tempDillatCoeff, _IFT_IsoInterfaceDamageMaterial_2_talpha);

    // Parse the table file
    IR_GIVE_FIELD(ir, tablename, _IFT_IsoInterfaceDamageMaterial_2_tablename);

    is.open(tablename.c_str(), std :: ifstream :: in);
    if ( !is.is_open() ) {
        OOFEM_ERROR2( "IsoInterfaceDamageMaterial :: initializeFrom: Can't open table file %s.", tablename.c_str() );
    }

    // Read first line
    if ( is >> nbrOfLinesToRead ) {
        printf("NumberofLinestoRead: %d\n", nbrOfLinesToRead);
    } else {
        OOFEM_ERROR("IsoInterfaceDamageMaterial :: initializeFrom: Error reading table file, first line should be "
                    "an integer stating how many strain damage pairs that exist in the file.");
    }

    damages.resize(nbrOfLinesToRead + 1);
    strains.resize(nbrOfLinesToRead + 1);

    // Insert a (0,0) pair
    strains(0) = damages(0) = 0.0;

    for ( int i = 0; i < nbrOfLinesToRead; i++ ) {
        if ( !( is >> strains(i + 1) >> damages(i + 1) ) ) {
            OOFEM_ERROR2("IsoInterfaceDamageMaterial :: initializeFrom: Error reading table file at line %d, expected a "
                         "strain damage pair.", i + 2);
        }

        if ( ( damages(i + 1) < damages(i) ) || ( strains(i + 1) < strains(i) ) ) {
            OOFEM_ERROR2("IsoInterfaceDamageMaterial :: initializeFrom: Error reading table file at line %d, strain "
                         "and damage must be given in an increasing order and be positive.", i + 2);
        }
    }

    // We add e0 to the strains since strains should be given as the increase in
    // strain relative to e0.
    strains.add(e0);

    return StructuralMaterial :: initializeFrom(ir);
}


void
IsoInterfaceDamageMaterial_2 :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->kn, _IFT_IsoInterfaceDamageMaterial_2_kn);
    input.setField(this->ks, _IFT_IsoInterfaceDamageMaterial_2_ks);
    input.setField(this->ft, _IFT_IsoInterfaceDamageMaterial_2_ft);
    input.setField(this->tablename, _IFT_IsoInterfaceDamageMaterial_2_tablename);
    input.setField(this->maxOmega, _IFT_IsoInterfaceDamageMaterial_2_maxOmega);
    input.setField(this->tempDillatCoeff, _IFT_IsoInterfaceDamageMaterial_2_talpha);
}


void
IsoInterfaceDamageMaterial_2 :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    kappa = macbra( strain.at(1) );
}

void
IsoInterfaceDamageMaterial_2 :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    if ( kappa > this->e0 ) {
        // Linear interpolation between able values.

        // If out of bounds damage is set to the last given damage value in the table
        if ( kappa >= strains.at( strains.giveSize() ) ) {
            omega = damages.at( damages.giveSize() );
        } else {
            // std::lower_bound uses binary search to find index with value bounding kappa from above
            double *index = std :: lower_bound(strains.givePointer(), strains.givePointer() + strains.giveSize(), kappa);

#if 0
            printf("e0 %lf\n", e0);
            printf( "sizeof %d\n", sizeof( strains.givePointer() ) );
            printf("pointer: %d\n", index);
            printf("value: %lf\n", * index);
            printf( "Index found: %d\n", index - strains.givePointer() );
            printf( "First index: %d\n", strains.givePointer() );
            printf( "Last index: %d\n", strains.givePointer() + strains.giveSize() );
#endif

            // Pointer arithmetic to find the values used in interpolation
            double x0 = strains(index - strains.givePointer() - 1);
            double x1 = strains( index - strains.givePointer() );
            double y0 = damages(index - strains.givePointer() - 1);
            double y1 = damages( index - strains.givePointer() );

            // Interpolation formula
            omega = y0 + ( y1 - y0 ) * ( kappa - x0 ) / ( x1 - x0 );
        }
    } else {
        omega = 0.0;
    }
}


IsoInterfaceDamageMaterialStatus_2 :: IsoInterfaceDamageMaterialStatus_2(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
    damage = tempDamage = 0.0;
}


IsoInterfaceDamageMaterialStatus_2 :: ~IsoInterfaceDamageMaterialStatus_2()
{ }


void
IsoInterfaceDamageMaterialStatus_2 :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
        fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
}


void
IsoInterfaceDamageMaterialStatus_2 :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempDamage = this->damage;
}

void
IsoInterfaceDamageMaterialStatus_2 :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
}


contextIOResultType
IsoInterfaceDamageMaterialStatus_2 :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IsoInterfaceDamageMaterialStatus_2 :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
} // end namespace oofem
