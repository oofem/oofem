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

IsoInterfaceDamageMaterial_2 :: IsoInterfaceDamageMaterial_2(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{
    maxOmega = 0.999999;
}


IsoInterfaceDamageMaterial_2 :: ~IsoInterfaceDamageMaterial_2() { }


void
IsoInterfaceDamageMaterial_2 :: giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    IsoInterfaceDamageMaterialStatus_2 *status = static_cast< IsoInterfaceDamageMaterialStatus_2 * >( this->giveStatus(gp) );
    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;

    // compute equivalent strain
    equivStrain = macbra( jump.at(1) );

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
        this->computeDamageParam(omega, tempKappa, jump, gp);
    }

    this->give3dStiffnessMatrix_Eng(de, ElasticStiffness, gp, tStep);
    // damage in tension only
    if ( equivStrain >= 0.0 ) {
        de.times(1.0 - omega);
    }

    answer.beProductOf(de, jump);

    // update gp
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
}


void
IsoInterfaceDamageMaterial_2 :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    double om, un;
    IsoInterfaceDamageMaterialStatus_2 *status = static_cast< IsoInterfaceDamageMaterialStatus_2 * >( this->giveStatus(gp) );

    // assemble eleastic stiffness
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = kn;
    answer.at(2, 2) = ks;
    answer.at(3, 3) = ks;

    if ( rMode == ElasticStiffness ) {
        return;
    }

    if ( rMode == SecantStiffness ) {
        // Secant stiffness
        om = status->giveTempDamage();
        un = status->giveTempJump().at(1);
        om = min(om, maxOmega);
        // damage in tension only
        if ( un >= 0 ) {
            answer.times(1.0 - om);
        }

        return;
    } else {
        // Tangent Stiffness
        FloatArray se;
        const FloatArray &e = status->giveTempJump();
        se.beProductOf(answer, e);

        om = status->giveTempDamage();
        un = status->giveTempJump().at(1);
        om = min(om, maxOmega);
        // damage in tension only
        if ( un >= 0 ) {
            answer.times(1.0 - om);
            return;
#if 0
            // Unreachable code - commented out to supress compiler warnings
            double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
            if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                answer.at(1, 1) -= se.at(1) * dom;
                answer.at(2, 1) -= se.at(2) * dom;
            }
#endif
        }
    }
}


int
IsoInterfaceDamageMaterial_2 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IsoInterfaceDamageMaterialStatus_2 *status = static_cast< IsoInterfaceDamageMaterialStatus_2 * >( this->giveStatus(gp) );
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
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


IRResultType
IsoInterfaceDamageMaterial_2 :: initializeFrom(InputRecord *ir)
{
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

    // Parse the table file
    IR_GIVE_FIELD(ir, tablename, _IFT_IsoInterfaceDamageMaterial_2_tablename);

    is.open(tablename.c_str(), std :: ifstream :: in);
    if ( !is.is_open() ) {
        OOFEM_ERROR("Can't open table file %s.", tablename.c_str() );
    }

    // Read first line
    if ( is >> nbrOfLinesToRead ) {
        printf("NumberofLinestoRead: %d\n", nbrOfLinesToRead);
    } else {
        OOFEM_ERROR("Error reading table file, first line should be "
                    "an integer stating how many strain damage pairs that exist in the file.");
    }

    damages.resize(nbrOfLinesToRead + 1);
    strains.resize(nbrOfLinesToRead + 1);

    // Insert a (0,0) pair
    strains(0) = damages(0) = 0.0;

    for ( int i = 0; i < nbrOfLinesToRead; i++ ) {
        if ( !( is >> strains(i + 1) >> damages(i + 1) ) ) {
            OOFEM_ERROR("Error reading table file at line %d, expected a "
                         "strain damage pair.", i + 2);
        }

        if ( ( damages(i + 1) < damages(i) ) || ( strains(i + 1) < strains(i) ) ) {
            OOFEM_ERROR("Error reading table file at line %d, strain "
                         "and damage must be given in an increasing order and be positive.", i + 2);
        }
    }

    // We add e0 to the strains since strains should be given as the increase in
    // strain relative to e0.
    strains.add(e0);

    return StructuralInterfaceMaterial :: initializeFrom(ir);
}


void
IsoInterfaceDamageMaterial_2 :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(this->kn, _IFT_IsoInterfaceDamageMaterial_2_kn);
    input.setField(this->ks, _IFT_IsoInterfaceDamageMaterial_2_ks);
    input.setField(this->ft, _IFT_IsoInterfaceDamageMaterial_2_ft);
    input.setField(this->tablename, _IFT_IsoInterfaceDamageMaterial_2_tablename);
    input.setField(this->maxOmega, _IFT_IsoInterfaceDamageMaterial_2_maxOmega);
}


void
IsoInterfaceDamageMaterial_2 :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    kappa = macbra( strain.at(1) );
}

void
IsoInterfaceDamageMaterial_2 :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    if ( kappa > this->e0 ) {
        // Linear interpolation between table values.

        // If out of bounds damage is set to the last given damage value in the table
        if ( kappa >= strains.at( strains.giveSize() ) ) {
            omega = damages.at( damages.giveSize() );
        } else {
            // std::lower_bound uses binary search to find index with value bounding kappa from above
            int index = (int)(std :: lower_bound(strains.givePointer(), strains.givePointer() + strains.giveSize(), kappa) - strains.givePointer());

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
            double x0 = strains(index - 1);
            double x1 = strains(index );
            double y0 = damages(index - 1);
            double y1 = damages(index );

            // Interpolation formula
            omega = y0 + ( y1 - y0 ) * ( kappa - x0 ) / ( x1 - x0 );
        }
    } else {
        omega = 0.0;
    }
}


IsoInterfaceDamageMaterialStatus_2 :: IsoInterfaceDamageMaterialStatus_2(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
    damage = tempDamage = 0.0;
}


IsoInterfaceDamageMaterialStatus_2 :: ~IsoInterfaceDamageMaterialStatus_2()
{ }


void
IsoInterfaceDamageMaterialStatus_2 :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
        fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
}


void
IsoInterfaceDamageMaterialStatus_2 :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempDamage = this->damage;
}

void
IsoInterfaceDamageMaterialStatus_2 :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
}


contextIOResultType
IsoInterfaceDamageMaterialStatus_2 :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IsoInterfaceDamageMaterialStatus_2 :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
} // end namespace oofem
