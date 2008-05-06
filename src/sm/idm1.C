/* $Header: /home/cvs/bp/oofem/sm/src/idm1.C,v 1.8.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

// file: idm1.C


#include "idm1.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mmaclosestiptransfer.h"
#include "datastream.h"
#include "contextioerr.h"

#ifndef __MAKEDEPEND
#include <math.h>
#endif


#ifdef IDM_USE_MMAClosestIPTransfer
MMAClosestIPTransfer IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMAShapeFunctProjection
MMAShapeFunctProjection IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMALeastSquareProjection
MMALeastSquareProjection IsotropicDamageMaterial1 :: mapper;
#endif

IsotropicDamageMaterial1 :: IsotropicDamageMaterial1(int n, Domain *d) : IsotropicDamageMaterial(n, d)
    //
    // constructor
    //
{
    // deleted by paren, where linearElasticMaterial instance declared
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}


IsotropicDamageMaterial1 :: ~IsotropicDamageMaterial1()
//
// destructor
//
{ }

IRResultType
IsotropicDamageMaterial1 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int equivStrainType;
    IsotropicDamageMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    IR_GIVE_FIELD(ir, e0, IFT_IsotropicDamageMaterial1_e0, "e0"); // Macro
    // in ef variable the wf (crack opening) is stored.
    IR_GIVE_FIELD(ir, ef, IFT_IsotropicDamageMaterial1_ef, "ef"); // Macro
    //if (e0 >= ef/le) _warning2 ("instanciateFrom: suspicious e0>=ef/le", 1);

    equivStrainType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, equivStrainType, IFT_IsotropicDamageMaterial1_equivstraintype, "equivstraintype"); // Macro
    if ( equivStrainType == 1 ) {
        this->equivStrainType = EST_Rankine;
    } else if ( equivStrainType == 2 ) {
        this->equivStrainType = EST_ElasticEnergy;
    } else {
        this->equivStrainType = EST_Mazars; // default
    }

    this->mapper.initializeFrom(ir);

    return IRRT_OK;
}


int
IsotropicDamageMaterial1 :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    IsotropicDamageMaterial :: giveInputRecordString(str, keyword);
    linearElasticMaterial->giveInputRecordString(str, false);
    sprintf(buff, " e0 %e ef %e equivstraintype %d", this->e0, this->ef, ( int ) this->equivStrainType);
    str += buff;

    return 1;
}




void
IsotropicDamageMaterial1 :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }

    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains, fullstrain;

        crossSection->giveFullCharacteristicVector(fullstrain, gp, strain);

        // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            double nu = lmat->give(NYxz);
            fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            double nu = lmat->give(NYxz);
            fullstrain.at(2) = -nu *fullstrain.at(1);
            fullstrain.at(3) = -nu *fullstrain.at(1);
        }

        this->computePrincipalValues(principalStrains, fullstrain, principal_strain);

        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        kappa = sqrt(posNorm);
    } else if ( this->equivStrainType == EST_Rankine ) {
        // EST_Rankine equiv strain measure
        int i;
        FloatMatrix de;
        FloatArray stress, fullStress, principalStress;
        double sum = 0.;

        lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
        stress.beProductOf(de, strain);
        crossSection->giveFullCharacteristicVector(fullStress, gp, stress);
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        for ( i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                sum += principalStress.at(i) * principalStress.at(i);
            }
        }

        kappa = sqrt(sum) / lmat->give('E');
    } else if ( this->equivStrainType == EST_ElasticEnergy ) {
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
        stress.beProductOf(de, strain);
        sum = dotProduct( strain, stress, strain.giveSize() );

        kappa = sqrt( sum / lmat->give('E') );
    } else {
        _error("computeEquivalentStrain: unknown EquivStrainType");
    }
}

void
IsotropicDamageMaterial1 :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    if ( kappa > this->e0 ) {
        // omega = 1.0-(this->e0/kappa)*exp(-(kappa-this->e0)/(this->ef-this->e0));
        IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
        int nite = 0;
        double R, Lhs, E, Ft, help;

        //
        // iteration to achieve objectivity
        // we are finding state, where elastic stress is equal to
        // stress from crack-opening relation (ef = wf characterizes the carc opening diagram)

        omega = 0.0;
        E = this->giveLinearElasticMaterial()->give('E');
        Ft = E * this->e0;
        //printf ("\nle=%f, kappa=%f", status->giveLe(), kappa);
        do {
            nite++;
            help = status->giveLe() * omega * kappa / this->ef;
            R = ( 1. - omega ) * E * kappa - Ft *exp(-help);
            Lhs = E * kappa - Ft *exp(-help) * status->giveLe() * kappa / this->ef;
            omega += R / Lhs;
            // printf ("\n%d: R=%f, omega=%f",nite, R, omega);
            if ( nite > 40 ) {
                _error("computeDamageParam: algorithm not converging");
            }
        } while ( fabs(R) >= IDM1_ITERATION_LIMIT );

        if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
            _error("computeDamageParam: internal error");
        }
    } else {
        omega = 0.0;
    }
}


void
IsotropicDamageMaterial1 :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp)
{
    int i, indx = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3), fullstrain;
    FloatMatrix principalDir(3, 3);
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    crossSection->giveFullCharacteristicVector(fullstrain, gp, strainVector);

    if ( ( kappa > this->e0 ) && ( status->giveDamage() == 0. ) ) {
        this->computePrincipalValDir(principalStrains, principalDir, fullstrain, principal_strain);
        // finfd index of max positive principal strain
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
        // remember le in cooresponding status
        status->setLe(le);

        if ( e0 >= ef / le ) {
            _warning2("instanciateFrom: suspicious e0>=ef/le", 1);
        }
    }
}


Interface *
IsotropicDamageMaterial1 :: giveInterface(InterfaceType type)
{
    if ( type == MaterialModelMapperInterfaceType ) {
        return ( MaterialModelMapperInterface * ) this;
    } else {
        return NULL;
    }
}

int
IsotropicDamageMaterial1 :: MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep)
{
    int result;
    FloatArray intVal, strainIncr(3);
    IntArray toMap(3);
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);


    toMap.at(1) = ( int ) IST_MaxEquivalentStrainLevel;
    toMap.at(2) = ( int ) IST_DamageTensor;
    toMap.at(3) = ( int ) IST_StrainTensor;
    this->mapper.init(oldd, toMap, gp, tStep);

    result = mapper.mapVariable(intVal, gp, IST_MaxEquivalentStrainLevel, tStep);
    if ( result ) {
        status->setTempKappa( intVal.at(1) );
    }

    result = mapper.mapVariable(intVal, gp, IST_DamageTensor, tStep);
    if ( result ) {
        status->setTempDamage( intVal.at(1) );
    }

#ifdef IDM_USE_MAPPEDSTRAIN
    result = mapper.mapVariable(intVal, gp, IST_StrainTensor, tStep);
    if ( result ) {
        status->letTempStrainVectorBe(intVal);
    }

#endif
    status->updateYourself(tStep);

    return result;
}




int
IsotropicDamageMaterial1 :: MMI_update(GaussPoint *gp,  TimeStep *tStep, FloatArray *estrain)
{
    int result = 1;
    FloatArray intVal, strain;
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);

    // now update all internal vars accordingly
    strain = status->giveStrainVector();
#ifdef IDM_USE_MAPPEDSTRAIN
    this->giveRealStressVector(intVal, ReducedForm, gp, strain, tStep);
#else
    this->giveRealStressVector(intVal, ReducedForm, gp, * estrain, tStep);
#endif
    this->updateYourself(gp, tStep);
    return result;
}


int
IsotropicDamageMaterial1 :: MMI_finish(TimeStep *tStep)
{
    this->mapper.finish(tStep);
    return 1;
}


IsotropicDamageMaterial1Status :: IsotropicDamageMaterial1Status(int n, Domain *d, GaussPoint *g) :
    IsotropicDamageMaterialStatus(n, d, g)
{
    le = 0.0;
}

void
IsotropicDamageMaterial1Status :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterialStatus :: initTempStatus();
}



void
IsotropicDamageMaterial1Status :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterialStatus :: updateYourself(atTime);
}



contextIOResultType
IsotropicDamageMaterial1Status :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IsotropicDamageMaterial1Status :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

