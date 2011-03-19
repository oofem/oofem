/* $Header: /home/cvs/bp/oofem/sm/src/idmnl1.C,v 1.7.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *           Copyright (C) 2010 Christian Hoover, Vit Smilauer
 *
 *
 *
 *   Czech Technical University, Faculty of Civil Engineering,
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

// file: idm2.C


#include "idm2.h"
#include "structuralcrosssection.h"

#ifndef __MAKEDEPEND
 #include <math.h>
#endif


namespace oofem {
IsotropicDamageMaterial2 :: IsotropicDamageMaterial2(int n, Domain *d) : IsotropicDamageMaterial1(n, d)
{ }


IsotropicDamageMaterial2 :: ~IsotropicDamageMaterial2()
{ }

IRResultType
IsotropicDamageMaterial2 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int equivStrainType;
    //
    //eps0 - maximum elastic strain
    //gf - Fracture energy to propogate the crack on unit
    //epsf - maximum strain at zero stress
    //
    IsotropicDamageMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    // This is a switch to determine what type of softening the user would like.
    //  Type == 0, (default) exponential softening is used
    //  Type == 1, linear softening is used
    //  Type == 2, bilinear softening is used
    SofteningType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, SofteningType, IFT_IsotropicDamageMaterial2_softeningtype, "stype"); // Macro

    IR_GIVE_FIELD(ir, eps0, IFT_IsotropicDamageMaterial2_eps0, "eps0"); // Macro
    // Gf is the fracture energy required to advance the crack one unit
    IR_GIVE_FIELD(ir, gf, IFT_IsotropicDamageMaterial2_gf, "gf"); // Macro
    if ( SofteningType == 2 ) {
        // epsf is for the bilinear law, and corresponds to the strain at zero stress
        IR_GIVE_FIELD(ir, epsk, IFT_IsotropicDamageMaterial2_epsk, "epsk"); // Macro
        // Gft is for the bilinear law, and corresponds to the total energy required to fail the specimen
        IR_GIVE_FIELD(ir, gft, IFT_IsotropicDamageMaterial2_gft, "gft"); // Macro
    }

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


void
IsotropicDamageMaterial2 :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    double epsf, sigmak;
    double wf;
    double R, Lhs, E, Ft, help;
    int nite;
    if ( kappa > eps0 ) {
        E = this->giveLinearElasticMaterial()->give('E', gp);
        IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
        wf = 2 * gf / E / eps0;
        ef = wf / status->giveLe();
        if ( ef < eps0 ) {
            OOFEM_WARNING5("ef %f < eps0 %f, this leads to material snapback in element %d, characteristic length %f", gp->giveElement()->giveNumber(), ef, eps0, status->giveLe() );
            OOFEM_WARNING4("Material number %d, Decrease eps0, or increase Gf from %f to Gf=%f", this->giveNumber(), gf, E*eps0*E*eps0*status->giveLe()/2./E );
        }

        switch ( this->SofteningType ) {
        case ( 0 ):
            // This is the origonal code that was used before I replaced it
            nite = 0;
            //
            // iteration to achieve objectivity
            // we are finding state, where elastic stress is equal to
            // stress from crack-opening relation (ef = wf characterizes the carc opening diagram)
            omega = 0.0;
            Ft = E * eps0;
            do {
                nite++;
                help = omega * kappa / ef;
                R = ( 1. - omega ) * E * kappa - Ft *exp(-help);
                Lhs = E * kappa - Ft *exp(-help) * kappa / ef;
                omega += R / Lhs;
                if ( nite > 40 ) {
                    _error("computeDamageParam: algorithm not converging");
                }
            } while ( fabs(R) >= IDM1_ITERATION_LIMIT );

            if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
                _error("computeDamageParam: damage parameter out of range, snap-back problems");
            }

            break;
        case ( 1 ):
            // This is the linear damage law
            if ( kappa < ef ) {
                omega = ( 1.0 - ( eps0 / kappa ) * ( ef - kappa ) / ( ef - eps0 ) );
            } else if ( kappa <= eps0 ) {
                omega = 0.0;
            } else {
                omega = maxOmega;
            }

            break;
        case ( 2 ):
            // This is the bi-linear damage law
            sigmak = E * eps0 * ( ef - epsk ) / ( ef - eps0 );
            epsf = 2 * ( gft - gf ) / sigmak / status->giveLe() + ef;
            if ( ( epsk > ef ) || ( epsk < eps0 ) ) {
                OOFEM_WARNING4("epsk %f is not between ef %f and eps0 %f", epsk, ef, eps0);
            }

            if ( gft < gf ) {
                OOFEM_WARNING3("the total fracture energy gft %f must be greater than the initial fracture energy gf %f", gft, gf);
            }

            if ( kappa <= epsk ) {
                omega = 1.0 - ( ( eps0 / kappa ) * ( epsk - kappa ) / ( epsk - eps0 ) + ( ( sigmak / ( E * kappa ) ) * ( kappa - eps0 ) / ( epsk - eps0 ) ) );
            } else if ( kappa > epsk && kappa <= epsf ) {
                omega = 1.0 - ( ( sigmak / ( E * kappa ) ) * ( epsf - kappa ) / ( epsf - epsk ) );
            } else if ( kappa <= e0 ) {
                omega = 0.0;
            } else {
                omega = maxOmega;
            }

            break;

        default:
            OOFEM_ERROR1("Unknown softening type");
        }
    } else {
        omega = 0.;
    }
}

void
IsotropicDamageMaterial2 :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp)
{
    int i, indx = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3), fullstrain;
    FloatMatrix principalDir(3, 3);
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    // const double eps0 = this->give(eps0_ID, gp);
    // const double gf = this->give(gf_ID, gp);
    double wf, E;
    E = this->giveLinearElasticMaterial()->give('E', gp);
    wf = 2 * gf / eps0 / E;
    // const double ef = this->give(ef_ID, gp);


    crossSection->giveFullCharacteristicVector(fullstrain, gp, strainVector);

    if ( ( kappa > eps0 ) && ( status->giveDamage() == 0. ) ) {
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

        // CGH 2/2/2010 This was changed for the nonlocal project
        // to ensure the maximum principle strain was always in the x-direction
        // This ensures that opening will always happen in the x-direction so that Gf

        crackPlaneNormal.at(1) = 1;
        crackPlaneNormal.at(2) = 0;
        crackPlaneNormal.at(3) = 0;
        le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
        // remember le in cooresponding status
        status->setLe(le);

        if ( eps0 >= wf / le ) {
            _warning2("instanciateFrom: suspicious eps0>=wf/le", 1);
        }
    }
}

// double
// IsotropicDamageMaterial2 :: give(int aProperty, GaussPoint* gp)
// {
//   double answer;
//   if (RandomMaterialExtensionInterface::give(aProperty, gp, answer)) {
//     return answer;
//   } else if (aProperty == eps0_ID) {
//     return this->eps0;
//   } else if (aProperty == gf_ID) {
//     return this->gf;
//     // } else if (aProperty == ef2_ID) {
//     //   return this->ef2;
//     // } else if (aProperty == sigmak_ID) {
//     //   return this->sigmak;
//   } else {
//     return IsotropicDamageMaterial::give(aProperty, gp);
//   }
// }
} // end namespace oofem
