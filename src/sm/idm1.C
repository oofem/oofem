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

#include "idm1.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"

namespace oofem {
#ifdef IDM_USE_MMAClosestIPTransfer
MMAClosestIPTransfer IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMAContainingElementProjection
MMAContainingElementProjection IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMAShapeFunctProjection
MMAShapeFunctProjection IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMALeastSquareProjection
MMALeastSquareProjection IsotropicDamageMaterial1 :: mapper;
#endif

IsotropicDamageMaterial1 :: IsotropicDamageMaterial1(int n, Domain *d) : IsotropicDamageMaterial(n, d),
    RandomMaterialExtensionInterface()
    //
    // constructor
    //
{
    // deleted by parent, where linearElasticMaterial instance declared
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    equivStrainType = EST_Unknown;
    softType = ST_Unknown;
    k = 0.;
    md = 1.;
    damageLaw = 0;
    e0 = 0.;
    ef = 0.;
    wf = 0.;
    gf = 0.;
    ecsMethod = ECSM_Unknown;
}


IsotropicDamageMaterial1 :: ~IsotropicDamageMaterial1()
//
// destructor
//
{ }

IRResultType
IsotropicDamageMaterial1 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";     // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int equivStrainType;
    IsotropicDamageMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    // specify the type of formula for equivalent strain
    IR_GIVE_OPTIONAL_FIELD(ir, equivStrainType, IFT_IsotropicDamageMaterial1_equivstraintype, "equivstraintype");     // Macro
    if ( equivStrainType == 1 ) {
        this->equivStrainType = EST_Rankine_Smooth;
    } else if ( equivStrainType == 2 ) {
        this->equivStrainType = EST_ElasticEnergy;
    } else if ( equivStrainType == 3 ) {
        this->equivStrainType = EST_Mises;
        IR_GIVE_FIELD(ir, k, IFT_IsotropicDamageMaterial1_k, "k");
    } else if ( equivStrainType == 4 ) {
        this->equivStrainType = EST_Rankine_Standard;
    } else if ( equivStrainType == 5 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStress;
    } else if ( equivStrainType == 6 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStrain;
    } else {
        this->equivStrainType = EST_Mazars;     // default
    }

    // specify the type of formula for damage evolution law
    IR_GIVE_OPTIONAL_FIELD(ir, damageLaw, IFT_IsotropicDamageMaterial1_softeningtype, "damlaw");
    if ( damageLaw != 7 ) {
      IR_GIVE_FIELD(ir, e0, IFT_IsotropicDamageMaterial1_e0, "e0");
    }
    //applies only in this class
    if ( (giveClassID() == IsotropicDamageMaterial1Class) ||   (giveClassID() == IDNLMaterialClass) ) {
        switch ( damageLaw ) {
        case 0:     // exponential softening - default
            if ( ir->hasField(IFT_IsotropicDamageMaterial1_wf, "wf") ) {
                this->softType = ST_Exponential_Cohesive_Crack;
                IR_GIVE_FIELD(ir, wf, IFT_IsotropicDamageMaterial1_wf, "wf");
            } else if ( ir->hasField(IFT_IsotropicDamageMaterial1_gf, "gf") ) {
                this->softType = ST_Exponential_Cohesive_Crack;
                IR_GIVE_FIELD(ir, gf, IFT_IsotropicDamageMaterial1_ef, "gf");
            } else {
                this->softType = ST_Exponential;
                IR_GIVE_FIELD(ir, ef, IFT_IsotropicDamageMaterial1_ef, "ef");
            }

            break;
        case 1:     // linear softening law
            if ( ir->hasField(IFT_IsotropicDamageMaterial1_wf, "wf") ) {
                this->softType = ST_Linear_Cohesive_Crack;
                IR_GIVE_FIELD(ir, wf, IFT_IsotropicDamageMaterial1_wf, "wf");
            } else if ( ir->hasField(IFT_IsotropicDamageMaterial1_gf, "gf") ) {
                this->softType = ST_Linear_Cohesive_Crack;
                IR_GIVE_FIELD(ir, gf, IFT_IsotropicDamageMaterial1_ef, "gf");
            } else {
                this->softType = ST_Linear;
                IR_GIVE_FIELD(ir, ef, IFT_IsotropicDamageMaterial1_ef, "ef");
            }

            break;
        case 2:     // bilinear softening
            if ( ir->hasField(IFT_IsotropicDamageMaterial1_gf, "gf") ) {
                this->softType = ST_BiLinear_Cohesive_Crack;
                IR_GIVE_FIELD(ir, gf, IFT_IsotropicDamageMaterial1_ef, "gf");
                // ek is for the bilinear law, and corresponds to the strain at the knee point
                IR_GIVE_FIELD(ir, ek, IFT_IsotropicDamageMaterial1_ek, "ek");         // Macro
                // Gft is for the bilinear law, and corresponds to the total energy required to fail the specimen
                IR_GIVE_FIELD(ir, gft, IFT_IsotropicDamageMaterial1_gft, "gft");         // Macro
            } else {
                OOFEM_ERROR("Bilinear softening for wf and ef not implemented");
            }

            break;
        case 3: // reserved for Hordijk's law
            break;
        case 4:
            this->softType = ST_Mazars;
            IR_GIVE_FIELD(ir, At, IFT_IsotropicDamageMaterial1_At, "at");
            IR_GIVE_FIELD(ir, Bt, IFT_IsotropicDamageMaterial1_Bt, "bt");
            break;
        case 5:
            this->softType = ST_Smooth;
            md = 1.;
            IR_GIVE_OPTIONAL_FIELD(ir, md, IFT_IsotropicDamageMaterial1_md, "md");
            break;
        case 6:
            this->softType = ST_Disable_Damage;
            break;
        case 7:
            this->softType = ST_SmoothExtended;
            double E, ep, ft; // auxiliary input variables
            IR_GIVE_FIELD(ir, E,  IFT_IsotropicLinearElasticMaterial_e, "e");
            IR_GIVE_FIELD(ir, ep, IFT_IsotropicDamageMaterial1_ep, "ep");
            IR_GIVE_FIELD(ir, ft, IFT_IsotropicDamageMaterial1_ft, "ft");
            this->md = 1. / log(E * ep / ft);
            this->e0 = ep * pow(md, 1. / md);
            IR_GIVE_FIELD(ir, e1, IFT_IsotropicDamageMaterial1_e1, "e1");
            this->s1 = e1 * exp( -pow(e1 / e0, md) );
            this->ef = -e1 / ( 1. - md * pow(e1 / e0, md) );
            IR_GIVE_FIELD(ir, e2, IFT_IsotropicDamageMaterial1_e2, "e2");
            IR_GIVE_FIELD(ir, nd, IFT_IsotropicDamageMaterial1_nd, "nd");
            break;
        default:
            OOFEM_ERROR2("Softening type number %d is unknown", damageLaw);
        }
    }

    if ( ( softType == ST_Exponential_Cohesive_Crack ) || ( softType == ST_Linear_Cohesive_Crack ) || ( softType == ST_BiLinear_Cohesive_Crack ) ) {
      int ecsm = 0;
      IR_GIVE_OPTIONAL_FIELD(ir, ecsm, IFT_IsotropicDamageMaterial1_md, "ecsm");
      switch ( ecsm ) {
      case 1: ecsMethod = ECSM_SquareRootOfArea; break;
      case 2: ecsMethod = ECSM_ProjectionCentered; break;
      case 3: ecsMethod = ECSM_Oliver1; break;
      case 4: ecsMethod = ECSM_Oliver1modified; break;
      default: ecsMethod = ECSM_Projection;
      }
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
            double nu = lmat->give(NYxz, gp);
            fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            double nu = lmat->give(NYxz, gp);
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
    } else if ( ( this->equivStrainType == EST_Rankine_Smooth ) || ( this->equivStrainType == EST_Rankine_Standard ) ) {
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
                if ( this->equivStrainType == EST_Rankine_Smooth ) {
                    sum += principalStress.at(i) * principalStress.at(i);
                } else if ( sum < principalStress.at(i) )   {
                    sum = principalStress.at(i);
                }
            } else if ( sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
            }
        }
        if ( this->equivStrainType == EST_Rankine_Smooth ) {
            sum = sqrt(sum);
        }
        kappa = sum / lmat->give('E', gp);

    } else if ( ( this->equivStrainType == EST_ElasticEnergy ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStrain ) ){
        // equivalent strain expressions based on elastic energy
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
        if ( this->equivStrainType == EST_ElasticEnergy ) {
            // standard elastic energy
            stress.beProductOf(de, strain);
            sum = strain.dotProduct(stress);
        } else if ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) {
            // elastic energy corresponding to positive part of stress
            FloatArray fullStress, principalStress;
            crossSection->giveFullCharacteristicVector(fullStress, gp, stress);
            this->computePrincipalValues(principalStress, fullStress, principal_stress);
            // TO BE FINISHED

        } else {
            // elastic energy corresponding to positive part of strain
            // TO BE DONE
        }

        kappa = sqrt( sum / lmat->give('E', gp) );
    } else if ( this->equivStrainType == EST_Mises ) {
        double nu = lmat->give(NYxz, NULL);
        FloatArray principalStrains, fullstrain;
        crossSection->giveFullCharacteristicVector(fullstrain, gp, strain);
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            fullstrain.at(2) = -nu *fullstrain.at(1);
            fullstrain.at(3) = -nu *fullstrain.at(1);
        }

        this->computePrincipalValues(principalStrains, fullstrain, principal_strain);
        double I1e, J2e;
        this->computeStrainInvariants(principalStrains, I1e, J2e);
        double a, b, c;
        a = ( k - 1 ) * I1e / ( 2 * k * ( 1 - 2 * nu ) );
        b = ( k - 1 ) * ( k - 1 ) * I1e * I1e / ( ( 1 - 2 * nu ) * ( 1 - 2 * nu ) );
        c = 12 * k * J2e / ( ( 1 + nu ) * ( 1 + nu ) );
        kappa = a + 1 / ( 2 * k ) * sqrt(b + c);
    } else {
        _error("computeEquivalentStrain: unknown EquivStrainType");
    }
}

void
IsotropicDamageMaterial1 :: computeStrainInvariants(const FloatArray &strainVector, double &I1e, double &J2e)
{
    I1e = strainVector.at(1) + strainVector.at(2) + strainVector.at(3);
    double s1 = strainVector.at(1) * strainVector.at(1);
    double s2 = strainVector.at(2) * strainVector.at(2);
    double s3 = strainVector.at(3) * strainVector.at(3);
    J2e = 1. / 2. * ( s1 + s2 + s3 ) - 1. / 6. * ( I1e * I1e );
}

void
IsotropicDamageMaterial1 :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    if ( this->softType == ST_Disable_Damage ) { //dummy material with no damage
        omega = 0.;
    } else if ( this->softType == ST_Exponential_Cohesive_Crack || this->softType == ST_Linear_Cohesive_Crack || this->give(gf_ID, gp) != 0. ) { // adjustment of softening law according to the element size, given crack opening or fracture energy
        computeDamageParamForCohesiveCrack(omega, kappa, gp);
    } else { // no adjustment according to element size, given fracturing strain
        omega = damageFunction(kappa, gp);
    }
}

void
IsotropicDamageMaterial1 :: computeDamageParamForCohesiveCrack(double &omega, double kappa, GaussPoint *gp)
{
    const double e0 = this->give(e0_ID, gp);  // e0 is the strain at the peak stress
    const double E = this->giveLinearElasticMaterial()->give('E', gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);     // wf is the crack opening
    double Le;
    omega = 0.0;

    if ( kappa > e0 ) {
        if ( this->gf != 0. ) { //cohesive crack model
            if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
                wf = this->gf / E / e0; // wf is the crack opening
            } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi)linear softening law
                wf = 2. * gf / E / e0; // wf is the crack opening
            } else {
                OOFEM_ERROR2("Gf unsupported for softening type softType = %d", softType);
            }
        }

        IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
        Le = status->giveLe();
        ef = wf / Le;    //ef is the fracturing strain
        if ( ef < e0 ) { //check that no snapback occurs
            double minGf = 0.;
            OOFEM_WARNING5("ef %e < e0 %e, this leads to material snapback in element %d, characteristic length %f", ef, e0, gp->giveElement()->giveNumber(), Le);
            if ( gf != 0. ) { //cohesive crack
                if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
                    minGf = E * e0 * e0 * Le;
                } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi)linear softening law
                    minGf = E * e0 * e0 * Le / 2.;
                } else {
                    OOFEM_WARNING2("Gf unsupported for softening type softType = %d", softType);
                }

                OOFEM_ERROR4("Material number %d, decrease e0, or increase Gf from %f to Gf=%f", this->giveNumber(), gf, minGf);
            }
            _error("\n");
        }

        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            if ( kappa < ef ) {
                omega = ( ef / kappa ) * ( kappa - e0 ) / ( ef - e0 );
            } else {
                omega = 1.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
            }
        } else if (  this->softType == ST_BiLinear_Cohesive_Crack ) {
            double ek = this->give(ek_ID, gp);
            double gft = this->give(gft_ID, gp);
            double ef = 2 * gf / E / e0 / Le; //the first part corresponds to linear softening
            double sigmak = E * e0 * ( ef - ek ) / ( ef - e0 );
            double epsf = 2 * ( gft - gf ) / sigmak / Le + ef;
            if ( ( ek > ef ) || ( ek < e0 ) ) {
                OOFEM_WARNING4("ek %f is not between e0 %f and ef %f", ek, e0, ef);
            }

            if ( gft < gf ) {
                OOFEM_ERROR3("The total fracture energy gft %f must be greater than the initial fracture energy gf %f", gft, gf);
            }

            if ( kappa <= ek ) {
                omega = 1.0 - ( ( e0 / kappa ) * ( ek - kappa ) / ( ek - e0 ) + ( ( sigmak / ( E * kappa ) ) * ( kappa - e0 ) / ( ek - e0 ) ) );
            } else if ( kappa > ek && kappa <= epsf ) {
                omega = 1.0 - ( ( sigmak / ( E * kappa ) ) * ( epsf - kappa ) / ( epsf - ek ) );
            } else if ( kappa <= e0 ) {
                omega = 0.0;
            } else {
                omega = maxOmega;
            }
        } else if (  this->softType == ST_Exponential_Cohesive_Crack ) {
            // exponential cohesive crack - iteration needed
            double R, Lhs, help;
            int nite = 0;
            // iteration to achieve objectivity
            // we are looking for a state in which the elastic stress is equal to
            // the stress from crack-opening relation
            // ef has now the meaning of strain
            do {
                nite++;
                help = omega * kappa / ef;
                R = ( 1. - omega ) * kappa - e0 *exp(-help);
                Lhs = kappa - e0 *exp(-help) * kappa / ef;
                omega += R / Lhs;
                if ( nite > 40 ) {
                    _error("computeDamageParamForCohesiveCrack: algorithm not converging");
                }
            } while ( fabs(R) >= e0 * IDM1_ITERATION_LIMIT );
        } else {
            OOFEM_ERROR1("Unknown softening type for cohesive crack model.");
        }

        if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
            _error("computeDamageParam: damage parameter out of range, snap-back problems");
        }
    }
}
double
IsotropicDamageMaterial1 :: damageFunction(double kappa, GaussPoint *gp)
{
    const double e0 = this->give(e0_ID, gp);
    double ef = 0.;
    if ( softType == ST_Linear || softType == ST_Exponential || softType == ST_SmoothExtended ) {
        ef = this->give(ef_ID, gp);         // ef is the fracturing strain
    }

    switch ( softType ) {
    case ST_Linear:
        if ( kappa <= e0 ) {
            return 0.0;
        } else if ( kappa < ef ) {
            return ( ef / kappa ) * ( kappa - e0 ) / ( ef - e0 );
        } else {
            return 1.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
        }

    case ST_Exponential:
        if ( kappa > e0 ) {
            return 1.0 - ( e0 / kappa ) * exp( -( kappa - e0 ) / ( ef - e0 ) );
        } else {
            return 0.0;
        }

    case ST_Mazars:
        return 1.0 - ( 1.0 - At ) * e0 / kappa - At *exp( -Bt *( kappa - e0 ) );

    case ST_Smooth:
        return 1.0 - exp( -pow(kappa / e0, md) );

    case ST_SmoothExtended:
        // a special damage law used in Grassl and Jirasek (2010)
        if ( kappa <= e1 ) {
            return 1.0 - exp( -pow(kappa / e0, md) );
        } else {
            return 1.0 - s1 *exp( -( kappa - e1 ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) ) / kappa;
        }

    default:
        printf("IsotropicDamageMaterial1::damageFunction ... undefined softening type %d\n", softType);
    }

    return 0.;         // to make the compiler happy
}

double
IsotropicDamageMaterial1 :: complianceFunction(double kappa, GaussPoint *gp)
{
    double om = damageFunction(kappa, gp);
    return om / ( 1. - om );
}

void
IsotropicDamageMaterial1 :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp)
{
    int i, indx = 1;
    double le;
    double E = this->giveLinearElasticMaterial()->give('E', gp);
    FloatArray principalStrains, crackPlaneNormal(3), fullstrain, crackVect(3);
    FloatMatrix principalDir(3, 3);
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    const double e0 = this->give(e0_ID, gp);
    const double ef = this->give(ef_ID, gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);


    if ( softType == ST_Disable_Damage ) {
        return;
    }

    if ( gf != 0. ) { //cohesive crack model
        if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
            wf = gf / E / e0; // wf is the crack opening
        } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi) linear softening law
            wf = 2. * gf / E / e0; // wf is the crack opening
        } else {
            OOFEM_ERROR2("Gf unsupported for softening type softType = %d", softType);
        }
    }

    crossSection->giveFullCharacteristicVector(fullstrain, gp, strainVector);

    if ( ( kappa > e0 ) && ( status->giveDamage() == 0. ) ) {
        this->computePrincipalValDir(principalStrains, principalDir, fullstrain, principal_strain);
        // find index of max positive principal strain
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        // find index with minimal value but non-zero for plane-stress condition - this is the crack direction
        indx = 1;
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) < principalStrains.at(indx) && fabs(principalStrains.at(i))>1.e-10 ) {
                indx = i;
            }
        }
        for ( i = 1; i <= 3; i++ ) {
            crackVect.at(i) = principalDir.at(i, indx);
        }
        status->setCrackVector(crackVect);

        // old approach (default projection method)
        // le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
        // new approach, with choice of method
        le = gp->giveElement()->giveCharacteristicSize(gp, crackPlaneNormal, ecsMethod);
        // remember le in corresponding status
        status->setLe(le);

        // compute and store the crack angle (just for postprocessing)
        double ca = 3.1415926 / 2.;
        if ( crackPlaneNormal.at(1) != 0.0 ) {
            ca = atan( crackPlaneNormal.at(2) / crackPlaneNormal.at(1) );
        }

        status->setCrackAngle(ca);
        
        
        
        if ( this->gf != 0. && e0 >= ( wf / le ) ) { // case for a given fracture energy
            _error6("Fracturing strain %e is lower than the elastic strain e0=%e, possible snap-back. Element number %d, wf %e, le %e", wf / le, e0, gp->giveElement()->giveLabel(), wf, le );
        } else if ( wf == 0. && e0 >= ef ) {
            _error5("Fracturing strain ef=%e is lower than the elastic strain e0=%f, possible snap-back. Increase fracturing strain to %f. Element number %d", ef, e0, e0, gp->giveElement()->giveLabel());
        } else if ( ef == 0. && e0 * le >= wf ) {
            _error5("Crack opening at zero stress wf=%f is lower than the elastic displacement w0=%f, possible snap-back. Increase crack opening wf to %f. Element number %d", wf, e0 * le, e0 * le,gp->giveElement()->giveLabel());
        }
    }
}

double
IsotropicDamageMaterial1 :: give(int aProperty, GaussPoint *gp)
{
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        return answer;
    } else if ( aProperty == e0_ID ) {
        return this->e0;
    } else if ( aProperty == ef_ID ) {
        return this->ef;
    } else if ( aProperty == ek_ID ) {
        return this->ek;
    } else if ( aProperty == gf_ID ) {
        return this->gf;
    } else if ( aProperty == wf_ID ) {
        return this->wf;
    } else if ( aProperty == gft_ID ) {
        return this->gft;
    } else {
        return IsotropicDamageMaterial :: give(aProperty, gp);
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


MaterialStatus *
IsotropicDamageMaterial1 :: CreateStatus(GaussPoint *gp) const
{
    IsotropicDamageMaterial1Status *answer = new IsotropicDamageMaterial1Status(1, IsotropicDamageMaterial1 :: domain, gp);
    return answer;
}

MaterialStatus *
IsotropicDamageMaterial1 :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status;

    status = gp->giveMaterialStatus();
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus(status);
            this->_generateStatusVariables(gp);
        }
    }

    return status;
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
    IsotropicDamageMaterialStatus(n, d, g), RandomMaterialStatusExtensionInterface()
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

Interface *
IsotropicDamageMaterial1Status :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return ( RandomMaterialStatusExtensionInterface * ) this;
    } else {
        return NULL;
    }
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
}     // end namespace oofem
