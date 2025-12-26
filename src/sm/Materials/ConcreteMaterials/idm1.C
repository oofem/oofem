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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "idm1.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "strainvector.h"
#include "stressvector.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"
#include "crosssection.h"
#include "oofemtxtinputrecord.h"

namespace oofem {
REGISTER_Material(IsotropicDamageMaterial1);

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
{
    // deleted by parent, where linearElasticMaterial instance declared
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}


IsotropicDamageMaterial1 :: ~IsotropicDamageMaterial1()
{
    if ( sourceElemSet ) {
        delete sourceElemSet;
    }
}

void
IsotropicDamageMaterial1 :: initializeFrom(InputRecord &ir)
{
    int equivStrainTypeRecord;
    IsotropicDamageMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    checkSnapBack = 1; //check by default
    IR_GIVE_OPTIONAL_FIELD(ir, checkSnapBack, _IFT_IsotropicDamageMaterial1_checkSnapBack);

    // specify the type of formula for equivalent strain
    equivStrainTypeRecord = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, equivStrainTypeRecord, _IFT_IsotropicDamageMaterial1_equivstraintype);
    if ( equivStrainTypeRecord == 0 ) {
        this->equivStrainType = EST_Mazars;
    } else if ( equivStrainTypeRecord == 1 ) {
        this->equivStrainType = EST_Rankine_Smooth;
    } else if ( equivStrainTypeRecord == 2 ) {
        this->equivStrainType = EST_ElasticEnergy;
    } else if ( equivStrainTypeRecord == 3 ) {
        this->equivStrainType = EST_Mises;
        IR_GIVE_FIELD(ir, k, _IFT_IsotropicDamageMaterial1_k);
    } else if ( equivStrainTypeRecord == 4 ) {
        this->equivStrainType = EST_Rankine_Standard;
    } else if ( equivStrainTypeRecord == 5 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStress;
    } else if ( equivStrainTypeRecord == 6 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStrain;
    } else if ( equivStrainTypeRecord == 7 ) {
        this->equivStrainType = EST_Griffith;
        IR_GIVE_OPTIONAL_FIELD(ir, griff_n, _IFT_IsotropicDamageMaterial1_n);
    } else {
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_equivstraintype, "Unknown equivStrainType");
    }

    // specify the type of formula for damage evolution law
    IR_GIVE_OPTIONAL_FIELD(ir, damageLaw, _IFT_IsotropicDamageMaterial1_damageLaw);
    if ( ( damageLaw != 6 ) && ( damageLaw != 7 )  ) {
        IR_GIVE_FIELD(ir, e0, _IFT_IsotropicDamageMaterial1_e0);
    }

    //applies only in this class
    switch ( damageLaw ) {
    case 0:     // exponential softening - default
        if ( ir.hasField(_IFT_IsotropicDamageMaterial1_wf) ) {
            this->softType = ST_Exponential_Cohesive_Crack;
            IR_GIVE_FIELD(ir, wf, _IFT_IsotropicDamageMaterial1_wf);
        } else if ( ir.hasField(_IFT_IsotropicDamageMaterial1_gf) ) {
            this->softType = ST_Exponential_Cohesive_Crack;
            IR_GIVE_FIELD(ir, gf, _IFT_IsotropicDamageMaterial1_gf);
        } else {
            this->softType = ST_Exponential;
            IR_GIVE_FIELD(ir, ef, _IFT_IsotropicDamageMaterial1_ef);
        }

        break;
    case 1:     // linear softening law
        if ( ir.hasField(_IFT_IsotropicDamageMaterial1_wf) ) {
            this->softType = ST_Linear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, wf, _IFT_IsotropicDamageMaterial1_wf);
        } else if ( ir.hasField(_IFT_IsotropicDamageMaterial1_gf) ) {
            this->softType = ST_Linear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, gf, _IFT_IsotropicDamageMaterial1_gf);
        } else {
            this->softType = ST_Linear;
            IR_GIVE_FIELD(ir, ef, _IFT_IsotropicDamageMaterial1_ef);
        }

        break;
    case 2:     // bilinear softening
        gf  = 0.;
        gft = 0.;
        ek  = 0.;
        wf  = 0.;
        wk  = 0.;
        sk  = 0.;
        if ( ir.hasField(_IFT_IsotropicDamageMaterial1_gf) ) {
            this->softType = ST_BiLinear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, gf, _IFT_IsotropicDamageMaterial1_gf);
            // Gft is for the bilinear law, and corresponds to the total energy required to fail the specimen
            IR_GIVE_FIELD(ir, gft, _IFT_IsotropicDamageMaterial1_gft);

            if ( ir.hasField(_IFT_IsotropicDamageMaterial1_ek) ) {
                // ek is for the bilinear law, and corresponds to the strain at the knee point
                IR_GIVE_FIELD(ir, ek, _IFT_IsotropicDamageMaterial1_ek);
            } else {
                double dumWf1, E;

                // wk is for the bilinear law, and corresponds to the crack opening at the knee point
                IR_GIVE_FIELD(ir, wk, _IFT_IsotropicDamageMaterial1_wk);
                IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);

                dumWf1 = 2. * gf / ( e0 * E );
                if ( dumWf1 < wk ) {
                    throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_wk, "Bilinear softening: wk is larger then gf allows");
                }
                sk = ( e0 * E ) * ( 1 - wk / dumWf1 );
                wf = 2. * ( gft - e0 * E * wk / 2.0 ) / sk;

                // for computational purposes
                gf = 0.;
                gft = 0.;
            }
        } else if ( ir.hasField(_IFT_IsotropicDamageMaterial1_wk) ) {
            double E;

            this->softType = ST_BiLinear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, wf, _IFT_IsotropicDamageMaterial1_wf);
            // wk is for the bilinear law, and corresponds to the crack opening at the knee point
            IR_GIVE_FIELD(ir, wk, _IFT_IsotropicDamageMaterial1_wk);
            // sk is for the bilinear law, and corresponds to the stress at the knee point
            IR_GIVE_FIELD(ir, sk, _IFT_IsotropicDamageMaterial1_sk);
            IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);

            if ( wk < 0.0 || wk > wf ) {
                throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_wk, "Bilinear softening: wk must be in interval <0;wf>");
            }
            if ( sk < 0.0 || sk > e0 * E ) {
                throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_sk, "Bilinear softening: sk must be in interval <0;ft>");
            }
        } else if ( ir.hasField(_IFT_IsotropicDamageMaterial1_wkwf) ) {
            double dummy, E;
            this->softType = ST_BiLinear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, wf, _IFT_IsotropicDamageMaterial1_wf);
            // wkwf is for the bilinear law, and corresponds to the ratio of crack opening at the knee point and max crack opening
            IR_GIVE_FIELD(ir, dummy, _IFT_IsotropicDamageMaterial1_wkwf);
            if ( dummy < 0.0 || dummy > 1.0 ) {
                throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_wkwf, "Bilinear softening: wk/wf ratio (wkwf) must be in interval <0;1>");
            } else {
                wk = dummy * wf;
            }
            // sk is for the bilinear law, and corresponds to the ratio of stress at the knee point an tensile strength
            IR_GIVE_FIELD(ir, dummy, _IFT_IsotropicDamageMaterial1_skft);
            if ( dummy < 0.0 || dummy > 1.0 ) {
                throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_skft, "Bilinear softening: sk/ft ratio (skft) must be in interval <0;1>");
            } else {
                IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);
                sk = dummy * e0 * E;
            }
        } else {
            throw ValueInputException(ir, "none", "Bilinear softening for ef not implemented");
        }

        // check if the model is reduced to linear softening
        if ( gf == 0.0 && gft != 0 ) {
            OOFEM_WARNING("Bilinear softening: parameters defined as for Linear_Cohesive_Crack");
            this->softType = ST_Linear_Cohesive_Crack;
            gf = gft;
        } else if ( gft < gf ) {
            throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_gft, "Bilinear softening: gft < gf");
        } else if ( wk == 0.0 && wf != 0 ) {
            OOFEM_WARNING("Bilinear softening: parameters defined as for Linear_Cohesive_Crack");
            this->softType = ST_Linear_Cohesive_Crack;
        } else if ( wf < wk ) {
            throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_wf, "Bilinear softening: wf < wk");
        } else if ( gf == 0 && sk == 0.0 ) {
            OOFEM_WARNING("Bilinear softening: parameters defined as for Linear_Cohesive_Crack");
            this->softType = ST_Linear_Cohesive_Crack;
            wf = wk;
        }

        break;
    case 3:
        this->softType = ST_Hordijk_Cohesive_Crack;
        c1 = 3.;
        IR_GIVE_OPTIONAL_FIELD(ir, c1, _IFT_IsotropicDamageMaterial1_c1);
        c2 = 6.93;
        IR_GIVE_OPTIONAL_FIELD(ir, c2, _IFT_IsotropicDamageMaterial1_c2);
        if ( ir.hasField(_IFT_IsotropicDamageMaterial1_wf) ) {
            IR_GIVE_FIELD(ir, wf, _IFT_IsotropicDamageMaterial1_wf);
        } else if ( ir.hasField(_IFT_IsotropicDamageMaterial1_gf) ) {
            IR_GIVE_FIELD(ir, gf, _IFT_IsotropicDamageMaterial1_gf);
            double E;
            IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);
            double aux = c1 * c1 * c1 / ( c2 * c2 * c2 * c2 );
            aux *= ( 6. - ( c2 * c2 * c2 + 3. * c2 * c2 + 6. * c2 + 6. ) * exp(-c2) );
            aux += ( 1. - exp(-c2) ) / c2;
            aux -= 0.5 * ( 1. + c1 * c1 * c1 ) * exp(-c2);
            wf = gf / ( aux * E * e0 );
        } else {
            throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_wf, "wf or gf must be specified for Hordijk softening law");
        }
        break;
    case 4:
        this->softType = ST_Mazars;
        IR_GIVE_FIELD(ir, At, _IFT_IsotropicDamageMaterial1_At);
        IR_GIVE_FIELD(ir, Bt, _IFT_IsotropicDamageMaterial1_Bt);
        break;
    case 5:
        this->softType = ST_Smooth;
        md = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, md, _IFT_IsotropicDamageMaterial1_md);
        break;
    case 6:
        this->softType = ST_Disable_Damage;
        break;
    case 7:
        this->softType = ST_SmoothExtended;
        double E; // auxiliary input variables
        IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);
        IR_GIVE_FIELD(ir, this->ep, _IFT_IsotropicDamageMaterial1_ep);
        IR_GIVE_FIELD(ir, this->ft, _IFT_IsotropicDamageMaterial1_ft);
        this->md = 1. / log(E * this->ep / this->ft);
        this->e0 = ep * pow(md, 1. / md);
        IR_GIVE_FIELD(ir, e1, _IFT_IsotropicDamageMaterial1_e1);
        this->s1 = e1 * exp( -pow(e1 / e0, md) );
        this->ef = -e1 / ( 1. - md * pow(e1 / e0, md) );
        IR_GIVE_FIELD(ir, e2, _IFT_IsotropicDamageMaterial1_e2);
        IR_GIVE_FIELD(ir, nd, _IFT_IsotropicDamageMaterial1_nd);
        break;
    case 8:     // power-exponential softening
        this->softType = ST_PowerExponential;
        IR_GIVE_FIELD(ir, ef, _IFT_IsotropicDamageMaterial1_ef);
        IR_GIVE_FIELD(ir, md, _IFT_IsotropicDamageMaterial1_md);
        break;

    case 9:     // double exponential softening
        this->softType = ST_DoubleExponential;
        IR_GIVE_FIELD(ir, ef, _IFT_IsotropicDamageMaterial1_ef);
        IR_GIVE_OPTIONAL_FIELD(ir, c2, _IFT_IsotropicDamageMaterial1_c2);
        IR_GIVE_FIELD(ir, e2, _IFT_IsotropicDamageMaterial1_e2);
        break;
    case 10:     // modified power-exponential softening
        this->softType = ST_ModPowerExponential;
        IR_GIVE_FIELD(ir, ef, _IFT_IsotropicDamageMaterial1_ef);
        IR_GIVE_FIELD(ir, md, _IFT_IsotropicDamageMaterial1_md);
        break;
    case 11:     // trilinear softening
        this->softType = ST_Trilinear_Cohesive_Crack;
        IR_GIVE_FIELD(ir, w_k, _IFT_IsotropicDamageMaterial1_w_k);
        IR_GIVE_FIELD(ir, w_r, _IFT_IsotropicDamageMaterial1_w_r);
        IR_GIVE_FIELD(ir, w_f, _IFT_IsotropicDamageMaterial1_w_f);
        IR_GIVE_FIELD(ir, f_k, _IFT_IsotropicDamageMaterial1_f_k);
        IR_GIVE_FIELD(ir, f_r, _IFT_IsotropicDamageMaterial1_f_r);
        break;

    default:
        throw ValueInputException(ir, _IFT_IsotropicDamageMaterial1_damageLaw, "Unknown value");
    }

    if ( ( softType == ST_Exponential_Cohesive_Crack ) || ( softType == ST_Linear_Cohesive_Crack ) || ( softType == ST_BiLinear_Cohesive_Crack ) || ( softType == ST_Trilinear_Cohesive_Crack )) {
        int ecsm = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, ecsm, _IFT_IsotropicDamageMaterial1_ecsm);
        switch ( ecsm ) {
        case 1: ecsMethod = ECSM_SquareRootOfArea;
            break;
        case 2: ecsMethod = ECSM_ProjectionCentered;
            break;
        case 3: ecsMethod = ECSM_Oliver1;
            break;
        case 4: ecsMethod = ECSM_Oliver1modified;
            break;
        default: ecsMethod = ECSM_Projection;
        }
    }

    if ( permStrain == 4 ) {
        IR_GIVE_FIELD(ir, ps_alpha, _IFT_IsotropicDamageMaterial1_alphaps);
        //IR_GIVE_FIELD(ir, ps_H, _IFT_IsotropicDamageMaterial1_h);
    }

    this->mapper.initializeFrom(ir);
}


void
IsotropicDamageMaterial1 :: giveInputRecord(DynamicInputRecord &input)
{
    // Done first so that record name is overwritten afterwards by femcomponent:
    this->linearElasticMaterial->giveInputRecord(input);
    this->mapper.giveInputRecord(input);
    IsotropicDamageMaterial :: giveInputRecord(input);
    RandomMaterialExtensionInterface :: giveInputRecord(input);

    input.setField(this->equivStrainType, _IFT_IsotropicDamageMaterial1_equivstraintype);
    input.setField(this->damageLaw, _IFT_IsotropicDamageMaterial1_damageLaw);
    if ( ( damageLaw != 6 ) && ( damageLaw != 7 ) ) {
        input.setField(this->e0, _IFT_IsotropicDamageMaterial1_e0);
    }
    switch ( damageLaw ) {
    case 0:     // exponential softening - default
        if ( this->softType == ST_Exponential_Cohesive_Crack ) {
            input.setField(this->wf, _IFT_IsotropicDamageMaterial1_wf);
            input.setField(this->gf, _IFT_IsotropicDamageMaterial1_gf);
        } else if ( this->softType == ST_Exponential ) {
            input.setField(this->ef, _IFT_IsotropicDamageMaterial1_ef);
        }
        break;
    case 1:
        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            if ( this->wf != 0.0 ) {
                input.setField(this->wf, _IFT_IsotropicDamageMaterial1_wf);
            } else if ( this->gf != 0.0 ) {
                input.setField(this->gf, _IFT_IsotropicDamageMaterial1_gf);
            } else {
                input.setField(this->ef, _IFT_IsotropicDamageMaterial1_ef);
            }
        }
        break;
    case 2:
        if ( this->softType == ST_BiLinear_Cohesive_Crack ) {
            input.setField(this->gf, _IFT_IsotropicDamageMaterial1_gf);
            input.setField(this->gft, _IFT_IsotropicDamageMaterial1_gft);
            if ( this->ek != 0.0 ) {
                input.setField(this->ek, _IFT_IsotropicDamageMaterial1_ek);
            } else {
                input.setField(this->wk, _IFT_IsotropicDamageMaterial1_wk);
            }
        } else if ( this->wk != 0.0 ) {
            input.setField(this->wf, _IFT_IsotropicDamageMaterial1_wf);
            input.setField(this->wk, _IFT_IsotropicDamageMaterial1_wk);
            input.setField(this->sk, _IFT_IsotropicDamageMaterial1_sk);
        }
        break;
    case 3:
        input.setField(this->wf, _IFT_IsotropicDamageMaterial1_wf);
        break;
    case 4:
        input.setField(this->At, _IFT_IsotropicDamageMaterial1_At);
        input.setField(this->Bt, _IFT_IsotropicDamageMaterial1_Bt);
        break;
    case 5:
        input.setField(this->md, _IFT_IsotropicDamageMaterial1_md);
        break;
    case 7:
        input.setField(this->ep, _IFT_IsotropicDamageMaterial1_ep);
        input.setField(this->ft, _IFT_IsotropicDamageMaterial1_ft);
        input.setField(this->e1, _IFT_IsotropicDamageMaterial1_e1);
        input.setField(this->e2, _IFT_IsotropicDamageMaterial1_e2);
        input.setField(this->nd, _IFT_IsotropicDamageMaterial1_nd);
        break;
    case 8:
        input.setField(this->ef, _IFT_IsotropicDamageMaterial1_ef);
        input.setField(this->md, _IFT_IsotropicDamageMaterial1_md);
        break;
    case 11: // Trilinear softening
        input.setField(this->w_k, _IFT_IsotropicDamageMaterial1_w_k);
        input.setField(this->w_r, _IFT_IsotropicDamageMaterial1_w_r);
        input.setField(this->w_f, _IFT_IsotropicDamageMaterial1_w_f);
        input.setField(this->f_k, _IFT_IsotropicDamageMaterial1_f_k);
        input.setField(this->f_r, _IFT_IsotropicDamageMaterial1_f_r);
        break;
    }
    if ( softType == ST_Exponential_Cohesive_Crack || softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack || softType == ST_Trilinear_Cohesive_Crack ) {
        input.setField(this->ecsMethod, _IFT_IsotropicDamageMaterial1_ecsm);
    }
}


double
IsotropicDamageMaterial1 :: computeEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto lmat = this->linearElasticMaterial;

    if ( strain.isEmpty() ) {
        return 0.;
    }

    FloatArray fullStrain;
    StructuralMaterial :: giveFullSymVectorForm( fullStrain, strain, gp->giveMaterialMode() );
    // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
    if ( gp->giveMaterialMode() == _PlaneStress ) {
        double nu = lmat->give(NYxz, gp);
        fullStrain.at(3) = -nu * ( fullStrain.at(1) + fullStrain.at(2) ) / ( 1. - nu );
    } else if ( gp->giveMaterialMode() == _1dMat ) {
        double nu = lmat->give(NYxz, gp);
        fullStrain.at(2) = -nu *fullStrain.at(1);
        fullStrain.at(3) = -nu *fullStrain.at(1);
    }

    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains;

        this->computePrincipalValues(principalStrains, fullStrain, principal_strain);

        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        return sqrt(posNorm);
    } else if ( ( this->equivStrainType == EST_Rankine_Smooth ) || ( this->equivStrainType == EST_Rankine_Standard ) ) {
        // EST_Rankine equiv strain measure
        double sum = 0.;
        FloatArray stress, fullStress, principalStress;
        FloatMatrix de;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        stress.beProductOf(de, strain);
        StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                if ( this->equivStrainType == EST_Rankine_Smooth ) {
                    sum += principalStress.at(i) * principalStress.at(i);
                } else if ( sum < principalStress.at(i) ) {
                    sum = principalStress.at(i);
                }
            } else if ( sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
            }
        }

        if ( this->equivStrainType == EST_Rankine_Smooth ) {
            sum = sqrt(sum);
        }

        return sum / lmat->give('E', gp);
    } else if ( ( this->equivStrainType == EST_ElasticEnergy ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStrain ) ) {
        // equivalent strain expressions based on elastic energy
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        if ( this->equivStrainType == EST_ElasticEnergy ) {
            // standard elastic energy
            stress.beProductOf(de, strain);
            sum = strain.dotProduct(stress);
        } else if ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) {
            // elastic energy corresponding to positive part of stress
            FloatArray fullStress, principalStress;
            StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
            this->computePrincipalValues(principalStress, fullStress, principal_stress);
            // TO BE FINISHED
            sum = 0.;
            OOFEM_ERROR("Elastic energy corresponding to positive part of stress not finished");
        } else {
            // elastic energy corresponding to positive part of strain
            // TO BE DONE
            sum = 0.;
            OOFEM_ERROR("Elastic energy corresponding to positive part of strain not finished");
        }

        return sqrt( sum / lmat->give('E', gp) );
    } else if ( this->equivStrainType == EST_Mises ) {
        double nu = lmat->give(NYxz, NULL);
        FloatArray principalStrains;

        this->computePrincipalValues(principalStrains, fullStrain, principal_strain);
        double I1e, J2e;
        this->computeStrainInvariants(principalStrains, I1e, J2e);
        double a, b, c;
        a = ( k - 1 ) * I1e / ( 2 * k * ( 1 - 2 * nu ) );
        b = ( k - 1 ) * ( k - 1 ) * I1e * I1e / ( ( 1 - 2 * nu ) * ( 1 - 2 * nu ) );
        c = 12 * k * J2e / ( ( 1 + nu ) * ( 1 + nu ) );
        return a + 1 / ( 2 * k ) * sqrt(b + c);
    } else if ( this->equivStrainType == EST_Griffith ) {
        double kappa1 = 0.0, kappa2 = 0.0;
        FloatArray stress, fullStress, principalStress;
        FloatMatrix de;
        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        stress.beProductOf(de, strain);
        StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        double sum = 0., maxStress = 0.;
        //Compute equivalent strain for Rankine's criterion
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 && sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
            }
        }
        kappa1 = sum / lmat->give('E', gp);
        //Compute equivalent strain for Griffith's criterion
        // Check zero division first
        maxStress = max( fabs( principalStress.at(1) ), fabs( principalStress.at(3) ) );
        if ( maxStress == 0. || fabs( principalStress.at(3) ) < 1.e-6 * maxStress || fabs( principalStress.at(1) + principalStress.at(3) ) < 1.e-6 * maxStress ) {
            //Skip evaluation
        } else if ( principalStress.at(1) / principalStress.at(3) >= -0.33333 ) {
            kappa2 = -( principalStress.at(1) - principalStress.at(3) ) * ( principalStress.at(1) - principalStress.at(3) ) / this->griff_n / ( principalStress.at(1) + principalStress.at(3) ) / lmat->give('E', gp);
        }
        return max(max(kappa1, 0.0), kappa2);
    } else {
        OOFEM_ERROR("unknown EquivStrainType");
        return 0.;
    }
}

//Computes derivative of the equivalent strain with regards to strain, used in tangent formulation
void
IsotropicDamageMaterial1 :: computeEta(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const
{
    LinearElasticMaterial *lmat = this->linearElasticMaterial;

    if ( strain.isEmpty() ) {
        answer.zero();
        return;
    }

    if ( this->equivStrainType == EST_Mazars ) {
        int dim = 0;
        double posNorm = 0.0;
        double nu = lmat->give(NYxz, gp);
        FloatArray principalStrains;
        FloatMatrix N;

        if ( strain.giveSize() == 6 || gp->giveMaterialMode() == _3dMat ) {
            dim = 3;
            this->computePrincipalValDir(principalStrains, N, strain, principal_strain);
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            dim = 1;
            StrainVector fullStrain(strain, _1dMat);
            fullStrain.computePrincipalValDir(principalStrains, N);
            principalStrains.resizeWithValues(3);
            principalStrains.at(2) = -nu *principalStrains.at(1);
            principalStrains.at(3) = -nu *principalStrains.at(1);
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            dim = 2;
            StrainVector fullStrain(strain, _PlaneStress);
            fullStrain.computePrincipalValDir(principalStrains, N);
            principalStrains.resizeWithValues(3);
            principalStrains.at(3) = -nu * ( principalStrains.at(1) + principalStrains.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            dim = 3;
            StrainVector fullStrain(strain, _PlaneStrain);
            fullStrain.computePrincipalValDir(principalStrains, N);
        } else {
            dim = 0;
            OOFEM_ERROR("Unknown material mode.");
        }
        FloatArray n(dim);
        FloatMatrix Eta(dim, dim);
        Eta.zero();

        for ( int i = 1; i <= 3; i++ ) {
            if ( i <= dim ) {
                if ( principalStrains.at(i) > 0.0 ) {
                    n.beColumnOf(N, i);
                    Eta.plusDyadSymmUpper( n, principalStrains.at(i) );
                }
            }

            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        Eta.symmetrized();

        double kappa = sqrt(posNorm);

        int numberOfEl = ( dim * ( dim - 1 ) / 2 + dim );
        answer.resize(numberOfEl);

        if ( kappa != 0 ) {
            Eta.times(1. / kappa);
        } else {
            answer.zero();
            return;
        }

        if ( strain.giveSize() == 6 || gp->giveMaterialMode() == _3dMat ) {
            answer.at(1) = Eta.at(1, 1);
            answer.at(2) = Eta.at(2, 2);
            answer.at(3) = Eta.at(3, 3);
            answer.at(4) = Eta.at(2, 3);
            answer.at(5) = Eta.at(1, 3);
            answer.at(6) = Eta.at(1, 2);
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            answer.at(1) = Eta.at(1, 1);
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            answer.at(1) = Eta.at(1, 1);
            answer.at(2) = Eta.at(2, 2);
            answer.at(3) = Eta.at(1, 2);
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            answer.resize(4);
            answer.at(1) = Eta.at(1, 1);
            answer.at(2) = Eta.at(2, 2);
            answer.at(3) = Eta.at(3, 3);
            answer.at(4) = Eta.at(1, 2);
        }
    } else if ( ( this->equivStrainType == EST_Rankine_Smooth ) || ( this->equivStrainType == EST_Rankine_Standard ) ) {
        int index = 0, dim = 0;
        double sum = 0.;
        FloatArray stress, principalStress, eta;
        FloatMatrix de, N, Eta;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        stress.beProductOf(de, strain);

        if ( strain.giveSize() == 6 || gp->giveMaterialMode() == _3dMat ) {
            this->computePrincipalValDir(principalStress, N, strain, principal_stress);
            dim = 3;
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            StressVector fullStress(stress, _1dMat);
            fullStress.computePrincipalValDir(principalStress, N);
            principalStress.resizeWithValues(3);
            dim = 1;
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            StressVector fullStress(stress, _PlaneStress);
            fullStress.computePrincipalValDir(principalStress, N);
            principalStress.resizeWithValues(3);
            dim = 2;
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            StressVector fullStress(stress, _PlaneStrain);
            fullStress.computePrincipalValDir(principalStress, N);
            dim = 3;
        } else {
            dim = 0;
            OOFEM_ERROR("Unknown material mode.");
        }

        FloatArray n(dim);
        n.zero();
        Eta.resize(dim, dim);
        Eta.zero();
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                if ( this->equivStrainType == EST_Rankine_Smooth ) {
                    sum += principalStress.at(i) * principalStress.at(i);
                    if ( i <= dim ) {
                        for ( int j = 1; j <= dim; j++ ) {
                            n.at(j) = N.at(j, i);
                        }
                    }

                    Eta.plusDyadSymmUpper( n, principalStress.at(i) );
                } else if ( sum < principalStress.at(i) ) {
                    sum = principalStress.at(i);
                    index = i;
                }
            } else if ( sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
                index = i;
            }
        }

        Eta.symmetrized();

        int numberOfEl = ( dim * ( dim - 1 ) / 2 + dim );
        eta.resize(numberOfEl);

        if ( this->equivStrainType == EST_Rankine_Smooth ) {
            sum = sqrt(sum);
            double kappa =  sum / lmat->give('E', gp);

            if ( kappa != 0 ) {
                Eta.times(1. / kappa);
            } else {
                answer.zero();
                return;
            }

            Eta.times( 1. / kappa / lmat->give('E', gp) / lmat->give('E', gp) );
        } else if ( this->equivStrainType == EST_Rankine_Standard ) {
            for ( int i = 1; i <= dim; i++ ) {
                n.at(i) = N.at(i, index);
            }

            Eta.beDyadicProductOf(n, n);
            Eta.times( 1. / lmat->give('E', gp) );
        }

        if ( strain.giveSize() == 6 || gp->giveMaterialMode() == _3dMat ) {
            eta.at(1) = Eta.at(1, 1);
            eta.at(2) = Eta.at(2, 2);
            eta.at(3) = Eta.at(3, 3);
            eta.at(4) = Eta.at(2, 3);
            eta.at(5) = Eta.at(1, 3);
            eta.at(6) = Eta.at(1, 2);
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            eta.at(1) = Eta.at(1, 1);
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            eta.at(1) = Eta.at(1, 1);
            eta.at(2) = Eta.at(2, 2);
            eta.at(3) = Eta.at(1, 2);
        }  else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            eta.resize(4);
            eta.at(1) = Eta.at(1, 1);
            eta.at(2) = Eta.at(2, 2);
            eta.at(3) = Eta.at(3, 3);
            eta.at(4) = 2. * Eta.at(1, 2);
        }

        answer.beProductOf(de, eta);
    } else if ( this->equivStrainType == EST_ElasticEnergy ) {
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        // standard elastic energy
        stress.beProductOf(de, strain);
        sum = strain.dotProduct(stress);

        double kappa = sqrt( sum / lmat->give('E', gp) );
        answer = stress;
        if ( kappa != 0 ) {
            answer.times(1. / lmat->give('E', gp) / kappa);
        } else {
            answer.times(0.);
        }
    } else {
        OOFEM_ERROR("unknown EquivStrainType");
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

double
IsotropicDamageMaterial1 :: computeDamageParam(double kappa, const FloatArray &strain, GaussPoint *gp) const
{
    if ( this->softType == ST_Disable_Damage ) { //dummy material with no damage
        return 0.;
    } else if ( isCrackBandApproachUsed() ) { // adjustment of softening law according to the element size, given crack opening or fracture energy
        return computeDamageParamForCohesiveCrack(kappa, gp);
    } else { // no adjustment according to element size, given fracturing strain
        return damageFunction(kappa, gp);
    }
}

double
IsotropicDamageMaterial1 :: computeDamageParamForCohesiveCrack(double kappa, GaussPoint *gp) const
{
    const double e0 = this->give(e0_ID, gp);  // e0 is the strain at the peak stress
    const double E = this->linearElasticMaterial->give('E', gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);     // wf is the crack opening
    double omega = 0.0;

    if ( kappa > e0 ) {
        if ( this->gf != 0. ) { //cohesive crack model
            if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
                wf = this->gf / E / e0; // wf is the crack opening
            } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi)linear softening law
                wf = 2. * gf / E / e0; // wf is the crack opening
            } else {
                OOFEM_ERROR("Gf unsupported for softening type softType = %d", softType);
            }
        } else if ( softType == ST_BiLinear_Cohesive_Crack ) {
            wf = this->wk / ( e0 * E - this->sk ) * ( e0 * E );
        } else if ( softType == ST_Trilinear_Cohesive_Crack ) {
            wf = this->w_f;
        }


        auto status = static_cast< IsotropicDamageMaterial1Status * >( this->giveStatus(gp) );
        double Le = status->giveLe();
        double ef = wf / Le;    //ef is the fracturing strain /// FIXME CHANGES BEHAVIOR!
        if ( ef < e0 ) { //check that no snapback occurs
            double minGf = 0.;
            OOFEM_WARNING("ef %e < e0 %e, this leads to material snapback in element %d, characteristic length %f", ef, e0, gp->giveElement()->giveGlobalNumber(), Le);
            if ( gf != 0. ) { //cohesive crack
                if ( softType == ST_Exponential_Cohesive_Crack ) { //exponential softening
                    minGf = E * e0 * e0 * Le;
                } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { //(bi)linear softening law
                    minGf = E * e0 * e0 * Le / 2.;
                } else {
                    OOFEM_WARNING("Gf unsupported for softening type softType = %d", softType);
                }

                if ( checkSnapBack ) {
                    OOFEM_ERROR("Material number %d, decrease e0, or increase Gf from %f to Gf=%f", this->giveNumber(), gf, minGf);
                }
            }

            if ( checkSnapBack ) { //given fracturing strain
                OOFEM_ERROR("Material number %d, increase ef %f to minimum e0 %f", this->giveNumber(), ef, e0);
            }
        }

        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            if ( kappa < ef ) {
                omega = ( ef / kappa ) * ( kappa - e0 ) / ( ef - e0 );
            } else {
                omega = 1.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
            }
        } else if (  this->softType == ST_BiLinear_Cohesive_Crack ) {
            double gft = this->give(gft_ID, gp);
            double ef, sigmak, epsf, ek;
            if ( gft > 0.0 ) {
                ek = this->give(ek_ID, gp);
                ef = 2 * gf / E / e0 / Le; //the first part corresponds to linear softening
                sigmak = E * e0 * ( ef - ek ) / ( ef - e0 );
                epsf = 2 * ( gft - gf ) / sigmak / Le + ef;

                if ( gft < gf ) {
                    OOFEM_ERROR("The total fracture energy gft %f must be greater than the initial fracture energy gf %f", gft, gf);
                }
            } else {
                ek     = this->wk / Le + ( this->sk ) / E;
                ef     = ( this->wk / ( e0 * E - this->sk ) * ( e0 * E ) ) / Le;
                sigmak = this->sk;
                epsf   = this->wf / Le;
            }
            if ( ( ek > ef ) || ( ek < e0 ) ) {
                OOFEM_WARNING("ek %f is not between e0 %f and ef %f", ek, e0, ef);
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
                R = ( 1. - omega ) * kappa - e0 *exp(-help); //residuum
                Lhs = kappa - e0 *exp(-help) * kappa / ef; //- dR / (d omega)
                omega += R / Lhs;
                if ( nite > 40 ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(R) >= e0 * IDM1_ITERATION_LIMIT );
        } else if ( this->softType == ST_Trilinear_Cohesive_Crack ) {
			double eps_k = this->w_k/Le + this->f_k/E;
			double eps_r = this->w_r/Le + this->f_r/E;
			double eps_f = this->w_f/Le ;
			double f_t = E*e0;
			if ( kappa > e0 && kappa <= eps_k ) {
				double slope=((f_k - f_t)/w_k);
				omega = E/(E+slope*Le) - (f_t)/(kappa*(E+slope*Le));
			} else if ( kappa > eps_k && kappa <= eps_r ) {
				double slope=((f_r - f_k)/(w_r-w_k));
				omega = E/(E+slope*Le) +(1.0/kappa)*(w_k*slope-f_k)/(Le*slope+E);
			} else if ( kappa > eps_r && kappa <= eps_f ) {
				double slope=(- f_r/(w_f-w_r));
				omega = E/(E+slope*Le) +(1.0/kappa)*(w_r*slope-f_r)/(Le*slope+E);
			} else {
				omega = 1.0;
			}
        } else {
            OOFEM_ERROR("Unknown softening type for cohesive crack model.");
        }

        if ( omega > 1.0 ) {
            OOFEM_WARNING("damage parameter is %f, which is greater than 1, snap-back problems", omega);
            omega = maxOmega;
            if ( checkSnapBack ) {
                OOFEM_ERROR("x");
            }
        }

        if ( omega < 0.0 ) {
            OOFEM_WARNING("damage parameter is %f, which is smaller than 0, snap-back problems", omega);
            omega = 0.0;
            if ( checkSnapBack ) {
                OOFEM_ERROR("x");
            }
        }
    }
    return omega;
}


double
IsotropicDamageMaterial1 :: damageFunction(double kappa, GaussPoint *gp) const
{
    const double e0 = this->give(e0_ID, gp);
    double ef = 0.;
    if ( softType == ST_Linear || softType == ST_Exponential || softType == ST_SmoothExtended || softType == ST_PowerExponential || softType == ST_ModPowerExponential || softType == ST_DoubleExponential ) {
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
            if ( permStrain == 4 ) {
                return ( kappa - e0 * exp( -( kappa - e0 ) / ( ef - e0 ) ) ) / ( kappa + ps_alpha );
            } else {
                return 1.0 - ( e0 / kappa ) * exp( -( kappa - e0 ) / ( ef - e0 ) );
            }
        }
        return 0.0;

    case ST_DoubleExponential:
        if ( kappa > e0 ) {
            return 1.0 - ( 1. - c2 ) * ( e0 / kappa ) * exp( -( kappa - e0 ) / ( ef - e0 ) ) - c2 * ( e0 / kappa ) * exp( -( kappa - e0 ) / ( e2 - e0 ) );
        } else {
            return 0.0;
        }

    case ST_PowerExponential:
        if ( kappa > e0 ) {
            return 1.0 - ( e0 / kappa ) * exp( -pow( ( kappa - e0 ) / ( ef - e0 ), md ) );
        } else {
            return 0.0;
        }

    case ST_ModPowerExponential:
        if ( kappa > e0 ) {
            return 1.0 - ( e0 / kappa ) * exp( -( pow(kappa, md) - pow(e0, md) ) / ( pow(ef, md) - pow(e0, md) ) );
        } else {
            return 0.0;
        }

    case ST_Mazars:
        return 1.0 - ( 1.0 - At ) * e0 / kappa - At *exp( -Bt * ( kappa - e0 ) );

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
        OOFEM_WARNING(":damageFunction ... undefined softening type %d\n", softType);
    }

    return 0.;         // to make the compiler happy
}

double
IsotropicDamageMaterial1 :: damageFunctionPrime(double kappa, GaussPoint *gp) const
{
    const double e0 = this->give(e0_ID, gp);
    double ef = 0.;
    const double E = this->linearElasticMaterial->give('E', gp);
    IsotropicDamageMaterial1Status *status = static_cast< IsotropicDamageMaterial1Status * >( this->giveStatus(gp) );
    const double Le = status->giveLe();

    if ( softType == ST_Linear || softType == ST_Exponential || softType == ST_SmoothExtended ) {
        ef = this->give(ef_ID, gp);         // ef is the fracturing strain
    }

    switch ( softType ) {
    case ST_Linear:
    {
        if ( kappa <= e0 ) {
            return 0.0;
        } else if ( kappa < ef ) {
            return ( ef * e0 ) / ( ef - e0 ) / ( kappa * kappa );
        } else {
            return 0.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
        }
    } break;
    case ST_Exponential:
    {
        if ( kappa > e0 ) {
            //	  return ( e0 / ( kappa * kappa ) ) * exp( -( kappa - e0 ) / ( ef - e0 )  + e0 / ( kappa * ( ef - e0 ) ) ) * exp( -( kappa - e0 ) / ( ef - e0 ) );
            return ( e0 / ( ef - e0 ) / kappa +  e0 / ( kappa * kappa ) )  * exp( -( kappa - e0 ) / ( ef - e0 ) );
        } else {
            return 0.0;
        }
    } break;
    case ST_Mazars:
    {
        return ( 1.0 - At ) * e0 / kappa / kappa - At *Bt *exp( -Bt * ( kappa - e0 ) );
    } break;
    case ST_Smooth:
    {
        return md / e0 *pow(kappa, md - 1.) * exp( -pow(kappa / e0, md) );
    } break;
    case ST_SmoothExtended:
    {
        // a special damage law used in Grassl and Jirasek (2010)
        if ( kappa <= 0 ) {
            return 0;
        } else if ( kappa <= e1 ) {
            return exp( -pow(kappa / e0, md) ) * md / pow(e0, md) * pow(kappa, md - 1.);
        } else {
            double a = ( ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) -  ef * nd * pow( ( kappa - e1 ) / e2, nd ) ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) );
            double answer =   s1 * exp( -( kappa - e1 ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) ) / kappa / kappa +  s1 *exp( -( kappa - e1 ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) ) / kappa * a;
            return answer;
        }
    } break;
    case ST_Linear_Cohesive_Crack:
    {
        double wf = 2. * gf / E / e0; // wf is the crack opening
        if ( kappa <= e0 ) {
            return 0.0;
        } else if ( kappa < wf / Le ) {
            return ( e0 / ( kappa * kappa ) / ( 1. - Le * e0 / wf ) );
        } else {
            return 0.0;
        }
    } break;
    case ST_Exponential_Cohesive_Crack:
    {
        if ( kappa > e0 ) {
            ef = gf / E / e0 / Le;
            double omega = status->giveTempDamage();
            double help = exp(omega * kappa / ef);
            double ret = -( ( omega * ef - ef ) * help - omega * e0 ) / ( ef * kappa * help - e0 * kappa );
            if ( std :: isnan(ret) ) {
                return 0.;
            }
            return ret;
        } else {
            return 0.0;
        }
    } break;
    case ST_Trilinear_Cohesive_Crack:
    {
    	double Le = status->giveLe();
    	const double E = this->linearElasticMaterial->give('E', gp);
    	double eps_k = this->w_k/Le + this->f_k/E;
    	double eps_r = this->w_r/Le + this->f_r/E;
    	double eps_f = this->w_f/Le ;
    	double f_t = E*e0;
        if ( kappa <= e0 ) {
            return 0.0;
        } else if ( kappa > e0 && kappa <= eps_k ) {
        	double slope=((f_k - f_t)/w_k);
        	return (f_t)/(E+slope*Le)*(1.0/(kappa*kappa));
        } else if ( kappa > eps_k && kappa <= eps_r ) {
        	double slope=((f_r - f_k)/(w_r-w_k));
        	return -(w_k*slope-f_k)/(Le*slope+E)*(1.0/(kappa*kappa));
        } else if ( kappa > eps_r && kappa <= eps_f ) {
        	double slope=(- f_r/(w_f-w_r));
        	return -(w_r*slope-f_r)/(Le*slope+E)*(1.0/(kappa*kappa));
        } else {
            return 0.0;
        }
    } break;
    default:
        OOFEM_ERROR("undefined softening type %d\n", softType);
    }

    return 0.;         // to make the compiler happy
}

double
IsotropicDamageMaterial1 :: complianceFunction(double kappa, GaussPoint *gp) const
{
    double om = damageFunction(kappa, gp);
    return om / ( 1. - om );
}

double
IsotropicDamageMaterial1 :: evaluatePermanentStrain(double kappa, double omega) const
{
    switch ( permStrain ) {
    case 1:
        return 0.;

        break;
    case 2:
        return 0.;

        break;
    case 3:
        return 0.;

        break;
    case 4:
        return ps_alpha * omega / ( 1. - omega );

        break;
    default:
        return 0.;
    }
    ;
}


void
IsotropicDamageMaterial1 :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp) const
{
    auto status = static_cast< IsotropicDamageMaterial1Status * >( this->giveStatus(gp) );
    int indx = 1;
    double le = 0.;
    double E = this->linearElasticMaterial->give('E', gp);
    FloatArray principalStrains, crackPlaneNormal, fullStrain, crackVect;
    FloatMatrix principalDir;

    const double e0 = this->give(e0_ID, gp);
    const double ef = this->give(ef_ID, gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);

    if ( softType == ST_Disable_Damage ) {
        return;
    }
    if ( softType == ST_Trilinear_Cohesive_Crack ) {
                wf = this->w_f;
            }

    if ( gf != 0. ) { //cohesive crack model
        if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
            wf = gf / E / e0; // wf is the crack opening
        } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi) linear softening law
            wf = 2. * gf / E / e0; // wf is the crack opening
        } else {
            OOFEM_ERROR("Gf unsupported for softening type softType = %d", softType);
        }
    } else if ( softType == ST_BiLinear_Cohesive_Crack ) {
        wf = this->wk / ( e0 * E - this->sk ) * ( e0 * E );
    }

    StructuralMaterial :: giveFullSymVectorForm( fullStrain, strainVector, gp->giveMaterialMode() );


    if ( ( kappa > e0 ) && ( ( status->giveDamage() == 0. ) || ( status->giveLe() == 0.0 ) ) ) {
        // zero Le can happen after adaptive update; need to recompute Le
        this->computePrincipalValDir(principalStrains, principalDir, fullStrain, principal_strain);
        // find index of max positive principal strain
        for ( int i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        crackPlaneNormal.beColumnOf(principalDir, indx);

        // find index with minimal value but non-zero for plane-stress condition - this is the crack direction
        indx = 1;
        for ( int i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) < principalStrains.at(indx) && fabs( principalStrains.at(i) ) > 1.e-10 ) {
                indx = i;
            }
        }

        // Use orientation of the worst inclusion for Griffith criterion in compression.
        if ( this->equivStrainType == EST_Griffith ) {
            FloatArray stress, fullStress, principalStress, crackV(3); // , crackPlaneN(3);
            FloatMatrix de;
            LinearElasticMaterial *lmat = this->linearElasticMaterial;
            lmat->giveStiffnessMatrix( de, SecantStiffness, gp, domain->giveEngngModel()->giveCurrentStep() );
            stress.beProductOf(de, strainVector);
            StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
            this->computePrincipalValDir(principalStress, principalDir, fullStress, principal_stress);
            if (  principalStress.at(1) <= 1.e-10 && principalStress.at(2) <= 1.e-10 && principalStress.at(3) <= 1.e-10 ) {
                int indexMax = principalStress.giveIndexMaxElem();
                int indexMin = principalStress.giveIndexMinElem();
                int indexMid = 0;
                if ( indexMin + indexMax == 3 ) {
                    indexMid = 3;
                } else if ( indexMin + indexMax == 4 ) {
                    indexMid = 2;
                } else if ( indexMin + indexMax == 5 ) {
                    indexMid = 1;
                }

                //inclination from maximum compressive stress (plane sig_1 and sig_3)
                double twoPsi = ( principalStress.at(indexMin) - principalStress.at(indexMax) ) / 2. / ( principalStress.at(indexMax) + principalStress.at(indexMin) );
                double psi = acos(twoPsi) / 2.;
                for ( int i = 1; i <= 3; i++ ) {
                    crackV.at(i) = principalDir.at(i, indexMin);
                    // crackPlaneN = principalDir.at(i, indexMax);
                }

                //rotate around indexMid axis
                //see http://en.wikipedia.org/wiki/Rotation_matrix and Rodrigues' rotation formula
                FloatMatrix ux(3, 3), dyadU(3, 3), unitMtrx(3, 3), rotMtrx(3, 3);
                FloatArray u(3);
                ux.zero();
                ux.at(1, 2) = -principalDir.at(3, indexMid);
                ux.at(1, 3) = principalDir.at(2, indexMid);
                ux.at(2, 1) = principalDir.at(3, indexMid);
                ux.at(2, 3) = -principalDir.at(1, indexMid);
                ux.at(3, 1) = -principalDir.at(2, indexMid);
                ux.at(3, 2) = principalDir.at(1, indexMid);
                u.at(1) = principalDir.at(1, indexMid);
                u.at(2) = principalDir.at(2, indexMid);
                u.at(3) = principalDir.at(3, indexMid);
                dyadU.beDyadicProductOf(u, u);
                unitMtrx.beUnitMatrix();
                unitMtrx.times( cos(psi) );
                ux.times( sin(psi) );
                dyadU.times( 1. - cos(psi) );
                rotMtrx.zero();
                rotMtrx.add(unitMtrx);
                rotMtrx.add(ux);
                rotMtrx.add(dyadU);
                crackV.rotatedWith(rotMtrx, 'n');
                for ( int i = 1; i <= 3; i++ ) {
                    principalDir.at(i, indx) = crackV.at(i);
                }

                crackPlaneNormal.rotatedWith(rotMtrx, 'n');
            }
        }


        crackVect.beColumnOf(principalDir, indx);

        status->setCrackVector(crackVect);

        if ( isCrackBandApproachUsed() ) { // le needed only if the crack band approach is used
            le = gp->giveElement()->giveCharacteristicSize(gp, crackPlaneNormal, ecsMethod);
            // store le in corresponding status
            status->setLe(le);
        }

        // compute and store the crack angle (just for postprocessing)
        double ca = M_PI / 2.;
        if ( crackPlaneNormal.at(1) != 0.0 ) {
            ca = atan( crackPlaneNormal.at(2) / crackPlaneNormal.at(1) );
        }

        status->setCrackAngle(ca);

        if ( this->gf != 0. && e0 >= ( wf / le ) ) { // case for a given fracture energy
            OOFEM_WARNING("Fracturing strain %e is lower than the elastic strain e0=%e, possible snap-back. Element number %d, wf %e, le %e", wf / le, e0, gp->giveElement()->giveLabel(), wf, le);
            if ( checkSnapBack ) {
                OOFEM_ERROR("x");
            }
        } else if ( wf == 0. && e0 >= ef ) {
            OOFEM_WARNING( "Fracturing strain ef=%e is lower than the elastic strain e0=%f, possible snap-back. Increase fracturing strain to %f. Element number %d", ef, e0, e0, gp->giveElement()->giveLabel() );
            if ( checkSnapBack ) {
                OOFEM_ERROR("x");
            }
        } else if ( ef == 0. && e0 * le >= wf ) {
            OOFEM_WARNING( "Crack opening at zero stress wf=%f is lower than the elastic displacement w0=%f, possible snap-back. Increase crack opening wf to %f. Element number %d", wf, e0 * le, e0 * le, gp->giveElement()->giveLabel() );
            if ( checkSnapBack ) {
                OOFEM_ERROR("x");
            }
        }
    }
}

double
IsotropicDamageMaterial1 :: give(int aProperty, GaussPoint *gp) const
{
    double answer;
    if ( static_cast< IsotropicDamageMaterial1Status * >( this->giveStatus(gp) )->_giveProperty(aProperty, answer) ) {
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
        return static_cast< MaterialModelMapperInterface * >(this);
    } else {
        return nullptr;
    }
}

void
IsotropicDamageMaterial1::saveContext(DataStream &stream, ContextMode mode)
{
    if ( ( mode & CM_Definition ) ) {
        DynamicInputRecord input;
        this->giveInputRecord(input);
        if ( !stream.write(input.giveRecordInTXTFormat()
                           ) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
IsotropicDamageMaterial1::restoreContext(DataStream &stream, ContextMode mode)
{
    if ( ( mode & CM_Definition ) ) {
        std::string input;
        if ( !stream.read(input) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        OOFEMTXTInputRecord ir(0, input);
        this->initializeFrom(ir);
    }
}


std::unique_ptr<MaterialStatus> 
IsotropicDamageMaterial1 :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<IsotropicDamageMaterial1Status>(gp);
}

MaterialStatus*
IsotropicDamageMaterial1 :: giveStatus(GaussPoint *gp) const
{
    if ( !gp->hasMaterialStatus()) {
            gp->setMaterialStatus(this->CreateStatus(gp));
            this->_generateStatusVariables(gp);
    }
    return static_cast<MaterialStatus*>(gp->giveMaterialStatus());
}


int
IsotropicDamageMaterial1 :: MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep)
{
    int result;
    FloatArray intVal;
    IntArray toMap(3);
    IsotropicDamageMaterial1Status *status = static_cast< IsotropicDamageMaterial1Status * >( this->giveStatus(gp) );


    toMap.at(1) = ( int ) IST_MaxEquivalentStrainLevel;
    toMap.at(2) = ( int ) IST_DamageTensor;
    toMap.at(3) = ( int ) IST_StrainTensor;


    if ( sourceElemSet == NULL ) {
        sourceElemSet = new Set(0, oldd);
        IntArray el;
        // compile source list to contain all elements on old odmain with the same material id
        for ( int i = 1; i <= oldd->giveNumberOfElements(); i++ ) {
            if ( oldd->giveElement(i)->giveCrossSection()->giveMaterial(gp)->giveNumber() == this->giveNumber() ) {
                // add oldd domain element to source list
                el.followedBy(i, 10);
            }
        }
        sourceElemSet->setElementList(el);
    }

    // Set up source element set if not set up by user
    this->mapper.init(oldd, toMap, gp, * sourceElemSet, tStep);

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
        FloatArray sr;
        this->giveReducedSymVectorForm( sr, intVal, gp->giveMaterialMode() );
        status->letTempStrainVectorBe(sr);
    }

#endif
    status->updateYourself(tStep);

    if ( result ) {
        FloatArray sr;
        this->giveReducedSymVectorForm( sr, intVal, gp->giveMaterialMode() );
        status->letTempStrainVectorBe(sr);
    }

    return result;
}




int
IsotropicDamageMaterial1 :: MMI_update(GaussPoint *gp,  TimeStep *tStep, FloatArray *estrain)
{
    int result = 1;
    FloatArray intVal;

    // now update all internal vars accordingly
#ifdef IDM_USE_MAPPEDSTRAIN
    IsotropicDamageMaterial1Status *status = static_cast< IsotropicDamageMaterial1Status * >( this->giveStatus(gp) );
    FloatArray strain = status->giveStrainVector();
    this->giveRealStressVector(intVal, gp, strain, tStep);
#else
    this->giveRealStressVector(intVal, gp, * estrain, tStep);
#endif
    gp->updateYourself(tStep);
    return result;
}


int
IsotropicDamageMaterial1 :: MMI_finish(TimeStep *tStep)
{
    this->mapper.finish(tStep);
    return 1;
}


IsotropicDamageMaterial1Status :: IsotropicDamageMaterial1Status(GaussPoint *g) :
    IsotropicDamageMaterialStatus(g), RandomMaterialStatusExtensionInterface()
{}

Interface *
IsotropicDamageMaterial1Status :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return static_cast< RandomMaterialStatusExtensionInterface * >(this);
    } else {
        return NULL;
    }
}
}     // end namespace oofem
