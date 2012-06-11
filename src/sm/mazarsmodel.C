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

#include "mazarsmodel.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
#define _MAZARS_MODEL_ITER_TOL 1.e-15
#define _MAZARS_MODEL_MAX_ITER 400.

MazarsMaterial :: MazarsMaterial(int n, Domain *d) : IsotropicDamageMaterial1(n, d)
    //
    // constructor
    //
{
    // force the loading / unloading control according to max. reached damage level
    llcriteria = idm_damageLevelCR;
}


MazarsMaterial :: ~MazarsMaterial()
//
// destructor
//
{ }

IRResultType
MazarsMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int ver;

    // Note: IsotropicDamageMaterial1 :: initializeFrom is not activated
    // because we do not always read ef and the equivalent strain type
    // cannot be selected
    // IsotropicDamageMaterial1 :: initializeFrom(ir);

    this->equivStrainType = EST_Mazars;

    IsotropicDamageMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);

    linearElasticMaterial->initializeFrom(ir);
    // E and nu are made available for direct access
    IR_GIVE_FIELD(ir, E, IFT_IsotropicLinearElasticMaterial_e, "e");
    IR_GIVE_FIELD(ir, nu, IFT_IsotropicLinearElasticMaterial_n, "n");

    ver = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, ver, IFT_MazarsMaterial_version, "version"); // Macro
    if ( ver == 1 ) {
        this->modelVersion = maz_modTension;
    } else if ( ver == 0 ) {
        this->modelVersion = maz_original;
    } else {
        _error("instanciateFrom: unknown version");
    }

    IR_GIVE_FIELD(ir, this->e0, IFT_MazarsMaterial_e0, "e0");
    IR_GIVE_FIELD(ir, this->Ac, IFT_MazarsMaterial_ac, "ac");

    this->Bc = ( Ac - 1.0 ) / ( Ac * e0 ); // default value, ensures smooth curve
    IR_GIVE_OPTIONAL_FIELD(ir, this->Bc, IFT_MazarsMaterial_bc, "bc");

    beta = 1.06;
    IR_GIVE_OPTIONAL_FIELD(ir, beta, IFT_MazarsMaterial_beta, "beta");

    if ( this->modelVersion == maz_original ) {
        IR_GIVE_FIELD(ir, this->At, IFT_MazarsMaterial_at, "at");
        IR_GIVE_FIELD(ir, this->Bt, IFT_MazarsMaterial_bt, "bt");
    } else if ( this->modelVersion == maz_modTension ) {
        // in case of modified model read ef instead of At, Bt
        IR_GIVE_FIELD(ir, this->ef, IFT_MazarsMaterial_ef, "ef");
    }

    // ask for optional "reference length"
    hReft = hRefc = 0.; // default values 0 => no adjustment for element size is used
    IR_GIVE_OPTIONAL_FIELD(ir, this->hReft, IFT_MazarsMaterial_hreft, "hreft");
    IR_GIVE_OPTIONAL_FIELD(ir, this->hRefc, IFT_MazarsMaterial_hrefc, "hrefc");

    this->mapper.initializeFrom(ir);

    return IRRT_OK;
}


void
MazarsMaterial :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    int i;
    double posNorm = 0.0;
    FloatArray principalStrains, strainb;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }

    crossSection->giveFullCharacteristicVector(strainb, gp, strain);
    // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
    int ndim = giveNumberOfSpatialDimensions(gp);
    if ( ndim == 2 ) {
        strainb.at(3) = -nu * ( strainb.at(1) + strainb.at(2) ) / ( 1. - nu );
    } else if ( ndim == 1 ) {
        strainb.at(2) = -nu *strainb.at(1);
        strainb.at(3) = -nu *strainb.at(1);
    }

    if ( ndim == 1 ) {
        principalStrains.resize(3);
        for ( i = 1; i <= 3; i++ ) {
            principalStrains.at(i) = strainb.at(i);
        }
    } else {
        this->computePrincipalValues(principalStrains, strainb, principal_strain);
    }

    /*
     * // simple check
     * double i1 = strainb.at(1)+strainb.at(2)+strainb.at(3) - principalStrains.at(1)-principalStrains.at(2)-principalStrains.at(3);
     * double i2 = strainb.at(1)*strainb.at(3)+strainb.at(2)*strainb.at(3)+strainb.at(1)*strainb.at(2) -
     * 0.25*(strainb.at(4)*strainb.at(4)+strainb.at(5)*strainb.at(5)+strainb.at(6)*strainb.at(6)) -
     * principalStrains.at(1)*principalStrains.at(3)+principalStrains.at(2)*principalStrains.at(3)+principalStrains.at(1)*principalStrains.at(2);
     * if ((fabs(i1) > 1.e-6) || (fabs(i2) > 1.e-6)) {
     * printf("v");
     * }
     * // end simple check
     */
    for ( i = 1; i <= 3; i++ ) {
        if ( principalStrains.at(i) > 0.0 ) {
            posNorm += principalStrains.at(i) * principalStrains.at(i);
        }
    }

    kappa = sqrt(posNorm);
}

/*
 * void
 * MazarsMaterial :: giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
 *                                                 MatResponseMode rMode,
 *                                                 GaussPoint *gp, TimeStep *atTime)
 * {
 *  //
 *  // return Elastic Stiffness matrix for normal Stresses
 *  LinearElasticMaterial *lMat = this->giveLinearElasticMaterial();
 *  FloatMatrix de;
 *  int i, j;
 *
 *  answer.resize(3, 3);
 *  lMat->give3dMaterialStiffnessMatrix(de, FullForm, rMode, gp, atTime);
 *  // copy first 3x3 submatrix to answer
 *  for ( i = 1; i <= 3; i++ ) {
 *      for ( j = 1; j <= 3; j++ ) {
 *          answer.at(i, j) = de.at(i, j);
 *      }
 *  }
 * }
 */

int
MazarsMaterial :: giveNumberOfSpatialDimensions(GaussPoint *gp)
{
    if ( gp->giveMaterialMode() == _1dMat ) {
        return 1;
    } else if ( gp->giveMaterialMode() == _PlaneStress ) {
        return 2;
    } else {
        return 3;
    }
}

void
MazarsMaterial :: giveNormalBlockOfElasticCompliance(FloatMatrix &answer, GaussPoint *gp)
{
    int ndim = giveNumberOfSpatialDimensions(gp);
    answer.resize(ndim, ndim);
    for ( int i = 1; i <= ndim; i++ ) {
        for ( int j = 1; j <= ndim; j++ ) {
            if ( i == j ) {
                answer.at(i, j) = 1. / E;
            } else {
                answer.at(i, j) = -nu / E;
            }
        }
    }
}

void
MazarsMaterial :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    int i; // positive_flag = 0, negat_count = 0;
    FloatMatrix de, ce;
    FloatArray sigp, epsti, epsi;
    double gt, gc, alpha_t, alpha_c, alpha, eqStrain2;

    if ( kappa <= this->e0 ) { // strain below damage threshold
        omega = 0.0;
        return;
    }

    // strain above damage threshold

    int ndim = giveNumberOfSpatialDimensions(gp);
    if ( !strain.isEmpty() ) {
        this->computePrincipalValues(epsi, strain, principal_strain);
    } else {
        epsi.resize(ndim);
        epsi.zero();
    }

    // construct the normal block of elastic compliance matrix
    giveNormalBlockOfElasticCompliance(ce, gp);
    // compute the normal block of elastic stiffness matrix
    de.beInverseOf(ce);

    // compute principal stresses
    sigp.beProductOf(de, epsi);
    // take positive part
    for ( i = 1; i <= ndim; i++ ) {
        if ( sigp.at(i) < 0. ) {
            sigp.at(i) = 0.;
        }
    }

    // compute principal strains due to positive stresses
    epsti.beProductOf(ce, sigp);

    // extend strains to 3 dimensions
    if ( ndim == 2 ) {
        epsi.resize(3);
        epsti.resize(3);
        epsi.at(3) = -nu * ( epsi.at(1) + epsi.at(2) ) / ( 1. - nu );
        epsti.at(3) = -nu * ( epsti.at(1) + epsti.at(2) ) / ( 1. - nu );
    } else if ( ndim == 1 ) {
        epsi.resize(3);
        epsti.resize(3);
        epsi.at(2) = epsi.at(3) = -nu *epsi.at(1);
        epsti.at(2) = epsti.at(3) = -nu *epsti.at(1);
    }

    /* the following section "improves" biaxial compression
     * // but it may lead to convergence problems, probably due to the condition > 1.e-6
     * // therefore it was commented out by Milan Jirasek on 22 Nov 2009
     *
     * positive_flag = negat_count = 0;
     *    for ( i = 1; i <= 3; i++ ) {
     *        if ( sigp.at(i) > 1.e-6 ) {
     *            positive_flag = 1;
     *            break;
     *        } else if ( sigp.at(i) < 1.e-6 ) {
     *            negat_count++;
     *        }
     *    }
     *
     *    if ( ( positive_flag == 0 ) && ( negat_count > 1 ) ) {
     *        // adjust equivalent strain to improve biaxial compression
     *        double f = 0.0, g = 0.0;
     *        for ( i = 1; i <= 3; i++ ) {
     *            if ( sigp.at(i) < 0 ) {
     *                f += sigp.at(i) * sigp.at(i);
     *                g += fabs( sigp.at(i) );
     *            }
     *        }
     *
     *        f = sqrt(f) / g;
     *        kappa *= f;
     *    }
     *
     *    // test adjusted kappa
     *    if ( kappa < this->e0 ) {
     *        omega = 0.0;
     *        return;
     *    }
     *    // end of section that improves biaxial compression
     */

    // evaluation of damage functions
    gt = computeGt(kappa, gp);
    gc = computeGc(kappa, gp);

    // evaluation of factors alpha_t and alpha_c
    alpha = 0.0;
    eqStrain2 = 0.0;
    for ( i = 1; i <= 3; i++ ) {
        if ( epsi.at(i) > 0.0 ) {
            eqStrain2 += epsi.at(i) * epsi.at(i);
            alpha += epsti.at(i) * epsi.at(i);
        }
    }

    if ( eqStrain2 > 0. ) {
        alpha /= eqStrain2;
    }

    if ( alpha > 1. ) {
        alpha = 1.;
    } else if ( alpha < 0. ) {
        alpha = 0.;
    }

    if ( this->beta == 1. ) {
        alpha_t = alpha;
        alpha_c = 1. - alpha;
    } else if ( alpha <= 0. )   {
        alpha_t = 0.;
        alpha_c = 1.;
    } else if ( alpha < 1. ) {
        alpha_t = pow(alpha, this->beta);
        alpha_c = pow( ( 1. - alpha ), this->beta );
    } else   {
        alpha_t = 1.;
        alpha_c = 0.;
    }

    omega = alpha_t * gt + alpha_c * gc;
    if ( omega > 1.0 ) {
        omega = 1.0;
    }
}

// evaluation of tensile damage
double MazarsMaterial :: computeGt(double kappa, GaussPoint *gp)
{
    double gt;
    if ( hReft <= 0. ) { // material law considered as independent of element size
        if ( this->modelVersion == maz_modTension ) { // exponential softening
            gt = 1.0 - ( this->e0 / kappa ) * exp( ( this->e0 - kappa ) / this->ef );
        } else { // maz_original
            gt = 1.0 - ( 1.0 - this->At ) * this->e0 / kappa - this->At *exp( -this->Bt * ( kappa - this->e0 ) );
        }

        return gt;
    }

    // tension objectivity (material law dependent on element size)
    int nite = 0;
    double aux, hCurrt, kappaRefT, dgt, R;
    MazarsMaterialStatus *status = ( MazarsMaterialStatus * ) this->giveStatus(gp);
    hCurrt = status->giveLe();
    kappaRefT = kappa;
    do {
        if ( this->modelVersion == maz_modTension ) {
            gt = 1.0 - ( this->e0 / kappaRefT ) * exp( ( this->e0 - kappaRefT ) / this->ef );
            dgt = this->e0 / kappaRefT / kappaRefT + this->e0 / kappaRefT / this->ef;
            dgt *= exp( ( this->e0 - kappaRefT ) / this->ef );
        } else { // maz_original
            gt   = 1.0 - ( 1.0 - this->At ) * this->e0 / kappaRefT - this->At *exp( -this->Bt * ( kappaRefT - this->e0 ) );
            dgt  = ( 1.0 - this->At ) * this->e0 / kappaRefT / kappaRefT + this->At *exp( -this->Bt * ( kappaRefT - this->e0 ) ) * this->Bt;
        }

        aux = 1 + gt * ( hReft / hCurrt - 1.0 );
        R   = kappaRefT * aux - kappa;
        if ( fabs(R) <= _MAZARS_MODEL_ITER_TOL ) {
            if ( ( gt < 0. ) || ( gt > 1.0 ) ) {
                _error("computeDamageParam: gt out of range ");
            }

            return gt * ( hReft * kappaRefT ) / ( hCurrt * kappa );
        }

        aux += dgt * kappaRefT * ( hReft / hCurrt - 1.0 );
        kappaRefT += -R / aux;
    } while ( nite++ < _MAZARS_MODEL_MAX_ITER );

    _error("computeDamageParam: tension objectivity iteration internal error - no convergence");
    return 0.; // just to make the compiler happy
}

// evaluation of compression damage
double MazarsMaterial :: computeGc(double kappa, GaussPoint *gp)
{
    double gc;
    if ( hRefc <= 0. ) { // material law considered as independent of element size
        gc = 1.0 - ( 1.0 - this->Ac ) * this->e0 / kappa - this->Ac *exp( -this->Bc * ( kappa - this->e0 ) );
        return gc;
    }

    // compression objectivity (material law dependent on element size)
    int nite = 0;
    double aux, hCurrc, kappaRefC, dgc, R;
    MazarsMaterialStatus *status = ( MazarsMaterialStatus * ) this->giveStatus(gp);
    hCurrc = status->giveLec();
    kappaRefC = kappa;
    do {
        gc   = 1.0 - ( 1.0 - this->Ac ) * this->e0 / kappaRefC - this->Ac *exp( -this->Bc * ( kappaRefC - this->e0 ) );
        aux = 1 + gc * ( hRefc / hCurrc - 1.0 );
        R   = kappaRefC * aux - kappa;
        if ( fabs(R) <= _MAZARS_MODEL_ITER_TOL ) {
            return gc * ( hRefc * kappaRefC ) / ( hCurrc * kappa );
        }

        dgc  = ( 1.0 - this->Ac ) * this->e0 / kappaRefC / kappaRefC + this->Ac *exp( -this->Bc * ( kappaRefC - this->e0 ) ) * this->Bc;
        aux += dgc * kappaRefC * ( hRefc / hCurrc - 1.0 );
        kappaRefC += -R / aux;
    } while ( nite++ < _MAZARS_MODEL_MAX_ITER );

    _error("computeDamageParam: compression objectivity iteration internal error - no convergence");
    return 0.; // just to make the compiler happy
}

void
MazarsMaterial :: initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp)
{
    int i, indmin = 1, indmax = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3), fullstrain;
    FloatMatrix principalDir(3, 3);
    MazarsMaterialStatus *status = ( MazarsMaterialStatus * ) this->giveStatus(gp);
    //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();

    // crossSection->giveFullCharacteristicVector(fullstrain, gp, totalStrainVector);

    if ( ( kappa > this->e0 ) && ( status->giveDamage() == 0. ) ) {
        //printf (":");


        if ( gp->giveMaterialMode() == _1dMat ) {
            crackPlaneNormal.zero();
            crackPlaneNormal.at(1) = 1.0;
            le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
            status->setLe(le);
            status->setLec(le);
            return;
        }


        if ( gp->giveMaterialMode() == _PlaneStress ) {
            principalDir.resize(2, 2);
        }

        this->computePrincipalValDir(principalStrains, principalDir, totalStrainVector, principal_strain);
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            if ( principalStrains.at(1) > principalStrains.at(2) ) {
                indmax = 1;
                indmin = 2;
            } else {
                indmax = 2;
                indmin = 1;
            }
        } else {
            // find index of max and min positive principal strains
            for ( i = 2; i <= 3; i++ ) {
                if ( principalStrains.at(i) > principalStrains.at(indmax) ) {
                    indmax = i;
                }

                if ( principalStrains.at(i) < principalStrains.at(indmin) ) {
                    indmin = i;
                }
            }
        }

        for ( i = 1; i <= principalStrains.giveSize(); i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indmax);
        }

        le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
        // remember le in cooresponding status for tension
        status->setLe(le);

        for ( i = 1; i <= principalStrains.giveSize(); i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indmin);
        }

        le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
        // remember le in cooresponding status for compression
        status->setLec(le);

        // printf ("les: %e %e\n", status->giveLe(), status->giveLec());
    }

    //status->setLe(hReft); status->setLec(hRefc);
}








MazarsMaterialStatus :: MazarsMaterialStatus(int n, Domain *d, GaussPoint *g) :
    IsotropicDamageMaterial1Status(n, d, g)
{
    lec = 0.0;
}

contextIOResultType
MazarsMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterial1Status :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& lec, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
MazarsMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterial1Status :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& lec, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
} // end namespace oofem
