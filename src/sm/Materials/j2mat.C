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

#include "j2mat.h"
#include "isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(J2Mat);

J2Mat :: J2Mat(int n, Domain *d) : MPlasticMaterial2(n, d)
{
    //
    // constructor
    //
    kinematicHardeningFlag = isotropicHardeningFlag = 0;
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    this->nsurf = 1;
}

J2Mat :: ~J2Mat()
{ }

IRResultType
J2Mat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    double value;

    result = MPlasticMaterial2 :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;
    
    result = linearElasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, value, _IFT_J2Mat_ry);
    k = value / sqrt(3.0);

    kinematicModuli = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, kinematicModuli, _IFT_J2Mat_khm);

    isotropicModuli = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, isotropicModuli, _IFT_J2Mat_ihm);

    if ( fabs(kinematicModuli) > 1.e-12 ) {
        kinematicHardeningFlag = 1;
    }

    if ( fabs(isotropicModuli) > 1.e-12 ) {
        isotropicHardeningFlag = 1;
    }

    int rma = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, rma, _IFT_J2Mat_rma);
    if ( rma == 0 ) {
        this->rmType = mpm_ClosestPoint;
    } else {
        this->rmType = mpm_CuttingPlane;
    }


    return IRRT_OK;
}


MaterialStatus *
J2Mat :: CreateStatus(GaussPoint *gp) const
{
    return new MPlasticMaterial2Status(1, this->giveDomain(), gp, this->giveSizeOfReducedHardeningVarsVector(gp));
}

int
J2Mat :: giveSizeOfFullHardeningVarsVector()
{
    /* Returns the size of hardening variables vector */
    int size = 0;

    if ( kinematicHardeningFlag ) {
        size += 6;                      /* size of full stress vector */
    }

    if ( isotropicHardeningFlag ) {
        size += 1;                      /* scalar value */
    }

    return size;
}

int
J2Mat :: giveSizeOfReducedHardeningVarsVector(GaussPoint *gp) const
{
    /* Returns the size of hardening variables vector */
    int size = 0;

    if ( kinematicHardeningFlag ) {
        size += StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    }

    if ( isotropicHardeningFlag ) {
        size += 1;                      /* scalar value */
    }

    return size;
}


double
J2Mat :: computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                             const FloatArray &strainSpaceHardeningVars)
{
    double f;
    FloatArray helpVector, backStress;

    if ( this->kinematicHardeningFlag ) {
        if ( stressVector.isNotEmpty() ) {
            helpVector = stressVector;
            this->giveStressBackVector(backStress, gp, strainSpaceHardeningVars);
            helpVector.add(backStress);
        } else {
            return -k;
        }
    } else {
        helpVector = stressVector;
    }

    f = sqrt( this->computeJ2InvariantAt(helpVector) );
    return f + sqrt(1. / 3.) * this->giveIsotropicHardeningVar(gp, strainSpaceHardeningVars) - this->k;
}


void
J2Mat :: computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                     const FloatArray &strainSpaceHardeningVars)
{
    /* stress gradient of yield function in full stress - strain space */

    double f, ax, ay, az, sx, sy, sz;
    FloatArray helpVector, backStress;

    answer.resize(6);
    answer.zero();
    if ( this->kinematicHardeningFlag ) {
        if ( stressVector.isNotEmpty() ) {
            this->giveStressBackVector(backStress, gp, strainSpaceHardeningVars);
            helpVector = stressVector;
            helpVector.add(backStress);
        } else {
            return;
        }
    } else {
        helpVector = stressVector;
    }

    f = sqrt( this->computeJ2InvariantAt(helpVector) );
    // check for yield value zero value
    if ( fabs(f) < 1.e-6 ) {
        return;
    }

    ax = helpVector.at(1);
    ay = helpVector.at(2);
    az = helpVector.at(3);

    sx = ( 2. / 3. ) * ax - ( 1. / 3. ) * ay - ( 1. / 3. ) * az;
    sy = ( 2. / 3. ) * ay - ( 1. / 3. ) * ax - ( 1. / 3. ) * az;
    sz = ( 2. / 3. ) * az - ( 1. / 3. ) * ay - ( 1. / 3. ) * ax;

    answer.at(1) = 0.5 * sx / f;
    answer.at(2) = 0.5 * sy / f;
    answer.at(3) = 0.5 * sz / f;
    answer.at(4) = helpVector.at(4) / f;
    answer.at(5) = helpVector.at(5) / f;
    answer.at(6) = helpVector.at(6) / f;
}


void
J2Mat :: computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &stress, const FloatArray &dlambda,
                                             const FloatArray &dplasticStrain, const IntArray &activeConditionMap)
{
    int size = this->giveSizeOfReducedHardeningVarsVector(gp);
    answer.resize(size);

    if ( this->kinematicHardeningFlag ) {
        int sizer = dplasticStrain.giveSize();
        double coeff = sqrt(2.) * ( 2. / 3. );
        for ( int i = 1; i <= sizer; i++ ) {
            answer.at(i) = dplasticStrain.at(i) * coeff;
        }
    }

    if ( isotropicHardeningFlag ) {
        answer.at(size) = sqrt(1. / 3.) * dlambda.at(1);
    }
}


void
J2Mat :: computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector,
                                const FloatArray &strainSpaceHardeningVariables)
{
    int kcount = 0, size = this->giveSizeOfReducedHardeningVarsVector(gp);
    FloatArray reducedKinematicGrad;

    if ( !hasHardening() ) {
        answer.clear();
        return;
    }

    answer.resize(size);

    /* kinematic hardening variables first */
    if ( this->kinematicHardeningFlag ) {
        this->computeReducedStressGradientVector(reducedKinematicGrad, ftype, isurf, gp, fullStressVector, strainSpaceHardeningVariables);
        kcount = reducedKinematicGrad.giveSize();
        for ( int i = 1; i <= kcount; i++ ) {
            answer.at(i) = ( -1.0 ) * this->kinematicModuli * reducedKinematicGrad.at(i);
        }
    }

    if ( this->isotropicHardeningFlag ) {
        answer.at(size) = ( -1.0 ) * this->isotropicModuli;
    }
}

void
J2Mat :: computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap,
                                                  const FloatArray &fullStressVector,
                                                  const FloatArray &strainSpaceHardeningVars,
                                                  const FloatArray &gamma)
{
    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    answer.resize(giveSizeOfReducedHardeningVarsVector(gp), rsize);
    answer.zero();

    if ( this->kinematicHardeningFlag ) {
        double coeff = sqrt(2.) * ( 2. / 3. ) * gamma.at(1);
        FloatMatrix h;

        this->computeReducedSSGradientMatrix(h, 1, gp, fullStressVector, strainSpaceHardeningVars);
        for ( int i = 1; i <= rsize; i++ ) {
            for ( int j = 1; j <= rsize; j++ ) {
                answer.at(i, j) = coeff * h.at(i, j);
            }
        }
    }
}

void
J2Mat :: computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf,
                                                const IntArray &activeConditionMap,
                                                const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVars,
                                                const FloatArray &gamma)
{
    int size = this->giveSizeOfReducedHardeningVarsVector(gp);
    answer.resize(size, 1);

    if ( this->kinematicHardeningFlag ) {
        int rsize;
        FloatArray loadGradSigVec;
        this->computeReducedStressGradientVector(loadGradSigVec, loadFunction, 1, gp, fullStressVector,
                                                 strainSpaceHardeningVars);
        rsize = loadGradSigVec.giveSize();
        for ( int i = 1; i <= rsize; i++ ) {
            answer.at(i, 1) = loadGradSigVec.at(i);
        }

        answer.times( sqrt(2.) * ( 2. / 3. ) );
    }

    if ( isotropicHardeningFlag ) {
        answer.at(size, 1) = sqrt(1. / 3.);
    }
}


void
J2Mat :: computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int isurf, GaussPoint *gp, const FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVars)
{
    int size;
    int imask, jmask;
    FloatArray helpVector, backStress, df(6);
    IntArray mask;
    double f, f32, f12, ax, ay, az;

    StructuralMaterial :: giveInvertedVoigtVectorMask( mask, gp->giveMaterialMode() );
    size = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );

    gradientMatrix.resize(size, size);
    gradientMatrix.zero();


    if ( fullStressVector.giveSize() != 0 ) {
        /* kinematic hardening variables first */
        if ( this->kinematicHardeningFlag ) {
            this->giveStressBackVector(backStress, gp, strainSpaceHardeningVars);
            helpVector = fullStressVector;
            helpVector.add(backStress);
        } else {
            helpVector = fullStressVector;
        }

        f = this->computeJ2InvariantAt(helpVector);
        f12 = sqrt(f);
        f32 = pow(f, 3. / 2.);

        ax = helpVector.at(1);
        ay = helpVector.at(2);
        az = helpVector.at(3);

        df.at(1) = ( 2. / 3. ) * ax - ( 1. / 3. ) * ay - ( 1. / 3. ) * az;
        df.at(2) = ( 2. / 3. ) * ay - ( 1. / 3. ) * ax - ( 1. / 3. ) * az;
        df.at(3) = ( 2. / 3. ) * az - ( 1. / 3. ) * ay - ( 1. / 3. ) * ax;
        df.at(4) = 2. * helpVector.at(4);
        df.at(5) = 2. * helpVector.at(5);
        df.at(6) = 2. * helpVector.at(6);

        for ( int i = 1; i <= 3; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( int j = i; j <= 3; j++ ) {
                if ( ( jmask = mask.at(j) ) == 0 ) {
                    continue;
                }

                if ( i == j ) {
                    gradientMatrix.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * ( 4. / 6 );
                } else {
                    gradientMatrix.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * ( -2. / 6 );
                    gradientMatrix.at(jmask, imask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * ( -2. / 6 );
                }
            }
        }

        for ( int i = 1; i <= 3; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( int j = 4; j <= 6; j++ ) {
                if ( ( jmask = mask.at(j) ) == 0 ) {
                    continue;
                }

                gradientMatrix.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
                gradientMatrix.at(jmask, imask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
            }
        }

        for ( int i = 4; i <= 6; i++ ) {
            if ( ( imask = mask.at(i) ) == 0 ) {
                continue;
            }

            for ( int j = i; j <= 6; j++ ) {
                if ( ( jmask = mask.at(j) ) == 0 ) {
                    continue;
                }

                if ( i == j ) {
                    gradientMatrix.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j) + 0.5 * ( 1. / f12 ) * 2.;
                } else {
                    gradientMatrix.at(imask, jmask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
                    gradientMatrix.at(jmask, imask) = -( 1. / 4. ) * ( 1. / f32 ) * df.at(i) * df.at(j);
                }
            }
        }
    }
}


void
J2Mat :: computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVariables)
{
    // something will be here for k1 vector
    int size = giveSizeOfReducedHardeningVarsVector(gp);
    FloatMatrix helpMat;
    gradientMatrix.resize(StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ), size);
    gradientMatrix.zero();

    if ( this->kinematicHardeningFlag ) {
        int kcount;
        this->computeReducedSSGradientMatrix(helpMat, i, gp, fullStressVector, strainSpaceHardeningVariables);
        helpMat.times( ( -1.0 ) * this->kinematicModuli );
        kcount = helpMat.giveNumberOfRows();
        for ( int ii = 1; ii <= kcount; ii++ ) {
            for ( int j = 1; j <= kcount; j++ ) {
                gradientMatrix.at(ii, j) = helpMat.at(ii, j);
            }
        }
    }
}


int
J2Mat :: hasHardening()
{
    return ( this->kinematicHardeningFlag || this->isotropicHardeningFlag );
}


double
J2Mat :: computeJ2InvariantAt(const FloatArray &stressVector)
{
    double answer;
    double v1, v2, v3;

    if ( stressVector.isEmpty() ) {
        return 0.0;
    }

    v1 = ( ( stressVector.at(1) - stressVector.at(2) ) * ( stressVector.at(1) - stressVector.at(2) ) );
    v2 = ( ( stressVector.at(2) - stressVector.at(3) ) * ( stressVector.at(2) - stressVector.at(3) ) );
    v3 = ( ( stressVector.at(3) - stressVector.at(1) ) * ( stressVector.at(3) - stressVector.at(1) ) );

    answer = ( 1. / 6. ) * ( v1 + v2 + v3 ) + stressVector.at(4) * stressVector.at(4) +
    stressVector.at(5) * stressVector.at(5) + stressVector.at(6) * stressVector.at(6);

    return answer;
}


void
J2Mat :: giveStressBackVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &strainSpaceHardeningVars)
{
    /* returns part of hardening vector corresponding to kinematic hardening */
    if ( this->kinematicHardeningFlag ) {
        IntArray mask;
        int isize;

        answer.resize(6);
        StructuralMaterial :: giveVoigtSymVectorMask( mask, gp->giveMaterialMode() );
        isize = mask.giveSize();
        //int rSize = this->giveSizeOfReducedHardeningVarsVector(gp);

        /* kinematic hardening variables are first */
        for ( int i = 1; i <= isize; i++ ) {
            answer.at( mask.at(i) ) = ( -1.0 ) * this->kinematicModuli * strainSpaceHardeningVars.at(i);
        }
    } else {
        answer.clear();
        return;
    }
}


double
J2Mat :: giveIsotropicHardeningVar(GaussPoint *gp, const FloatArray &strainSpaceHardeningVars)
{
    /* returns value in  hardening vector corresponding to isotropic hardening */
    if ( !isotropicHardeningFlag ) {
        return 0.;
    } else {
        int rSize = this->giveSizeOfReducedHardeningVarsVector(gp);

        return ( -1.0 ) * this->isotropicModuli * strainSpaceHardeningVars.at(rSize);
    }
}
} // end namespace oofem
