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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#include "sm/Elements/LatticeElements/latticestructuralelement.h"
#include "gausspoint.h"
#include "crosssection.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"

namespace oofem {
LatticeStructuralElement :: LatticeStructuralElement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{ }


void
LatticeStructuralElement :: computeGlobalRotationMatrix(FloatMatrix &answer, const FloatArray &psi)
{
    if ( psi.giveSize() != 3 ) {
        OOFEM_ERROR("psi param size mismatch");
    }

    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 1.;

    double psiSize = psi.computeNorm();
    if ( psiSize <= 1.e-40 ) {
        return;
    }

    FloatMatrix S(3, 3), SS(3, 3);
    computeSMtrx(S, psi);
    SS.beProductOf(S, S);
    S.times(sin(psiSize) / psiSize);
    SS.times( ( 1. - cos(psiSize) ) / ( psiSize * psiSize ) );

    answer.add(S);
    answer.add(SS);
}


void
LatticeStructuralElement :: computeSMtrx(FloatMatrix &answer, const FloatArray &vec)
{
    if ( vec.giveSize() != 3 ) {
        OOFEM_ERROR("vec param size mismatch");
    }

    answer.resize(3, 3);
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 0.;
    answer.at(1, 2) = -vec.at(3);
    answer.at(1, 3) =  vec.at(2);
    answer.at(2, 1) =  vec.at(3);
    answer.at(2, 3) = -vec.at(1);
    answer.at(3, 1) = -vec.at(2);
    answer.at(3, 2) =  vec.at(1);
}


void
LatticeStructuralElement :: logRotationVector(const FloatMatrix &R, FloatArray &theta)
{
    theta.resize(3);

    double tr = R.at(1, 1) + R.at(2, 2) + R.at(3, 3);
    double cosAngle = 0.5 * ( tr - 1.0 );
    if ( cosAngle >  1.0 ) cosAngle =  1.0;
    if ( cosAngle < -1.0 ) cosAngle = -1.0;
    double angle = acos(cosAngle);

    // Antisymmetric part: w = sin(angle) * axis
    FloatArray w(3);
    w.at(1) = 0.5 * ( R.at(3, 2) - R.at(2, 3) );
    w.at(2) = 0.5 * ( R.at(1, 3) - R.at(3, 1) );
    w.at(3) = 0.5 * ( R.at(2, 1) - R.at(1, 2) );

    if ( angle < 1.e-7 ) {
        theta = w;
    } else if ( angle < M_PI - 1.e-3 ) {
        theta = w;
        theta.times( angle / sin(angle) );
    } else {
        int k = 1;
        if ( R.at(2, 2) > R.at(k, k) ) k = 2;
        if ( R.at(3, 3) > R.at(k, k) ) k = 3;
        double nk = sqrt( 0.5 * ( R.at(k, k) + 1.0 ) );
        FloatArray axis(3);
        axis.at(k) = nk;
        for ( int i = 1; i <= 3; i++ ) {
            if ( i != k ) {
                axis.at(i) = 0.25 * ( R.at(i, k) + R.at(k, i) ) / nk;
            }
        }
        axis.normalize();
        if ( axis.dotProduct(w) < 0.0 ) {
            axis.times(-1.0);
        }
        theta = axis;
        theta.times(angle);
    }
}

void
LatticeStructuralElement :: initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority)
{
    StructuralElement :: initializeFrom(ir, priority);
}


void LatticeStructuralElement :: giveSectionScaleFactors3d(FloatArray &q, GaussPoint *gp)
{


    q.resize(6);
    double A  = this->giveArea(gp);
    double Aq1 = this->giveShearArea1(gp);
    double Aq2 = this->giveShearArea2(gp);
    double I1 = this->giveI1(gp);
    double I2 = this->giveI2(gp);
    double J = this->giveJ(gp);

    q.at(1) = A;
    q.at(2) = Aq1;
    q.at(3) = Aq2;
    q.at(4) = J;
    q.at(5) = I1;
    q.at(6) = I2;
}

void LatticeStructuralElement :: giveSectionScaleFactors2d(FloatArray &q, GaussPoint *gp)
{
    q.resize(3);
    double A  = this->giveArea(gp);
    double I2 = this->giveI2(gp);


    q.at(1) = A;
    q.at(2) = A;
    q.at(3) = I2;
}

void LatticeStructuralElement :: convertStressToResultants3d(FloatArray &S, const FloatArray &sigma, GaussPoint *gp)
{

    S = sigma;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors3d(q, gp);
        for (int i = 1; i <= 6; ++i) {
            S.at(i) *= q.at(i);
        }
    }
}

void LatticeStructuralElement :: convertStressToResultants2d(FloatArray &S, const FloatArray &sigma, GaussPoint *gp)
{
    S = sigma;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors2d(q, gp);

        for (int i = 1; i <= 3; ++i) {
            S.at(i) *= q.at(i);
        }
    }
}





void LatticeStructuralElement :: convertTangentToResultantTangent3d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp)
{

    DS = Dsig;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors3d(q, gp);

        // DS = Q * Dsig. Scale rows
        for (int i = 1; i <= 6; ++i) {
            for (int j = 1; j <= 6; ++j) {
                DS.at(i,j) *= q.at(i);
            }
        }
    }
}


void LatticeStructuralElement :: convertTangentToResultantTangent2d(FloatMatrix &DS, const FloatMatrix &Dsig, GaussPoint *gp)
{
    DS = Dsig;

    auto *mat = static_cast<LatticeStructuralMaterial*>(this->giveCrossSection()->giveMaterial(gp));
    if ( mat->giveLatticeResponseType() == LatticeStructuralMaterial::LRT_StressBased ) {
        FloatArray q;
        this->giveSectionScaleFactors2d(q, gp);
        // DS = Q * Dsig. Scale rows
        for (int i = 1; i <= 3; ++i) {
            for (int j = 1; j <= 3; ++j) {
                DS.at(i,j) *= q.at(i);
            }
        }
    }
    }


void
LatticeStructuralElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;

    if ( !this->isActivated(tStep) ) {
        this->computeBmatrixAt(gp, b);
        answer.resize(b.giveNumberOfRows());
        answer.zero();
        return;
    }

    this->computeBmatrixAt(gp, b);
    this->computeVectorOf(VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }
    answer.beProductOf(b, u);
    answer.times(1./giveLength());
    }



void
LatticeStructuralElement :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralElement :: printOutputAt(file, tStep);
}



} // end namespace oofem
