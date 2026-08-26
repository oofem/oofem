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

#include "domain.h"
#include "lattice3dnl.h"
#include "lattice3d.h"
#include "Materials/latticematstatus.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "mathfem.h"
#include "latticestructuralelement.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "lm/CrossSections/latticecrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "sm/Materials/structuralmaterial.h"
#endif

namespace oofem {
REGISTER_Element(Lattice3dNL);

Lattice3dNL :: Lattice3dNL(int n, Domain *aDomain) : Lattice3d(n, aDomain)
{
}

Lattice3dNL :: ~Lattice3dNL()
{}


    void
    Lattice3dNL::computeNLBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, TimeStep *tStep)
    {
        answer.resize(6, 12);
        answer.zero();

        FloatArray coordA(3), coordB(3), l1(3), l2(3), help(3), gpCoords(3);
        coordA = giveNode(1)->giveCoordinates();
        coordB = giveNode(2)->giveCoordinates();
        giveGpCoordinates(gpCoords, aGaussPoint);

        l1 = gpCoords-coordA;
        l2 = coordB-gpCoords;


        FloatMatrix rotationMatrix(3, 3);
        FloatArray spin(3);
        FloatArray uTotal(12);
        this->computeVectorOf(VM_Total, tStep, uTotal);
        spin.at(1) = uTotal.at(4);
        spin.at(2) = uTotal.at(5);
        spin.at(3) = uTotal.at(6);
        this->computeGlobalRotationMatrix(rotationMatrix, spin);
        FloatArray c1(3);
        c1.beProductOf(rotationMatrix, l1);


        spin.at(1) = uTotal.at(10);
        spin.at(2) = uTotal.at(11);
        spin.at(3) = uTotal.at(12);
        this->computeGlobalRotationMatrix(rotationMatrix, spin);
        FloatArray c2(3);
        c2.beProductOf(rotationMatrix, l2);

        //Normal displacement jump in x-direction
        //First node
        answer.at(1, 1) = -1.;
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 0.;
        answer.at(1, 4) = 0.;
        answer.at(1, 5) = -c1.at(3);
        answer.at(1, 6) = c1.at(2);
        //Second node
        answer.at(1, 7) = 1.;
        answer.at(1, 8) = 0.;
        answer.at(1, 9) = 0.;
        answer.at(1, 10) = 0.;
        answer.at(1, 11) = -c2.at(3);
        answer.at(1, 12) = c2.at(2);

        //Shear displacement jump in y-plane
        //first node
        answer.at(2, 1) = 0.;
        answer.at(2, 2) = -1.;
        answer.at(2, 3) =  0.;
        answer.at(2, 4) = c1.at(3);
        answer.at(2, 5) = 0;
        answer.at(2, 6) = -c1.at(1);
        //Second node
        answer.at(2, 7) = 0.;
        answer.at(2, 8) = 1.;
        answer.at(2, 9) =  0.;
        answer.at(2, 10) = c2.at(3);
        answer.at(2, 11) = 0;
        answer.at(2, 12) = -c2.at(1);

        //Shear displacement jump in z-plane
        //first node
        answer.at(3, 1) = 0.;
        answer.at(3, 2) = 0.;
        answer.at(3, 3) = -1.;
        answer.at(3, 4) = -c1.at(2);
        answer.at(3, 5) = c1.at(1);
        answer.at(3, 6) = 0.;
        //Second node
        answer.at(3, 7) = 0.;
        answer.at(3, 8) = 0.;
        answer.at(3, 9) =  1.;
        answer.at(3, 10) = -c2.at(2);
        answer.at(3, 11) = c2.at(1);
        answer.at(3, 12) = 0.;

        //Rotation around x-axis
        //First node
        answer.at(4, 1) = 0.;
        answer.at(4, 2) = 0;
        answer.at(4, 3) = 0.;
        answer.at(4, 4) = -1.;
        answer.at(4, 5) = 0.;
        answer.at(4, 6) = 0.;
        //Second node
        answer.at(4, 7) = 0.;
        answer.at(4, 8) = 0.;
        answer.at(4, 9) = 0.;
        answer.at(4, 10) = 1.;
        answer.at(4, 11) = 0.;
        answer.at(4, 12) = 0.;

        //Rotation around y-axis
        //First node
        answer.at(5, 1) = 0.;
        answer.at(5, 2) = 0.;
        answer.at(5, 3) = 0.;
        answer.at(5, 4) = 0.;
        answer.at(5, 5) = -1.;
        answer.at(5, 6) = 0.;
        //Second node
        answer.at(5, 7) = 0.;
        answer.at(5, 8) = 0.;
        answer.at(5, 9) =  0.;
        answer.at(5, 10) = 0.;
        answer.at(5, 11) = 1.;
        answer.at(5, 12) = 0.;

        //Rotation around z-axis
        //First node
        answer.at(6, 1) = 0.;
        answer.at(6, 2) = 0.;
        answer.at(6, 3) = 0.;
        answer.at(6, 4) = 0.;
        answer.at(6, 5) = 0.;
        answer.at(6, 6) = -1.;
        //Second node
        answer.at(6, 7) = 0.;
        answer.at(6, 8) = 0.;
        answer.at(6, 9) =  0.;
        answer.at(6, 10) = 0.;
        answer.at(6, 11) = 0.;
        answer.at(6, 12) = 1.;

        return;
    }



    void Lattice3dNL::updateYourself(TimeStep *tStep)
    // Updates the receiver at end of step.
    {
      Lattice3d::updateYourself(tStep);
    }

    void
    Lattice3dNL::initForNewStep()
    //initializes receiver to new time step or can be used
    //if current time step must be restarted
    {
        Lattice3d::initForNewStep();
    }


    void
    Lattice3dNL::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                             TimeStep *tStep)
    {
        FloatMatrix d, ds, bt, db, b, contrib;
        this->length = giveLength();

        FloatArray u;
        this->computeVectorOf(VM_Total, tStep, u);

        FloatArray coordA(3), coordB(3), coordGp(3);
        coordA = giveNode(1)->giveCoordinates();
        coordB = giveNode(2)->giveCoordinates();

        answer.resize(12, 12);
        answer.zero();

        IntegrationRule *iRule = integrationRulesArray [ 0 ].get();
        const int nGp = iRule->giveNumberOfIntegrationPoints();
        for ( int i = 0; i < nGp; ++i ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            this->giveGpCoordinates(coordGp, gp);
            this->computeNLBmatrixAt(gp, b, tStep);
            this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
            convertTangentToResultantTangent3d(ds, d, gp);

            FloatMatrix r(6, 6), rT(6, 6), dR(6, 6), rTDR(6, 6);
            computeCurrentGtoLStrainRotationMatrix(r, u, coordA, coordB, coordGp);

            rT.beTranspositionOf(r);
            dR.beProductOf(ds, r);
            rTDR.beProductOf(rT, dR);

            db.beProductOf(rTDR, b);
            db.times(1. / length);
            bt.beTranspositionOf(b);
            contrib.beProductOf(bt, db);
            answer.add(contrib);
        }

        return;
    }





    void
    Lattice3dNL::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    {
        //Compute normalised displacement jumps in global coordinate system first and then rotate it to the local coordiante system.

        FloatArray u;
        this->computeVectorOf(VM_Total, tStep, u);

        FloatArray coordA(3), coordB(3), coordGp(3), l1(3), l2(3), help(3);
        coordA = giveNode(1)->giveCoordinates();
        coordB = giveNode(2)->giveCoordinates();
        this->giveGpCoordinates(coordGp, gp);
        l1 = coordGp-coordA;
        l2 = coordB-coordGp;


        this->length = giveLength();

        FloatArray spin1(3), spin2(3);
        spin1.at(1) = u.at(4);  spin1.at(2) = u.at(5);  spin1.at(3) = u.at(6);
        spin2.at(1) = u.at(10); spin2.at(2) = u.at(11); spin2.at(3) = u.at(12);

        FloatMatrix R1, R2;
        this->computeGlobalRotationMatrix(R1, spin1);
        this->computeGlobalRotationMatrix(R2, spin2);

        FloatArray c1, c2;
        c1.beProductOf(R1, l1);
        c2.beProductOf(R2, l2);

        answer.resize(6);
        answer.at(1) = u.at(7) - u.at(1) - c1.at(1) - c2.at(1) + l1.at(1) + l2.at(1);
        answer.at(2) = u.at(8) - u.at(2) - c1.at(2) - c2.at(2) + l1.at(2) + l2.at(2);
        answer.at(3) = u.at(9) - u.at(3) - c1.at(3) - c2.at(3) + l1.at(3) + l2.at(3);

        // Rotational strain: exact relative rotation in the global frame, log(R2 R1^T),
        FloatMatrix R1t, Rrel;
        R1t.beTranspositionOf(R1);
        Rrel.beProductOf(R2, R1t);
        FloatArray theta;
        this->logRotationVector(Rrel, theta);
        answer.at(4) = theta.at(1);
        answer.at(5) = theta.at(2);
        answer.at(6) = theta.at(3);

        answer.times(1. / this->length);

        //Rotate strain vector to local coordinate system

        FloatMatrix rotMatrix(6, 6);
        computeCurrentGtoLStrainRotationMatrix(rotMatrix, u, coordA, coordB, coordGp);
        answer.rotatedWith(rotMatrix, 'n');
    }



    bool
    Lattice3dNL::computeGtoLStrainRotationMatrix(FloatMatrix &answer)
    {
        FloatMatrix lcs;
        answer.resize(6, 6);
        answer.zero();

        this->giveLocalCoordinateSystem(lcs);
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                answer.at(i, j) = lcs.at(i, j);
                answer.at(i + 3, j + 3) = lcs.at(i, j);
            }
        }

        return 1;
    }


    bool
    Lattice3dNL::computeGtoLRotationMatrix(FloatMatrix &answer)
    {
        return false;
    }


    void
    Lattice3dNL::giveInternalForcesVector(FloatArray &answer,
                                               TimeStep *tStep, int useUpdatedGpRecord)
    {
        FloatMatrix b, bt, bf, d;
        FloatArray u, stress, strain;
        this->computeVectorOf(VM_Total, tStep, u);
        this->length   = giveLength();

        FloatArray coordA(3), coordB(3);
        coordA = giveNode(1)->giveCoordinates();
        coordB = giveNode(2)->giveCoordinates();

        FloatMatrix rotationMatrixOne(3, 3), rotationMatrixTwo(3, 3);
        FloatArray spinOne(3), spinTwo(3);
        spinOne.at(1) = u.at(4);  spinOne.at(2) = u.at(5);  spinOne.at(3) = u.at(6);
        spinTwo.at(1) = u.at(10); spinTwo.at(2) = u.at(11); spinTwo.at(3) = u.at(12);
        this->computeGlobalRotationMatrix(rotationMatrixOne, spinOne);
        this->computeGlobalRotationMatrix(rotationMatrixTwo, spinTwo);

        answer.resize(12);
        answer.zero();

        IntegrationRule *iRule = this->integrationRulesArray [ 0 ].get();
        const int nGp = iRule->giveNumberOfIntegrationPoints();
        for ( int i = 0; i < nGp; ++i ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);

            FloatArray coordGp(3), l1(3), l2(3);
            this->giveGpCoordinates(coordGp, gp);
            l1 = coordGp - coordA;
            l2 = coordB - coordGp;

            if ( useUpdatedGpRecord == 1 ) {
                LatticeMaterialStatus *lmatStat = dynamic_cast < LatticeMaterialStatus * > ( this->giveCrossSection()->giveMaterial(gp)->giveStatus(gp) );
                stress = lmatStat->giveLatticeStress();
            } else {
                if ( !this->isActivated(tStep) ) {
                    strain.zero();
                }
                this->computeStrainVector(strain, gp, tStep);
                this->computeStressVector(stress, strain, gp, tStep);
            }

            FloatArray s;
            convertStressToResultants3d(s, stress, gp);

            FloatArray c1(3), c2(3);
            c1.beProductOf(rotationMatrixOne, l1);
            c2.beProductOf(rotationMatrixTwo, l2);

            FloatMatrix rotMatrix(6, 6);
            computeCurrentGtoLStrainRotationMatrix(rotMatrix, u, coordA, coordB, coordGp);

            //Calculate sectional forces in global coordinate system
            s.rotatedWith(rotMatrix, 't');

            answer.at(1)  += -s.at(1);
            answer.at(2)  += -s.at(2);
            answer.at(3)  += -s.at(3);
            answer.at(4)  +=  s.at(2) * c1.at(3) - s.at(3) * c1.at(2) - s.at(4);
            answer.at(5)  += -s.at(1) * c1.at(3) + s.at(3) * c1.at(1) - s.at(5);
            answer.at(6)  +=  s.at(1) * c1.at(2) - s.at(2) * c1.at(1) - s.at(6);
            answer.at(7)  +=  s.at(1);
            answer.at(8)  +=  s.at(2);
            answer.at(9)  +=  s.at(3);
            answer.at(10) +=  s.at(2) * c2.at(3) - s.at(3) * c2.at(2) + s.at(4);
            answer.at(11) += -s.at(1) * c2.at(3) + s.at(3) * c2.at(1) + s.at(5);
            answer.at(12) +=  s.at(1) * c2.at(2) - s.at(2) * c2.at(1) + s.at(6);
        }
    }


    void Lattice3dNL::computeCurrentGtoLStrainRotationMatrix(FloatMatrix &GtoLCurrent,
                                                                  const FloatArray &u,
                                                                  const FloatArray &coordA,
                                                                  const FloatArray &coordB,
                                                                  const FloatArray &coordGP)
    {

    // Position of GP along the element in the reference configuration
    FloatArray e = coordB - coordA;
    FloatArray r = coordGP - coordA;

    double num = r.at(1) * e.at(1) + r.at(2) * e.at(2) + r.at(3) * e.at(3);
    double den = e.at(1) * e.at(1) + e.at(2) * e.at(2) + e.at(3) * e.at(3);

    double t = num / den;
    if ( t < 0.0 ) {
        t = 0.0;
    }
    if ( t > 1.0 ) {
        t = 1.0;
    }

    // Interpolate nodal spin in global coordinates
    FloatArray spinOne(3), spinTwo(3), spinAverage(3), help(3);

    spinOne.at(1) = u.at(4);
    spinOne.at(2) = u.at(5);
    spinOne.at(3) = u.at(6);

    spinTwo.at(1) = u.at(10);
    spinTwo.at(2) = u.at(11);
    spinTwo.at(3) = u.at(12);

    spinAverage = spinOne;
    spinAverage.times(1.0 - t);

    help = spinTwo;
    help.times(t);

    spinAverage.add(help);

    // Initial global-to-local strain rotation matrix
    FloatMatrix GtoL0(6, 6);
    computeGtoLStrainRotationMatrix(GtoL0);

    // Extract initial local basis vectors in global coordinates
    FloatArray e10(3), e20(3), e30(3);
    for ( int j = 1; j <= 3; ++j ) {
        e10.at(j) = GtoL0.at(1, j);
        e20.at(j) = GtoL0.at(2, j);
        e30.at(j) = GtoL0.at(3, j);
    }

    // Rotate the initial basis by the average spin in GLOBAL coordinates
    FloatMatrix Qrot(3, 3);
    this->computeGlobalRotationMatrix(Qrot, spinAverage);

    FloatArray e1(3), e2(3), e3(3);
    e1.beProductOf(Qrot, e10);
    e2.beProductOf(Qrot, e20);
    e3.beProductOf(Qrot, e30);

    // Re-orthonormalise and enforce right-handedness
    e1.normalize();

    e3.beVectorProductOf(e1, e2);
    e3.normalize();

    e2.beVectorProductOf(e3, e1);
    e2.normalize();

    // Build current global-to-local 3x3 direction cosine matrix
    FloatMatrix Qloc(3, 3);
    Qloc.zero();
    for ( int j = 1; j <= 3; ++j ) {
        Qloc.at(1, j) = e1.at(j);
        Qloc.at(2, j) = e2.at(j);
        Qloc.at(3, j) = e3.at(j);
    }

    // Assemble the 6x6 strain/resultant rotation matrix
    GtoLCurrent.resize(6, 6);
    GtoLCurrent.zero();
    for ( int i = 1; i <= 3; ++i ) {
        for ( int j = 1; j <= 3; ++j ) {
            GtoLCurrent.at(i, j) = Qloc.at(i, j);
            GtoLCurrent.at(i + 3, j + 3) = Qloc.at(i, j);
        }
    }

    // Ensure local axis 1 points along the current element axis
    FloatArray x1_cur(3), x2_cur(3), vaxis(3);
    x1_cur = coordA;
    x2_cur = coordB;

    x1_cur.at(1) += u.at(1);
    x1_cur.at(2) += u.at(2);
    x1_cur.at(3) += u.at(3);

    x2_cur.at(1) += u.at(7);
    x2_cur.at(2) += u.at(8);
    x2_cur.at(3) += u.at(9);

    vaxis.beDifferenceOf(x2_cur, x1_cur);

    double norm = sqrt(vaxis.at(1) * vaxis.at(1)
                     + vaxis.at(2) * vaxis.at(2)
                     + vaxis.at(3) * vaxis.at(3));

    if ( norm > 0.0 ) {
        vaxis.times(1.0 / norm);

        FloatArray e1_from_R(3);
        e1_from_R.at(1) = GtoLCurrent.at(1, 1);
        e1_from_R.at(2) = GtoLCurrent.at(1, 2);
        e1_from_R.at(3) = GtoLCurrent.at(1, 3);

        double dot = e1_from_R.at(1) * vaxis.at(1)
                   + e1_from_R.at(2) * vaxis.at(2)
                   + e1_from_R.at(3) * vaxis.at(3);

        if ( dot < 0.0 ) {
            for ( int i = 1; i <= 3; ++i ) {
                for ( int j = 1; j <= 3; ++j ) {
                    GtoLCurrent.at(i, j) *= -1.0;
                    GtoLCurrent.at(i + 3, j + 3) *= -1.0;
                }
            }
        }
    }
    }
} // end namespace oofem
