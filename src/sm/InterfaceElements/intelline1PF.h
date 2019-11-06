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

#ifndef intelline1pf_h
#define intelline1pf_h

#include "sm/Elements/Interfaces/structuralinterfaceelement.h"
#include "sm/Elements/phasefieldelement.h"

#define _IFT_IntElLine1PF_Name "IntElLine1PF"
#define _IFT_IntElLine1PF_prescribedDamage "prescribeddamage"

namespace oofem {
class FEI2dLineLin;

/**
 * This class implements a two dimensional interface element.
 * Even if geometry approx is quadratic, the element is assumed straight
 * If not straight, the rotation matrix depends on actual integration point
 * and stiffness and strain computations should be modified.
 */
class IntElLine1PF : public StructuralInterfaceElement, PhaseFieldElement
{
protected:
    static FEI2dLineLin interp;

public:
    IntElLine1PF(int n, Domain *d);

    FEInterpolation *giveInterpolation() const override;

    int computeNumberOfDofs() override { return 12; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    double computeAreaAround(GaussPoint *gp) override;
    void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    FloatArrayF<2> computeCovarBaseVectorAt(GaussPoint *gp) const;

    int testElementExtension(ElementExtension ext) override { return 0; }

    //virtual Interface *giveInterface(InterfaceType) { return NULL; }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_IntElLine1PF_Name; }
    const char *giveClassName() const override { return "IntElLine1PF"; }
    void initializeFrom(InputRecord &ir) override;

    void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep) override
    {
        answer = this->giveInterfaceCrossSection()->giveEngTraction_2d(jump, gp, tStep);
    }

    void giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep) override
    {
        answer = this->giveInterfaceCrossSection()->give2dStiffnessMatrix_Eng(rMode, ip, tStep);
    }

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;


    // PF
    int giveNumberOfDofManagers() const override { return 4; }
    StructuralInterfaceElement *giveElement() override { return this; }
	void giveDofManDofIDMask_u(IntArray &answer) override;
    void giveDofManDofIDMask_d(IntArray &answer) override;

    virtual void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    virtual void computeStiffnessMatrix_ud(FloatMatrix &, MatResponseMode, TimeStep *);
    virtual void computeStiffnessMatrix_dd(FloatMatrix &, MatResponseMode, TimeStep *);
    //virtual void computeStiffnessMatrix_du(FloatMatrix &, MatResponseMode, TimeStep *);

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    virtual void giveInternalForcesVectorUD(FloatArray &fu, FloatArray &fd, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual double computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    virtual void computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer );
    virtual void computeBd_vectorAt(GaussPoint *gp, FloatArray &B);
    virtual void computeNd_vectorAt(const FloatArray &lCoords, FloatArray &N);
    virtual double computeFreeEnergy(GaussPoint *gp, TimeStep *tStep);

    double computeOldG(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double neg_MaCauley(double par) const
    {
        return 0.5 * ( std::abs(par) - par );
    }

    double neg_MaCauleyPrime(double par) const
    {
        return 0.5 * ( std::abs(par)/(par + 1.0e-12) - 1.0 ); // 0.5*(sign - 1) taken from Ragnars code
    }

    double MaCauley(double par) const
    {
        return 0.5 * ( std::abs(par) + par );
    }

    double MaCauleyPrime(double par) const
    {
        return 0.5 * ( std::abs(par)/(par + 1.0e-12) + 1.0 ); // 0.5*(sign + 1) taken from Ragnars code
    }


    void solveForLocalDamage(FloatMatrix &answer, TimeStep *tStep);
    void computeGMatrix(FloatMatrix &answer, const double damage, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);

protected:
    FloatArray unknownVectorU;
    FloatArray unknownVectorD;
    FloatArray deltaUnknownVectorD;

    FloatArray deltaAlpha; 
    FloatArray alpha;
    FloatArray oldAlpha;
    double prescribed_damage = 0.;

    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeGaussPoints() override;

    virtual int giveApproxOrder() { return 1; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_quad_1_interface; };
};
} // end namespace oofem
#endif
