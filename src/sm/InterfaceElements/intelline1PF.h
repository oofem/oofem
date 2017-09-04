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

#include "Elements/Interfaces/structuralinterfaceelement.h"
#include "Elements/phasefieldelement.h"
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
    virtual ~IntElLine1PF() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual int computeNumberOfDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual double computeAreaAround(GaussPoint *gp);
    virtual void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeCovarBaseVectorAt(GaussPoint *gp, FloatArray &G);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    //virtual Interface *giveInterface(InterfaceType) { return NULL; }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_IntElLine1PF_Name; }
    virtual const char *giveClassName() const { return "IntElLine1PF"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
    {
        this->giveInterfaceCrossSection()->giveEngTraction_2d(answer, gp, jump, tStep);
    }

    virtual void giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep)
    {
        this->giveInterfaceCrossSection()->give2dStiffnessMatrix_Eng(answer, rMode, ip, tStep);
    }

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);



    // PF
    int giveNumberOfDofManagers() const { return 4; }
    StructuralInterfaceElement *giveElement() { return this; }
	virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_d(IntArray &answer);

    virtual void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    virtual void computeStiffnessMatrix_ud(FloatMatrix &, MatResponseMode, TimeStep *);
    virtual void computeStiffnessMatrix_dd(FloatMatrix &, MatResponseMode, TimeStep *);
    //virtual void computeStiffnessMatrix_du(FloatMatrix &, MatResponseMode, TimeStep *);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void giveInternalForcesVectorUD(FloatArray &fu, FloatArray &fd, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual double computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    virtual void computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer );
    virtual void computeBd_vectorAt(GaussPoint *gp, FloatArray &B);
    virtual void computeNd_vectorAt(const FloatArray &lCoords, FloatArray &N);
    virtual double computeFreeEnergy(GaussPoint *gp, TimeStep *tStep);
    
    double computeOldG(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double neg_MaCauley(double par)
    {
        return 0.5 * ( abs(par) - par );
    }

    double neg_MaCauleyPrime(double par)
    {
        return 0.5 * ( abs(par)/(par + 1.0e-12) - 1.0 ); // 0.5*(sign - 1) taken from Ragnars code
    }

    double MaCauley(double par)
    {
        return 0.5 * ( abs(par) + par );
    }

    double MaCauleyPrime(double par)
    {
        return 0.5 * ( abs(par)/(par + 1.0e-12) + 1.0 ); // 0.5*(sign + 1) taken from Ragnars code
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
    double prescribed_damage;

    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    virtual int giveApproxOrder() { return 1; }
    Element_Geometry_Type giveGeometryType() const { return EGT_quad_1_interface; };
};
} // end namespace oofem
#endif
