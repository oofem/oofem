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
#if 0


#ifndef phasefieldelement_h
#define phasefieldelement_h

#include "coupledfieldselement.h"
#include "fei2dquadlin.h"

namespace oofem {
/**
 * Abstract class for phase field formulation
 */
class PhaseFieldElement : public CoupledFieldsElement
{
protected:
    IntArray loc_u, loc_d;
    int nlGeo;

public:
    PhaseFieldElement(int i, Domain *aDomain) : CoupledFieldsElement(i, aDomain) {};
    virtual ~PhaseFieldElement() {}

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask_u(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_d(IntArray &answer) = 0;

    // definition & identification
    virtual const char *giveClassName() const { return "PhaseFieldElement"; }

protected:

    virtual void computeNStressAt(GaussPoint *, FloatArray &) = 0;
    virtual void computeBStressAt(GaussPoint *, FloatArray &) = 0;
    virtual double computeVolumeAround(GaussPoint *) = 0;

    void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeStiffnessMatrixGen(FloatMatrix &, MatResponseMode, TimeStep *, 
        void (*Nfunc)(GaussPoint*, FloatMatrix), void (*Bfunc)(GaussPoint*, FloatMatrix),
        void (*NStiffness)(GaussPoint*, FloatMatrix), void (*BStiffness)(GaussPoint*, FloatMatrix),
        double (*volumeAround)(GaussPoint*)
        );

    void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ud(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_dd(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_du(FloatMatrix &, MatResponseMode, TimeStep *);

    void Duu_B(FloatMatrix &, MatResponseMode, GaussPoint *, TimeStep *);
    void Duu_d(FloatMatrix &, MatResponseMode, GaussPoint *, TimeStep *);
    double computeG(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeGPrim(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);

    void computeBStress_u(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord);
    void computeNStress_d(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord);
    void computeBStress_d(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord);
    

    void computeDisplacementUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN);
    void computeDamageUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN);
    
    //Interpolation matrices
    virtual void computeBd_matrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    //virtual void computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
    virtual void computeNd_matrixAt(const FloatArray &lCoords, FloatMatrix &N);

    static FEI2dQuadLin interpolation_u;
    static FEI2dQuadLin interpolation_d;
};
} // end namespace oofem

#endif
#endif