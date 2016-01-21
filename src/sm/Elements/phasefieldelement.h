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

#ifndef phasefieldelement_h
#define phasefieldelement_h

#include "../sm/Elements/PlaneStress/qplanstrss.h"

namespace oofem {
/**
 * Abstract class for phase field formulation
 */
class PhaseFieldElement
{
protected:
    IntArray loc_u, loc_d;

public:
    PhaseFieldElement( int i, Domain *aDomain );
    virtual ~PhaseFieldElement() {}

    virtual NLStructuralElement *giveElement( ) = 0;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask_u(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_d(IntArray &answer) = 0;

    virtual const char *giveClassName() const { return "PhaseFieldElement"; }

    int computeNumberOfDofs();
    void computeLocationArrayOfDofIDs(const IntArray &dofIdArray, IntArray &answer);

    double computeFreeEnergy( GaussPoint *gp, TimeStep *tStep );

    ///@todo move these to a cross section model later
    double internalLength;
    double giveInternalLength( ) { return internalLength; }
    double criticalEnergy;
    double giveCriticalEnergy() { return criticalEnergy; }
    double relaxationTime;
    double giveRelaxationTime( ) { return relaxationTime; }

protected:

    virtual void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ud(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_dd(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_du(FloatMatrix &, MatResponseMode, TimeStep *);

    
    double computeG(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeGPrim(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void computeBStress_u(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord);
    void computeNStress_d( FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord );

    void computeDisplacementUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN);
    void computeDamageUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN);
    
    //Interpolation matrices
    virtual void computeBd_matrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeNd_matrixAt(const FloatArray &lCoords, FloatMatrix &N);

};
} // end namespace oofem

#endif