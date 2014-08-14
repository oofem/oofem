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
#if 1
#ifndef coupledfieldselement_h
#define coupledfieldselement_h

#include "../sm/Elements/nlstructuralelement.h"

namespace oofem {
/**
 * Abstract class for gradient formulation of coupled damage-plasticity model(GradDp).
 * Yield function is formulated in the effective stress space and damage is driven by the nonlocal(over-nonlocal) cumulated plastic strain.
 * The new nonlocal degrees of freedom (with the meaning of the nonlocal cumulated plastic strain) are introduced with lower order of approximation functions than the displacement field to avoid spurious stress oscillations
 */
class CoupledFieldsElement : public NLStructuralElement
{
protected:
    int nlGeo;

public:
    CoupledFieldsElement(int i, Domain *aDomain);
    virtual ~CoupledFieldsElement() {}

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask (int inode, EquationID ut, IntArray& answer) const = 0;

protected:
    
    virtual void computeNStressAt(GaussPoint *, FloatArray &) = 0;
    virtual void computeBStressAt(GaussPoint *, FloatArray &) = 0;

    virtual double computeVolumeAround(GaussPoint *) = 0;
    void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *) = 0;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) = 0;



    void computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType valueMode, TimeStep *stepN, FloatArray &answer);
    void computeLocationArrayOfDofIDs(const IntArray &dofIdArray, IntArray &answer);

    void computeStiffnessMatrixGen(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep, 
        void (*Nfunc)(GaussPoint*, FloatMatrix), 
        void (*Bfunc)(GaussPoint*, FloatMatrix),
        void (*NStiffness)(FloatMatrix, MatResponseMode, GaussPoint*, TimeStep*), 
        void (*BStiffness)(FloatMatrix, MatResponseMode, GaussPoint*, TimeStep*),
        double (*volumeAround)(GaussPoint*) );

    
    void giveInternalForcesVectorGen(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord, 
        void (*Nfunc)(GaussPoint*, FloatMatrix), void (*Bfunc)(GaussPoint*, FloatMatrix, int, int), //(GaussPoint*, FloatMatrix)
        void (*NStress)(GaussPoint*, FloatArray), void (*BStress)(GaussPoint*, FloatArray),
        double (*volumeAround)(GaussPoint*)
        );

    //void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);

    
};
} // end namespace oofem

#endif
#endif