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
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser Base Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser Base Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Base Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef basemixedpressureelement_h
#define basemixedpressureelement_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/nlstructuralelement.h"

namespace oofem {
/**
 * Base class for mixed u-p formulation.
 * @author Martin Horak
 */
class BaseMixedPressureElement
{
protected:
    IntArray displacementDofsOrdering, pressureDofsOrdering;
    IntArray locationArray_u, locationArray_p;


public:
    BaseMixedPressureElement();
    virtual ~BaseMixedPressureElement() { }

    virtual void initializeFrom(InputRecord &ir);

protected:

    /// Pure virtual functions
    virtual NLStructuralElement *giveElement() = 0;
    virtual void computeVolumetricBmatrixAt(GaussPoint *gp, FloatArray &Bvol, NLStructuralElement *element) = 0;

    virtual void computePressureNMatrixAt(GaussPoint *, FloatArray &) = 0;

    virtual int giveNumberOfPressureDofs() = 0;
    virtual int giveNumberOfDisplacementDofs() = 0;
    virtual int giveNumberOfDofs() = 0;

    virtual void giveDofManDofIDMask_u(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_p(IntArray &answer) = 0;
    /// End of pure virtual functions



    virtual void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_up(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_pp(FloatMatrix &, MatResponseMode, TimeStep *);




    void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    void computePressure(double &answer,  GaussPoint *gp, TimeStep *tStep);


    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveInternalForcesVector_p(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);

    void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    void computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);

    virtual IntArray &giveDisplacementDofsOrdering() { return displacementDofsOrdering; }
    virtual IntArray &giveMicromorphicDofsOrdering() { return pressureDofsOrdering; }
    // void giveLocationArrayOfDofIDs( IntArray &answer, const UnknownNumberingScheme &s, const IntArray &dofIdArray );
    void giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_p, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u, const IntArray &dofIdArray_p);
    virtual void postInitialize();
    virtual void updateInternalState(TimeStep *tStep);
};
} // end namespace oofem

#endif
