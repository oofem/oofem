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

#ifndef q4axisymm_h
#define q4axisymm_h

#include "Elements/structuralelement.h"
#include "Elements/structural2delement.h"
#include "zznodalrecoverymodel.h"

///@name Input fields for Q4Axisymm
//@{
#define _IFT_Q4Axisymm_Name "q4axisymm"
#define _IFT_Q4Axisymm_nipfish "nipfish"
//@}

namespace oofem {
class FEI2dQuadQuad;

/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral -
 * elasticity finite element for axisymmetric 3d continuum.
 * Each node has 2 degrees of freedom.
 */
class Q4Axisymm : public AxisymElement, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dQuadQuad interp;
    int numberOfFiAndShGaussPoints;

public:
    Q4Axisymm(int n, Domain * d);
    virtual ~Q4Axisymm();

    virtual FEInterpolation *giveInterpolation() const;

    // definition & identification
    virtual Interface *giveInterface(InterfaceType);
    virtual const char *giveInputRecordName() const { return _IFT_Q4Axisymm_Name; }
    virtual const char *giveClassName() const { return "Q4axisymm"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    
protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    
    
#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type);
#endif
};
} // end namespace oofem
#endif // q4axisymm_h
