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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#ifndef prescribedgradientdd_h
#define prescribedgradientdd_h

#include "prescribedgradient.h"
#include "floatarrayf.h"

///@name Input fields for PrescribedGradientDD
//@{
#define _IFT_PrescribedGradientDD_Name "prescribedgradientdd"
#define _IFT_PrescribedGradientDD_ConcreteBoundary "conboundset"
#define _IFT_PrescribedGradientDD_ReinfXBound "reinfxbound"
#define _IFT_PrescribedGradientDD_ReinfYBound "reinfybound"
//@}

namespace oofem {
/**
 * Prescribes a displacement gradient with Dirichlet-Dirichlet boundary condition, cf.
 * Sciegaj, A., Larsson, F., Lundgren, K., Nilenius, F., & Runesson, K. (2018).
 * Two‐scale finite element modelling of reinforced concrete structures: Effective response and subscale fracture development.
 * International Journal for Numerical Methods in Engineering, 114(10), 1074–1102. https://doi.org/10.1002/nme.5776
 * Works with 2D RVEs comprising solid elements (concrete), reinforcement (beam elements) and interface elements in between.
 * Used in multiscale analyses of reinforced concrete structures. Currently, only orthogonal reinforcement in X and Y direction is supported.
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedGradientDD : public PrescribedGradient
{
public:

    PrescribedGradientDD(int n, Domain *d) : PrescribedGradient(n, d) { }
    virtual ~PrescribedGradientDD() { }

    double give(Dof *dof, ValueModeType mode, double time) override;

    bcType giveType() const override { return DirichletBT; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void updateCoefficientMatrix(FloatMatrix &C) override;

    void scale(double s) override { mGradient.times(s); }
    double domainSize(Domain *d, int set) override;

    const char *giveClassName() const override { return "PrescribedGradient"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedGradient_Name; }

protected:
    int conBoundSet; //element boundaries set for the concrete solid
    int reinfXBound; //set containing end (boundary) nodes of horizontal rebars
    int reinfYBound; //set containing end (boundary) nodes of vertical rebars
};
} // end namespace oofem

#endif // prescribedgradientdd_h
