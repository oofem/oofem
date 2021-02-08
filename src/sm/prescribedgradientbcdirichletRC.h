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

#ifndef prescribedgradientbcdirichlet_h
#define prescribedgradientbcdirichlet_h

#include "prescribedgradient.h"
#include "floatarrayf.h"

///@name Input fields for PrescribedGradientBCDirichletRC
//@{
#define _IFT_PrescribedGradientBCDirichletRC_Name "prescribedgradientbcdirichletrc"
#define _IFT_PrescribedGradientBCDirichletRC_ConcreteBoundary "conboundset"
#define _IFT_PrescribedGradientBCDirichletRC_ReinfXBound "reinfxbound"
#define _IFT_PrescribedGradientBCDirichletRC_ReinfYBound "reinfybound"
//@}

namespace oofem {
/**
 * Prescribes a displacement gradient with Dirichlet boundary condition, cf.
 * Sciegaj, A., Larsson, F., Lundgren, K., Nilenius, F., & Runesson, K. (2018).
 * Two‐scale finite element modelling of reinforced concrete structures: Effective response and subscale fracture development.
 * International Journal for Numerical Methods in Engineering, 114(10), 1074–1102. https://doi.org/10.1002/nme.5776
 * Works with 2D RVEs comprising solid elements (concrete), reinforcement (beam elements) and interface elements in between.
 * Used in multiscale analyses of reinforced concrete structures. Currently, only orthogonal reinforcement in X and Y direction is supported.
 *
 * Current implementation supports both Dirichlet-Dirichlet and Dirichlet-Neumann BCs.
 * This BC is applied to the set of nodes lying at the concrete boundary (either node set or elementboundaries set will work).
 * If the optional reinfxbound and reinfybound sets are specified, the gradient is prescribed on both concerete and steel
 * (Dirichlet-Dirichlet, see example RVE input in sm/deepbeamfe2_01.in.rve).
 * If only the conboundset input record is specified, the gradient is prescribed only on concrete
 * (Dirichlet-Neumann, see example RVE input in sm/deepbeam2_02.in.rve).
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedGradientBCDirichletRC : public PrescribedGradient
{
public:

    PrescribedGradientBCDirichletRC(int n, Domain *d) : PrescribedGradient(n, d) { }
    virtual ~PrescribedGradientBCDirichletRC() { }

    double give(Dof *dof, ValueModeType mode, double time) override;

    bcType giveType() const override { return DirichletBT; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void updateCoefficientMatrix(FloatMatrix &C) override;

    void scale(double s) override { mGradient.times(s); }
    double domainSize(Domain *d, int set) override;

    const char *giveClassName() const override { return "PrescribedGradientBCDirichletRC"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedGradientBCDirichletRC_Name; }

protected:
    int conBoundSet; //element boundaries set for the concrete solid
    int reinfXBound=0; //set containing end (boundary) nodes of horizontal rebars
    int reinfYBound=0; //set containing end (boundary) nodes of vertical rebars
};
} // end namespace oofem

#endif // prescribedgradientbcdirichlet_h
