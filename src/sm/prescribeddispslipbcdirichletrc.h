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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef prescribeddispslipbcdirichletrc_h
#define prescribeddispslipbcdirichletrc_h

#include "prescribeddispsliphomogenization.h"
#include "boundarycondition.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for PrescribedDispSlipBCDirichletRC
//@{
#define _IFT_PrescribedDispSlipBCDirichletRC_Name "prescribeddispslipbcdirichletrc"
#define _IFT_PrescribedDispSlipBCDirichletRC_ConcreteBoundary "conboundset"
#define _IFT_PrescribedDispSlipBCDirichletRC_ReinfXBound "reinfxbound"
#define _IFT_PrescribedDispSlipBCDirichletRC_ReinfYBound "reinfybound"
//@}

namespace oofem {
/**
 * Prescribes macroscopic displacement gradient, reinforcement slip field and slip gradient
 * on the reinforced concrete RVE using Dirichlet boundary conditions, cf.
 * Sciegaj, A., Larsson, F., Lundgren, K., Nilenius, F., & Runesson, K. (2019). A multiscale model for reinforced concrete with macroscopic variation of reinforcement slip. Computational Mechanics, 63(2), 139–158. https://doi.org/10.1007/s00466-018-1588-3
 * Macroscopic slip and slip gradient fields are optional (see PrescribedDispSlipHomogenization).
 * If not specified, this BC will prescribe only the displacement gradient, i.e., it will work the same was as a PrescribedGradient bc.
 * Works with 2D RVEs comprising solid elements (concrete), reinforcement (beam/truss elements) and interface elements in between.
 * Used in multiscale analyses of reinforced concrete structures. Currently, only orthogonal reinforcement in X and Y direction is supported.
 *
 * This BC can be applied to the set of:
 *      - nodes lying at the concrete boundary. In order to prescribe macroscopic displacement gradient on both concrete and steel,
 *          the node set this BC is applied to must include the end (boundary) nodes of the reinforcement bars as well. conboundset should contain only elementboundaries of concrete elements, e.g.
 *          PrescribedDispSlipBCDirichletRC 1 loadTimeFunction 1 dofs 2 1 2 ccoord 3 0.0 0.0 0.0 dispGrad 2 2 {0 0; 0 0} set 5 conboundset 6
 *          where set 5 contains all boundary nodes (concrete+steel) and set 6 contains boundary of concrete elements
 *      - elementboundaries of concrete elements. In this case displacement gradient is not prescribed on the reinforcement.
 *          PrescribedDispSlipBCDirichletRC 1 loadTimeFunction 1 dofs 2 1 2 ccoord 3 0.0 0.0 0.0 dispGrad 2 2 {0 0; 0 0} set 5 conboundset 5
 *          where set 5 contains boundary of concrete elements
 *
 * If the macroscopic slip/slip gradient field is to be prescribed onto the RVE, the optional node sets reinfxbound and reinfybound must be specified.
 * These node sets contain the boundary nodes of the horizontal and vertical reinforcement, respectively.
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedDispSlipBCDirichletRC : public BoundaryCondition, public PrescribedDispSlipHomogenization
{
public:
    PrescribedDispSlipBCDirichletRC(int n, Domain *d) : BoundaryCondition(n, d) { }
    virtual ~PrescribedDispSlipBCDirichletRC() { }

    double give(Dof *dof, ValueModeType mode, double time) override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    virtual void updateCoefficientMatrix(FloatMatrix &C);
    double domainSize(Domain *d, int set) override;

    void computeStress(FloatArray &sigma, TimeStep *tStep) override;
    void computeTransferStress(FloatArray &bStress, TimeStep *tStep) override;
    void computeReinfStress(FloatArray &rStress, TimeStep *tStep) override;

    void computeTangent(FloatMatrix &tangent, TimeStep *tStep) override;

    void scale(double s) override;

    const char *giveClassName() const override { return "PrescribedDispSlipBCDirichletRC"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedDispSlipBCDirichletRC_Name; }

protected:
    bool dispGradON = false; /// on/off flag specifying whether the displacement gradient should be applied via Neumann BCs
    bool slipON = false; /// on/off flag specifying whether the slip field should be applied via Neumann BCs
    bool slipGradON = false; /// on/off flag specifying whether the slip gradient should be applied via Neumann BCs

    int conBoundSet; //element boundaries set for the concrete solid
    int reinfXBound=0; //set containing end (boundary) nodes of horizontal rebars
    int reinfYBound=0; //set containing end (boundary) nodes of vertical rebars

    double giveOnConcrete(Dof *dof, int pos, const FloatArray u);
    double giveOnSteel(Dof *dof, int pos, const FloatArray u, const FloatArray us);
};
} // end namespace oofem

#endif // prescribeddispslipbcdirichletrc_h
