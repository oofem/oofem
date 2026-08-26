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

#ifndef latticeneumanncoupling_h
#define latticeneumanncoupling_h

#include "activebc.h"
#include "floatarray.h"
#include "intarray.h"

///@name Input fields for LatticeNeumannCoupling
//@{
#define _IFT_LatticeNeumannCoupling_Name "latticeneumanncoupling"
#define _IFT_LatticeNeumannCoupling_smnodes "smnodes"
#define _IFT_LatticeNeumannCoupling_tmnodes "tmnodes"
#define _IFT_LatticeNeumannCoupling_direction "direction"
//@}

namespace oofem {
class Node;
class TimeStep;

/**
 * Neumann-type hydro-mechanical coupling boundary condition (Approach 1 of
 * Grassl, Fahy, Gallipoli, Wheeler 2015). Reads the pore-fluid pressure P_f
 * from coupled transport lattice nodes (in the transport slave problem of a
 * StaggeredProblem) and applies it as an external force on mechanical lattice
 * nodes: f = P_f * distance * directionVector. Replaces the old
 * LatticeNeumannCouplingNode, whose DofManager-level load hook no longer
 * exists; the contribution is now injected through assembleVector.
 */
class LatticeNeumannCoupling : public ActiveBoundaryCondition
{
protected:
    /// Mechanical (this domain) nodes receiving the coupling force.
    IntArray smNodes;
    /// Transport (slave problem) nodes whose P_f drives each smNode.
    IntArray tmNodes;
    /// Load direction (size matches the mechanical node DOFs: D_u, D_v, R_w).
    FloatArray directionVector;

public:
    LatticeNeumannCoupling(int n, Domain *d) : ActiveBoundaryCondition(n, d) { }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorms = nullptr, void *lock = nullptr) override;

    const char *giveClassName() const override { return "LatticeNeumannCoupling"; }
    const char *giveInputRecordName() const override { return _IFT_LatticeNeumannCoupling_Name; }

protected:
    /// Returns the transport node (in the coupled transport slave problem) for the given id.
    Node *giveCoupledTransportNode(int tmNodeId);
};
} // namespace oofem
#endif // latticeneumanncoupling_h
