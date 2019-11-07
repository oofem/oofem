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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef latticeneumanncouplingnode_h
#define latticeneumanncouplingnode_h

#include "node.h"
#include "domain.h"
#include "floatarray.h"


///@name Input fields for LatticeCouplingNode
//@{
#define _IFT_LatticeNeumannCouplingNode_Name "latticeneumanncouplingnode"
#define _IFT_LatticeNeumannCouplingNode_direction "direction"
#define _IFT_LatticeNeumannCouplingNode_couplingnodes "couplingnodes"
//@}

namespace oofem {
class Node;
class Dof;
class NodalLoad;
class TimeStep;
class FloatArray;
class IntArray;

/**
 * Class implementing lattice coupling node
 * This DOF manager is based on the standard node but allows to read
 * in associated DOFs and elements which are used to calculated coupling forces
 */

class LatticeNeumannCouplingNode : public Node
{
protected:
    /// Array storing nodal coordinates.
    FloatArray directionVector;
    IntArray couplingNodes;

public:

    LatticeNeumannCouplingNode(int n, Domain *aDomain);                       // constructor

    ~LatticeNeumannCouplingNode();                                            // destructor

    virtual void  computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode);

    const char *giveClassName() const override { return "LatticeNeumannCouplingNode"; }
    void initializeFrom(InputRecord &ir) override;

    IntArray *giveCouplingNodes();

    void computeLoadCouplingContribution(FloatArray &answer, TimeStep *stepN);

    void printYourself() override;
};
} // end namespace oofem

#endif // latticecouplingnode_h
