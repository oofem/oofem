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

#ifndef latticedirichletcouplingnode_h
#define latticedirichletcouplingnode_h

#include "node.h"
#include "domain.h"
#include "floatarray.h"


///@name Input fields for LatticeCouplingNode
//@{
#define _IFT_LatticeDirichletCouplingNode_Name "latticedirichletcouplingnode"
#define _IFT_LatticeDirichletCouplingNode_couplingelements "couplingelements"
//@}

namespace oofem {
class Node;
class Dof;
class NodalLoad;
class TimeStep;
class IntArray;
class ParamKey;

/**
 * Class implementing lattice coupling node
 * This DOF manager is based on the standard node but allows to read
 * in associated DOFs and elements which are used to calculated coupling forces
 */

class LatticeDirichletCouplingNode : public Node
{
protected:
    /// Array storing nodal coordinates.
    IntArray couplingElements; //IntArray used here cause the element numbers are all whole numbers (integers) instead of using a float array .
    int couplingType;
    double activedirichletbc;


    static ParamKey IPK_LatticeDirichletCouplingNode_couplingelements;
public:

    LatticeDirichletCouplingNode(int n, Domain *aDomain);                       // constructor

    ~LatticeDirichletCouplingNode();                                            // destructor

    const char *giveClassName() const override { return "LatticeDirichletCouplingNode"; }
    void initializeFrom(InputRecord &ir, int priority) override;
    void postInitialize() override;

    IntArray *giveCouplingNodes();

    IntArray *giveCouplingElements();

    virtual void giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *stepN, bool padding = false);

    double giveUnknown(ValueModeType mode, TimeStep *stepN);

    double computeUnknownCouplingContribution(TimeStep *stepN);

    void printOutputAt(FILE *stream, TimeStep *stepN) override;

    void printYourself() override;
};
} // end namespace oofem

#endif // latticedirichletcouplingnode_h
