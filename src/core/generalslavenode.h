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
 *               Copyright (C) 1993 - 2022   Borek Patzak
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

#ifndef generalslavenode_h
#define generalslavenode_h

#include "node.h"
#include <vector>
///@name Input fields for SlaveNode
//@{
#define _IFT_GeneralSlaveNode_Name "generalslavenode"
#define _IFT_GeneralSlaveNode_masterSizes "mastersizes"
#define _IFT_GeneralSlaveNode_masterWeights "masterweights"
#define _IFT_GeneralSlaveNode_masterList "masterlist"
//@}

namespace oofem {
/**
 * Class implementing slave node connected to other nodes (masters) using predetermined weights.
 * The node can not have its own boundary or initial conditions,
 * they are determined completely from master dof conditions except for dofs of master type.
 * This approach is more general than SlaveNode because it enables individual slave node's dofs
 * to depend on different masters with different sizes
 * The GeneralSlaveNode is suitable, e.g., for setting periodic boundary conditions
 * @author Martin Horák
 */
class OOFEM_EXPORT GeneralSlaveNode : public Node
{
protected:

    /// Master nodes for all dofs.
    std::vector< IntArray >dofs_masterList;
    std::vector< IntArray >dofs_dofsList;
    std::vector< FloatArray >dofs_weightsList;

    IntArray masterSizes;

public:
    /**
     * Constructor. Creates a general slave  node with number n, belonging to aDomain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    GeneralSlaveNode(int n, Domain *aDomain) : Node(n, aDomain) {}
    /// Destructor.
    virtual ~GeneralSlaveNode(void) {}

    void initializeFrom(InputRecord &ir) override;
    virtual void postInitialize() override;
    virtual bool isDofTypeCompatible(dofType type) const override {
        return ( type == DT_master || type == DT_slave );
    }
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f) override;
    virtual const char *giveClassName() const override {
        return "GeneralSlaveNode";
    }
    virtual const char *giveInputRecordName() const override {
        return _IFT_GeneralSlaveNode_Name;
    }
};
} // end namespace oofem
#endif // generalslavenode_h
