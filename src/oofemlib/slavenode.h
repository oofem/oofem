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

#ifndef slavenode_h
#define slavenode_h

#include "node.h"

///@name Input fields for SlaveNode
//@{
#define _IFT_SlaveNode_Name "slavenode"
#define _IFT_SlaveNode_masterDofManagers "masterdofman"
#define _IFT_SlaveNode_weights "weights"
//@}

namespace oofem {
/**
 * Class implementing slave node connected to other nodes (masters) using predetermined weights.
 * Hanging node possess no degrees of freedom - all values are interpolated from corresponding master dofs.
 *
 * The contributions of hanging node are localized directly to master related equations.
 * The node can not have its own boundary or initial conditions,
 * they are determined completely from master dof conditions except for dofs of master type.
 * @see{HangingNode}
 */
class OOFEM_EXPORT SlaveNode : public Node
{
protected:
    /// Master nodes for all dofs.
    IntArray masterDofManagers;
    /// Weights for each master node.
    FloatArray masterWeights;

public:
    /**
     * Constructor. Creates a hanging node with number n, belonging to aDomain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    SlaveNode(int n, Domain * aDomain) : Node(n, aDomain) { }
    /// Destructor.
    virtual ~SlaveNode(void) { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void postInitialize();
    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_slave ); }
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    virtual const char *giveClassName() const { return "SlaveNode"; }
    virtual const char *giveInputRecordName() const { return _IFT_SlaveNode_Name; }
};
} // end namespace oofem
#endif // slavenode_h
