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

#ifndef qcnode_h
#define qcnode_h

#include "node.h"

///@name Input fields for HangingNode
//@{
#define _IFT_qcNode_Name "qcnode"
//@}

namespace oofem {

class ParamKey;

/**
 * Class implementing hanging node connected to other nodes (masters) using interpolation.
 * Hanging node possess no degrees of freedom - all values are interpolated from corresponding master dofs.
 *
 * The introduction of hanging nodes allows, for example, to include reinforcing bar elements inside
 * arbitrary FE-mesh of concrete specimen or facilitates the local refinement of FE-mesh.
 *
 * The contributions of hanging node are localized directly to master related equations.
 * The hanging node can not have its own boundary or initial conditions,
 * they are determined completely from master dof conditions.
 * The local coordinate system in slave is not supported in current implementation, the global cs applies.
 * On the other hand, hanging node can be loaded independently of master.
 *
 * If no master element number is supplied (or negative) then it will locate it using the global coordinates.
 * and if no master region number is supplied (or zero), it will look for elements in all regions.
 */
class OOFEM_EXPORT qcNode : public Node
{
protected:
    /// Number of the master element.
    int masterElement=-1;
    /// Region of the master element (used for automatic detection).
    int masterRegion=0;
    /// Type of qcNode (0 deactive, 1 master, 2 hanging)
    int qcNodeTypeLabel;
#ifdef __OOFEG
    /// Flag whether node is fully initialized already.
    bool initialized;
#endif

    /// Static input field keys
    static ParamKey IPK_qcNode_masterElement;
    static ParamKey IPK_qcNode_masterRegion;

public:
    /**
     * Constructor. Creates a hanging node with number n, belonging to aDomain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    qcNode(int n, Domain * aDomain);
    /// Destructor.
    virtual ~qcNode(void) { }

    void initializeFrom(InputRecord &ir, int priority) override;
    void postInitialize() override;
    void postInitializeAsHangingNode();
    int checkConsistency() override;
    bool isDofTypeCompatible(dofType type) const override { return ( type == DT_master || type == DT_slave ); }

    virtual bool initializeAsRepnode();
    virtual void setAsRepnode();
    virtual void setAsHanging();
    int giveQcNodeType() override { return this->qcNodeTypeLabel; }
    virtual int giveMasterElementNumber() { return this->masterElement; }

    void printOutputAt(FILE *stream, TimeStep *tStep) override;

    const char *giveClassName() const override { return "qcNode"; }
    const char *giveInputRecordName() const override { return _IFT_qcNode_Name; }
};
} // end namespace oofem
#endif // qcnode_h
