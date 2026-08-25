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

#ifndef continuumframenode_h
#define continuumframenode_h

#include "hangingnode.h"

///@name Input fields for ContinuumFrameNode
//@{
#define _IFT_ContinuumFrameNode_Name "continuumframenode"
//@}

namespace oofem {
/**
 * Class implementing a node with rotational DOFs embedded in a continuum (solid) mesh.
 *
 * For its translational DOFs the node behaves exactly like a HangingNode.
 * However, rotational DOFs are constrained to the infinitesimal
 * continuum rotations.
 *
 * Author: Peter Grassl
 */
class OOFEM_EXPORT ContinuumFrameNode : public HangingNode
{
public:
    /**
     * Constructor. Creates a continuum frame node with number n, belonging to aDomain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    ContinuumFrameNode(int n, Domain * aDomain);
    /// Destructor.
    virtual ~ContinuumFrameNode(void) { }

    void postInitialize() override;

    const char *giveClassName() const override { return "ContinuumFrameNode"; }
    const char *giveInputRecordName() const override { return _IFT_ContinuumFrameNode_Name; }
};
} // end namespace oofem
#endif // continuumframenode_h
