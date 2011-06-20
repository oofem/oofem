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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef hangingnode_h
#define hangingnode_h

#include "node.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

namespace oofem {
class FloatArray;
class IntArray;

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
 * @todo{Move all the grunt work to the FEI classes.}
 *
 * @note{When using slavemask it is necessary to be careful to avoid any discontinuities.}
 * @note{Nodal loads also require some thought.}
 */
class HangingNode : public Node
{
protected:
    /**
     * Type of interpolation from hangingnodes = 100*(number of master nodes) + 10*(order of polynomial approximation) + dimension
     * 211 - linear truss         321 - quadratic truss
     * 312 - linear triangle      622 - quadratic triangle
     * 412 - linear rectangle     822 - quadratic rectangle
     * 413 - linear tetrahedron  1023 - quadratic tetrahedron
     * 813 - linear hexahedron   2023 - quadratic hexahedron
     */
    int typeOfContrib;
    /// Local coordinates.
    FloatArray locoords;
    /// Count of Master Nodes.
    int countOfMasterNodes;
    /// Array of numbers of master DofManagers.
    IntArray masterDofMngr;
    /// Array of pointers to master Nodes.
    Node **masterNode;
    /// Vector of master contribution coefficients - SUMA of contributions == 1.0 ; if node has LCS, contributions are specified in this local coordinate system.
    FloatArray masterContribution;
#ifdef __OOFEG
    /// Flag whether consistency check already completed.
    bool consistencyChecked;
#endif

private:
    void deallocAuxArrays(void);
    /// Get natural coordinates from the hanging node coordinates.
    int computeNaturalCoordinates(void);

public:
    /**
     * Constructor. Creates a hanging node with number n, belonging to aDomain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    HangingNode(int n, Domain *aDomain);
    /// Destructor.
    ~HangingNode(void) { }

    IRResultType initializeFrom(InputRecord *ir);
    int checkConsistency();
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);
    
    /**
     * Compute vector of master contribution coefficients - SUMA of contributions == 1.0
     */
    int computeMasterContribution();

    const char *giveClassName() const { return "HangingNode"; }
    classType giveClassID() const { return HangingNodeClass; }
    bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_slave ); }
};
} // end namespace oofem
#endif // hangingnode_h
