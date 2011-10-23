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

#ifndef refinedelement_h
#define refinedelement_h

#include "alist.h"
#include "intarray.h"
#include "valuemodetype.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;
class Node;

class RefinedElement
{
protected:
    int elementId;
    AList< IntArray >fineNodeList;
    IntArray boundaryFlag;

public:
    RefinedElement(Domain *d, int elem, int level);
    ~RefinedElement();

    IntArray *giveFineNodeArray(int node);
    IntArray *giveBoundaryFlagArray(void) { return & boundaryFlag; }

    void giveBoundaryFlagArray(int inode, Element *element, IntArray &answer);

    bool giveBoundaryLoadArray1D(int inode, Element *element, IntArray &boundaryLoadArray);
    bool giveBoundaryLoadArray2D(int inode, Element *element, AList< IntArray > &boundaryLoadList);
    bool giveBoundaryLoadArray3D(int inode, Element *element, AList< IntArray > &boundaryLoadList);

    bool giveBcDofArray1D(int inode, Element *element, IntArray *sideBcDofId, int &sideNumBc, TimeStep *tStep);
    bool giveBcDofArray2D(int inode, Element *element, AList< IntArray > &sideBcDofIdList, IntArray &sideNumBc, TimeStep *tStep);
    bool giveBcDofArray3D(int inode, Element *element, AList< IntArray > &sideBcDofIdList, IntArray &sideNumBc,
                          AList< IntArray > &faceBcDofIdList, IntArray &faceNumBc, TimeStep *tStep);

protected:
    /**
     * Extract from dofArray of slave_node those Dofs that have compatible BCs with master_node
     * @param master_node Node to which Dof compatibility will be compared.
     * @param slave_node Node with original Dofs.
     * @param dofArray Array of ids of Dofs of slave_node to chose from.
     * @param dofs Number of Dofs in dofArray.
     * @param answer Array of ids of Dofs in dofArray with compatible BCs.
     * @param mode Mode of Dof values.
     * @param tStep Active time step.
     * @return Number of Dofs with compatible BCs.
     */
    int giveCompatibleBcDofArray(Node *master_node, Node *slave_node, IntArray &dofArray, int dofs, IntArray *answer,
                                 ValueModeType mode, TimeStep *tStep);

    /// Prints simple error message and exits.
    void error(const char *file, int line, const char *format, ...) const;
};
} // end namespace oofem
#endif // refinedelement_h
