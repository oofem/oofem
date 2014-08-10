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

#ifndef elementside_h
#define elementside_h

#include "dofmanager.h"

#define _IFT_ElementSide_Name "elementside"

namespace oofem {
class Dof;
class NodalLoad;
class TimeStep;
class FloatArray;
class IntArray;

/**
 * Class implementing element side having some DOFs in finite element mesh.
 * ElementSide possess degrees of freedom (see base class DofManager).
 * ElementSide id usually attribute of few elements and is managed by domain.
 *
 * Tasks:
 * - managing its degrees of freedom (giveDof).
 * - calculating its  load vector.
 * - printing and updating at end of step.
 * - managing its swapping to and from disk.
 */
class OOFEM_EXPORT ElementSide : public DofManager
{
public:
    /**
     * Constructor. Creates a element side belonging to domain.
     * @param n Side number in domain aDomain
     * @param aDomain Domain to which side belongs
     */
    ElementSide(int n, Domain * aDomain);
    /// Destructor.
    virtual ~ElementSide();

    // miscellaneous
    virtual const char *giveClassName() const { return "ElementSide"; }
    virtual const char *giveInputRecordName() const { return _IFT_ElementSide_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void printYourself();

    /**
     * Computes receiver transformation matrix from global cs. to dofManager specific
     * coordinate system (in which governing equations are assembled, for example the
     * local coordinate system in node).
     * @param answer Computed transformation matrix. It has generally dofIDArry.size rows and
     * if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
     * This is because this transformation should generally include not only transformation to
     * dof manager local coordinate system, but receiver dofs can be expressed using
     * dofs of another dofManager (In this case, square answer is produced only if all
     * dof transformation is required).
     * @param dofIDArry Array containing DofIDItem values for which transformation matrix is
     * assembled. If dofIDArry is NULL, then all receiver dofs are assumed.
     */
    virtual void computeTransformation(FloatMatrix &answer, const IntArray *dofIDArry);
    virtual bool requiresTransformation() { return false; }
    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_simpleSlave ); }
};
} // end namespace oofem
#endif // elementside_h
