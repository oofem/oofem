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

#ifndef macrolspace_h
#define macrolspace_h

#include "../sm/Elements/3D/lspace.h"

///@name Input fields for MacroLspace
//@{
#define _IFT_MacroLSpace_Name "macrolspace"
#define _IFT_MacroLspace_microMasterNodes "micromasternodes"
#define _IFT_MacroLspace_microBoundaryNodes "microboundarynodes"
#define _IFT_MacroLspace_stiffMatrxFileName "stiffmatrxfilename"
//@}

namespace oofem {
class Domain;
class EngngModel;
class MicroMaterial;

/**
 * This class implements a macroelement. It is derived from eight-node brick element.
 * The stiffness matrix is computed from underlying RVE and is condensed to 24 DoFs to corner nodes.
 *
 * @author Vit Smilauer
 */
class MacroLSpace : public LSpace
{
public:
    MacroLSpace(int n, Domain * d);
    virtual ~MacroLSpace();

    virtual const char *giveInputRecordName() const { return _IFT_MacroLSpace_Name; }
    virtual const char *giveClassName() const { return "MacroLSpace"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
    { OOFEM_ERROR("Macro space element doesn't support computing local unknown vector (yet)\n"); }

    /// Related to setting the boundary conditions of micro problem.
    virtual void changeMicroBoundaryConditions(TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    /**
     * Evaluates shape function at a given pointnodal representation of real internal forces obtained from microProblem.
     * @param answer Array of shape function values at given node.
     * @param coords Coordinates of nodes defining the interpolation geometry.
     * @param gcoords Global coordinates of point of interest.
     */
    virtual void evalInterpolation(FloatArray &answer, const std::vector< FloatArray > &coords, const FloatArray &gcoords);

    virtual void updateYourself(TimeStep *tStep);

protected:
    /// Array containing the node mapping from microscale (which microMasterNodes corresponds to which macroNode)
    IntArray microMasterNodes;
    IntArray microBoundaryNodes;
    IntArray microDOFs;
    bool firstCall;
    MicroMaterial *microMaterial;
    Domain *microDomain;
    EngngModel *microEngngModel;
    /// Information of iteration number.
    int iteration;
    /// Stores node number on the boundary in the triplets.
    IntArray microBoundaryDofManager;
    FloatMatrix stiffMatrix;
    /// Process with external file for the storage of stiffness matrix 0-None, 1-read, 2-write.
    int stiffMatrxFileNoneReadingWriting;
    /// Array containg the force vector from nodes (if condensation is skipped, use this vector).
    FloatArray internalMacroForcesVector;
    /// Last time step when stiffness matrix was assembled.
    TimeStep *lastStiffMatrixTimeStep;
};
} // end namespace oofem
#endif //macrolspace_h
