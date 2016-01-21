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

#ifndef SOLUTIONBASEDSHAPEFUNCTION_H_
#define SOLUTIONBASEDSHAPEFUNCTION_H_

#include "activebc.h"
#include "dofiditem.h"
#include "floatarray.h"

#define _IFT_SolutionbasedShapeFunction_Name "solutionbasedshapefunction"
#define _IFT_SolutionbasedShapeFunction_ShapeFunctionFile "shapefunctionfile"
#define _IFT_SolutionbasedShapeFunction_Externalset "externalset"
#define _IFT_SolutionbasedShapeFunction_UseCorrectionFactors "usecorrectionfactors"
#define _IFT_SolutionbasedShapeFunction_DumpSnapshots "dumpsnapshots"

namespace oofem {
class EngngModel;
class Element;
class Node;

struct SurfaceDataStruct {
    bool isPlus;
    bool isMinus;
    bool isZeroBoundary;
    bool isFree;
    DofIDItem DofID;
    DofManager *DofMan;
    double value;
};

struct modeStruct {
    double am, ap;
    EngngModel *myEngngModel;
    std :: vector< SurfaceDataStruct * >SurfaceData;
};

class OOFEM_EXPORT SolutionbasedShapeFunction : public ActiveBoundaryCondition
{
private:
    Node *myNode;
    IntArray myDofIDs;
    int externalSet;
    int order;
    double TOL;
    std :: string filename;
    bool useConstantBase;
    bool isLoaded;

    bool useCorrectionFactors;
    bool dumpSnapshot;

    double bigNorm;
    int worstIndex;

    bool isCoeff(ActiveDof *dof);

    std :: vector< modeStruct * >modes;

    void updateModelWithFactors(modeStruct *m);
    TimeStep *thisTimestep;

    FloatArray maxCoord, minCoord;

    void setBoundaryConditionOnDof(Dof *d, double value);

    void setLoads(EngngModel *myEngngModel, int d);
    void loadProblem();
    void init();

    void computeCorrectionFactors(modeStruct &myMode, IntArray *Dofs, double *am, double *ap);

    /**
     * @brief giveValueAtPoint
     * @param answer        Values as point corresponding to the coordinate and dofID
     * @param coords        Coordinated at which evaluation takes place
     * @param dofID         Specifies DofIDs to evaluate
     * @param myEngngModel  EngngModel
     */
    void giveValueAtPoint(FloatArray &answer, const FloatArray &coords, IntArray &dofID, EngngModel &myEngngModel);

    void giveCorrectedValueAtPoint(FloatArray &answer, const FloatArray &coords, IntArray &dofID, EngngModel &myEngngModel);

    void splitBoundaryNodeIDs(modeStruct &mode, Element &e, IntArray &boundary, IntArray &pList, IntArray &mList, IntArray &zList, FloatMatrix &nodeValues);

    void computeBaseFunctionValueAt(FloatArray &answer, FloatArray &coords, IntArray &dofIDs, EngngModel &myEngngModel);

    void initializeSurfaceData(modeStruct *mode);

    void copyDofManagersToSurfaceData(modeStruct *mode, IntArray nodeList, bool isPlus, bool isMinus, bool isZero);

    void whichBoundary(FloatArray &coord, bool &isPlus, bool &isMinus, bool &isZero);

public:
    SolutionbasedShapeFunction(int n, Domain * d);
    virtual ~SolutionbasedShapeFunction();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual bool requiresActiveDofs() { return true; }
    virtual int giveNumberOfInternalDofManagers() { return 1; }
    virtual DofManager *giveInternalDofManager(int i);

    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual double giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof);

    virtual bool hasBc(Dof *dof, TimeStep *tStep) { return false; }

    virtual bool isPrimaryDof(ActiveDof *dof) { return false; }

    virtual void computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs);

    virtual int giveNumberOfMasterDofs(ActiveDof *dof);
    virtual Dof *giveMasterDof(ActiveDof *dof, int mdof);

    virtual const char *giveClassName() const { return "SolutionbasedShapeFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_SolutionbasedShapeFunction_Name; }

    double checkIncompressibility(EngngModel &myEngngModel);
};
}

#endif /* SOLUTIONBASEDSHAPEFUNCTION_H_ */
