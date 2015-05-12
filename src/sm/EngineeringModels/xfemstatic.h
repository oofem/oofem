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

#ifndef XFEMSTATIC_H_
#define XFEMSTATIC_H_

#include "../sm/EngineeringModels/nlinearstatic.h"
#include "fracturemanager.h" 
#include <map>

///@name Input fields for XFEMStatic
//@{
#define _IFT_XFEMStatic_Name "xfemstatic"
#define _IFT_XFEMStatic_ForceRemap "forceremap"
//@}

//#define USE_FRACTURE_MANAGER

namespace oofem {
/**
 * Solver for XFEM simulations. The class inherits from NonLinearStatic and
 * adds support for evolving XFEM interfaces.
 *
 * The routines have to a large extent been copied from StaticFracture in Jim's repo.
 *
 * @author Erik Svenning
 */

class XFEMStatic : public NonLinearStatic
{
public:
    XFEMStatic(int i, EngngModel * _master = NULL);
    virtual ~XFEMStatic();
    virtual int requiresUnknownsDictionaryUpdate() { return updateStructureFlag; }
    virtual bool requiresEquationRenumbering(TimeStep *) { return updateStructureFlag; }

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);
    virtual void updateLoadVectors(TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);

    void initializeDofUnknownsDictionary(TimeStep *tStep);
    void setTotalDisplacementFromUnknownsInDictionary(ValueModeType mode, TimeStep *tStep);
    virtual void updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep);

    void setUpdateStructureFlag(bool flag) { updateStructureFlag = flag; }
    bool needsStructureUpdate() { return updateStructureFlag; }

    void buildDofMap();
    void setValsFromDofMap(FloatArray &oArray, const FloatArray &iArray);

protected:
    bool updateStructureFlag;

    bool mForceRemap;

    bool mSetValsFromDofMap;

    // Map for updating initialLoadVector after resize
    // Maps <domain index, node index, local dof index> to global equation number
    std :: map< std :: vector< int > , int >mDofEqnNumMap;

    // Jim
    FractureManager *fMan;
};
} /* namespace oofem */
#endif /* XFEMSTATIC_H_ */
