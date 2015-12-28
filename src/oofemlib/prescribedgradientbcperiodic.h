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


#ifndef prescribedgradientperiodic_h_
#define prescribedgradientperiodic_h_

#include "prescribedgradienthomogenization.h"
#include "activebc.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include <memory>
#include <map>

#define _IFT_PrescribedGradientBCPeriodic_Name "prescribedgradientperiodic"
#define _IFT_PrescribedGradientBCPeriodic_masterSet "masterset"
#define _IFT_PrescribedGradientBCPeriodic_jump "jump"

namespace oofem {
class Node;

/**
 * Prescribes an average displacement gradient based on microperiodicity.
 * @note Nodes on either master or slave side may NOT have local coordinate systems.
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT PrescribedGradientBCPeriodic : public ActiveBoundaryCondition, public PrescribedGradientHomogenization
{
protected:
    std :: unique_ptr< Node > strain;
    IntArray strain_id;

    std :: map< int, int > slavemap;
    FloatArray jump;
    
    int masterSet;

    /**
     * This is the central support function, which finds the corresponding slave nodes for each master node.
     */
    void findSlaveToMasterMap();

public:
    PrescribedGradientBCPeriodic(int n, Domain *d);
    virtual ~PrescribedGradientBCPeriodic();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual void postInitialize();

    virtual int giveNumberOfInternalDofManagers();
    virtual DofManager *giveInternalDofManager(int i);

    virtual int giveNumberOfMasterDofs(ActiveDof *dof);
    virtual Dof *giveMasterDof(ActiveDof *dof, int mdof);

    virtual void computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs);
    virtual void computeField(FloatArray &sigma, TimeStep *tStep);
    virtual void computeTangent(FloatMatrix &E, TimeStep *tStep);
    double giveUnknown(double val, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual double giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual bool isPrimaryDof(ActiveDof *dof);
    virtual double giveBcValue(Dof *dof, ValueModeType mode, TimeStep *tStep);
    virtual bool hasBc(Dof *dof, TimeStep *tStep);
    bool isStrainDof(Dof *dof);

    virtual bool requiresActiveDofs() { return true; }

    virtual const char *giveClassName() const { return "PrescribedGradientBCPeriodic"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCPeriodic_Name; }

};
} // end namespace oofem

#endif // prescribedgradientperiodic_h_
