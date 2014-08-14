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

#ifndef staticstructuralstaggered_h
#define staticstructuralstaggered_h

#include "../sm/EngineeringModels/staticstructural.h"
#include "staggeredsolver.h"

#define _IFT_StaticStructuralStaggered_Name "staticstructuralstaggered"

namespace oofem {
class SparseMtrx;

/**
 * Solves a static structural problem using a staggered approach. Experimental version.
 * @author Jim Brouzoulis
 */
class StaticStructuralStaggered : public StaticStructural
{
protected:

    
    std :: vector< CustomEquationNumbering > UnknownNumberingSchemeList;    
    std :: vector< SparseMtrx *> stiffnessMatrixList;
    std :: vector< FloatArray > fIntList;
    std :: vector< FloatArray > fExtList;

    
public:
    StaticStructuralStaggered(int i, EngngModel * _master = NULL);
    virtual ~StaticStructuralStaggered();
    virtual IRResultType initializeFrom(InputRecord *ir);

    //virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);

    virtual void terminate(TimeStep *tStep);

    //virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    
    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);

    virtual int forceEquationNumbering();

    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual fMode giveFormulation() { return TL; }
   

    // identification
    virtual const char *giveInputRecordName() const { return _IFT_StaticStructuralStaggered_Name; }
    virtual const char *giveClassName() const { return "StaticStructuralStaggered"; }

    // Experimental Jim
    virtual void updateInternalForcesForStaggeredSolver(FloatArray &answer, TimeStep* tStep, Domain* d, const CustomEquationNumbering& s);
    virtual void updateExternalForcesForStaggeredSolver(FloatArray &answer, TimeStep* tStep, Domain* d, const CustomEquationNumbering& s);    
    virtual void updateTangentStiffnessForStaggeredSolver(SparseMtrx *answer, TimeStep* tStep, Domain* d, const CustomEquationNumbering& s);
    SparseMtrxType giveSparseMtrxType() { return this->sparseMtrxType; };
    
    
    
    virtual int giveNewEquationNumber(int domain, DofIDItem);
    virtual int giveNewPrescribedEquationNumber(int domain, DofIDItem);
    virtual int giveNumberOfEquations(int id, CustomEquationNumbering &num);
};
} // end namespace oofem
#endif // staticstructural_h
