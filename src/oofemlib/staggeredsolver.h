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

#ifndef staggeredsolver_h
#define staggeredsolver_h

#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "linesearch.h"
#include "nrsolver.h"
#include "unknownnumberingscheme.h"

#include <set>
#include <vector>


///@name Input fields for StaggeredSolver
//@{
#define _IFT_StaggeredSolver_Name "staggeredsolver"
#define _IFT_StaggeredSolver_DofIdList "dofidlist"
#define _IFT_StaggeredSolver_DofIdListPositions "idpos"
//@}

namespace oofem {
class Domain;
class EngngModel;


class CustomEquationNumbering : public UnknownNumberingScheme
{
protected:
    bool prescribed;
    int numEqs;
    int numPresEqs;
    int number;

public:
    CustomEquationNumbering() : UnknownNumberingScheme(), prescribed(false), numEqs(0), numPresEqs(0) { dofIdArray.resize(0); }

    IntArray dofIdArray; // should be private
    virtual bool isDefault() const { return true; }
    void setDofIdArray(IntArray &array) { this->dofIdArray = array; }
    void setNumber(int num) { this->number = num; }
    int getNumber() { return this->number; }
    
    virtual int giveDofEquationNumber(Dof *dof) const {
        DofIDItem id = dof->giveDofID();
        //printf("asking for num %d \n", (int)id);
            if ( this->dofIdArray.contains( (int)id ) ) {
            return prescribed ? dof->__givePrescribedEquationNumber() : dof->__giveEquationNumber();    
        } else {
            return 0;  
        }
    }
    

    virtual int giveRequiredNumberOfDomainEquation() const { return numEqs; }

    int giveNewEquationNumber() { return ++numEqs; }
    int giveNewPrescribedEquationNumber() { return ++numPresEqs; }
    int giveNumEquations() { return this->numEqs; }
    int giveNumPresEquations() { return this->numPresEqs; }
};

/**
 * The staggered solver will perform Newton iterations on subsets of DofIDs, in a staggered manner.
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT StaggeredSolver : public NRSolver
{
private:
    IntArray totalIdList;
    IntArray idPos;
  
    // Lists for each dof group
    std :: vector< CustomEquationNumbering > UnknownNumberingSchemeList;
    std :: vector< std :: unique_ptr< SparseMtrx > > stiffnessMatrixList;
    std :: vector< FloatArray > fIntList;
    std :: vector< FloatArray > fExtList;
    
    std :: vector< IntArray > locArrayList;
    std :: vector< FloatArray > ddX;
    std :: vector< FloatArray > dX;
    std :: vector< FloatArray > X;
    

    void giveTotalLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s, Domain *d);
    bool checkConvergenceDofIdArray(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, TimeStep *tStep, IntArray &dofIdArray);

    void instanciateYourself();

public:
    StaggeredSolver(Domain * d, EngngModel * m);
    virtual ~StaggeredSolver() {}

    // Overloaded methods:
    virtual NM_Status solve(SparseMtrx &k, FloatArray &R, FloatArray *R0, FloatArray *iR,
                            FloatArray &X, FloatArray &dX, FloatArray &F,
                            const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *);

    virtual IRResultType initializeFrom(InputRecord *ir);
    
    virtual const char *giveClassName() const { return "StaggeredSolver"; }
    virtual const char *giveInputRecordName() const { return _IFT_StaggeredSolver_Name; }
};
} // end namespace oofem
#endif // staggeredsolver_h
