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

#ifndef problemsequence_h
#define problemsequence_h

#include "engngm.h"
#include "inputrecord.h"

///@name Input fields for ProblemSequence
//@{
#define _IFT_ProblemSequence_Name "problemsequence"
#define _IFT_ProblemSequence_engineeringModels "engngms"
//@}

namespace oofem {
class Function;

/**
 * Meta-engineering problem used to solve a sequence off different problems, all using the same domain.
 * For example, a static deformation followed
 *
 * @todo This class is far from finished and is still lacking vital functionality. It does not work yet.
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT ProblemSequence : public EngngModel
{
protected:
    /// List of engineering models to solve sequentially.
    std :: vector< std :: unique_ptr< EngngModel > >emodelList;
    std :: vector< std :: string >inputStreamNames;

    /// Keeps track of the active model in the analysis sequence.
    int activeModel;

public:
    /// Constructor
    ProblemSequence(int i, EngngModel * _master = NULL);
    /// Destructor.
    virtual ~ProblemSequence();

    ProblemSequence(const ProblemSequence &) = delete;
    ProblemSequence & operator=(const ProblemSequence &) = delete;

    EngngModel & giveActiveModel() { return *emodelList[activeModel]; }

    virtual void solveYourself();

    //virtual void initializeYourself(TimeStep *tStep);
    virtual int instanciateYourself(DataReader *dr, InputRecord *ir, const char *outFileName, const char *desc);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkProblemConsistency();

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // identification
    virtual const char *giveClassName() const { return "ProblemSequence"; }
    virtual const char *giveInputRecordName() const { return _IFT_ProblemSequence_Name; }

#ifdef __OOFEG
    virtual void drawYourself(oofegGraphicContext &gc);
    virtual void drawElements(oofegGraphicContext &gc);
    virtual void drawNodes(oofegGraphicContext &gc);
    virtual void showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep) { }
#endif

    virtual EngngModel *giveSlaveProblem(int i) { return NULL; }
    virtual int giveNumberOfSlaveProblems() { return 0; }

    virtual int instanciateDefaultMetaStep(InputRecord *ir) { return 1; }
};
} // end namespace oofem
#endif // problemsequence_h
