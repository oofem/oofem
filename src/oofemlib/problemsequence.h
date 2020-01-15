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
    ProblemSequence(int i, EngngModel *master=nullptr);
    /// Destructor.
    virtual ~ProblemSequence();

    ProblemSequence(const ProblemSequence &) = delete;
    ProblemSequence & operator=(const ProblemSequence &) = delete;

    EngngModel & giveActiveModel() { return *emodelList[activeModel]; }

    void solveYourself() override;

    //virtual void initializeYourself(TimeStep *tStep);
    int instanciateYourself(DataReader &dr, InputRecord &ir, const char *outFileName, const char *desc) override;
    void initializeFrom(InputRecord &ir) override;
    int checkProblemConsistency() override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    // identification
    const char *giveClassName() const override { return "ProblemSequence"; }
    const char *giveInputRecordName() const { return _IFT_ProblemSequence_Name; }

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &gc) override;
    void drawElements(oofegGraphicContext &gc) override;
    void drawNodes(oofegGraphicContext &gc) override;
    void showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep) override { }
#endif

    EngngModel *giveSlaveProblem(int i) override { return NULL; }
    int giveNumberOfSlaveProblems() override { return 0; }

    int instanciateDefaultMetaStep(InputRecord &ir) override { return 1; }
};
} // end namespace oofem
#endif // problemsequence_h
