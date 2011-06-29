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

#ifndef externalanalysis_h
#define externalanalysis_h

#include "inputrecord.h"

namespace oofem {
class EngngModel;
class Domain;
class TimeStep;

/**
 * States of the external analysis.
 */
enum ExternalAnalysisState {
    EAS_OK, ///< External analysis terminated successfully.
    EAS_Failed, ///< External analysis terminated with an error.
    EAS_NotConverged, ///< External analysis requires additional iterations.
};

/**
 * Different types of analysis.
 * @todo{This could possibly be moved to the engineering model class.}
 */
enum AnalysisType {
    AT_Static,  ///< Stiffness with no acceleration.
    AT_Dynamic, ///< Stiffness + mass.
    AT_AdaptiveStatic,  ///< Stiffness with no acceleration with dynamically changing mesh.
    AT_AdaptiveDynamic, ///< Stiffness + mass with dynamically changing mesh.
};

/**
 * Abstract class for coupling to external analysis tools for staggered solvers, typical for fluid structure interaction.
 * The external analysis is an attribute of an engineering model and is coupled to a specific domain.
 * @note{Implementations will need to specifically support recursive analysis types}
 * @note{This could possible also be used for normal staggered analysis within OOFEM modules}
 * @note{Under development; design is subject to change}
 * @author Mikael Öhman
 */
class ExternalAnalysis
{
protected:
    /// Master engineering problem.
    EngngModel *m;
    /// Domain index, default is 1.
    int di;

public:
    ExternalAnalysis(EngngModel *m): m(m), di(1) { }
    ~ExternalAnalysis() { }

    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    /// Checks internal consistency.
    virtual bool checkConsistency() {
        /*AnalysisType at = this->d->giveEngngModel()->giveAnalysisType();
        if ( this->supportsAnalysisType(at) ) {
            OOFEM_WARNING3("%s :: AnalysisType %d not supported.", giveClassName(), at);
            return false;
        }*/
        return true;
    }
    /// Sets a new domain index. Useful during adaptive refinements or similar dynamic domain changes.
    void setDomain(int di) { this->di = di; }

    /**
     * Exports the solution from the given domain and runs the external analysis.
     * @return State of the external analysis.
     */
    virtual ExternalAnalysisState solveYourselfAt(TimeStep *tStep) = 0;
    /**
     * Updates state when equilibrium has been reached.
     * This routine also involves writing result files and optionally saving context.
     * @param tStep Time step which has been finalized.
     */
    virtual void updateYourself(TimeStep *tStep);

    /**
     * Determines if the external analysis can be coupled to the specific analysis type.
     * @return True iff type is supported.
     */
    virtual bool supportsAnalysisType(AnalysisType at) = 0;
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "ExternalAnalysis"; }
};
} // end namespace oofem
#endif // externalanalysis_h
