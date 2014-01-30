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

#ifndef GNUPLOTEXPORTMODULE_H_
#define GNUPLOTEXPORTMODULE_H_

#include "exportmodule.h"

///@name Input fields for MatlabExportModule
//@{
#define _IFT_GnuplotExportModule_Name "gnuplot"
// Sum of reaction forces for each Dirichlet BC
#define _IFT_GnuplotExportModule_ReactionForces "reactionforces"
//@}

namespace oofem {
class EnrichmentItem;
class Crack;
/**
 * (Under development) The Gnuplot export module enables OOFEM to export some
 * data in a format that can be directly plotted with Gnuplot.
 *
 * @author Erik Svenning
 *
 * Created on: Jan 29, 2014
 */
class OOFEM_EXPORT GnuplotExportModule : public ExportModule {
public:
	GnuplotExportModule(int n, EngngModel *e);
	virtual ~GnuplotExportModule();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
    virtual void initialize();
    virtual void terminate();

    virtual const char *giveClassName() const { return "GnuplotExportModule"; };
    virtual const char *giveInputRecordName() const { return _IFT_GnuplotExportModule_Name; }

    /**
     * XFEM output
     */
    void outputXFEM(EnrichmentItem &iEI);
    void outputXFEM(Crack &iCrack);

protected:
	bool mExportReactionForces;

    /**
     * Stores the sum of reaction forces for each BC.
     */
    std::vector< std::vector<FloatArray> > mReactionForceHistory;
    std::vector< std::vector<double> > mDispHist;

    void outputReactionForces(TimeStep *tStep);

};
} // end namespace oofem
#endif /* GNUPLOTEXPORTMODULE_H_ */
