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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef solutionstatusexportmodule_h_
#define solutionstatusexportmodule_h_

#include "exportmodule.h"

///@name Input fields for SolutionStatusExportModule
//@{
#define _IFT_SolutionStatusExportModule_Name "solutionstatus"
#define _IFT_SolutionStatusExportModule_format "fmt" ///< Filename where rules are defined (normally the input file).
//@}

namespace oofem {
class Domain;
class Element;
class DofManager;


/**
 * Configurable solution status export module. Creates and continuously updates the status file
 according to simulation progress.
 */
class OOFEM_EXPORT SolutionStatusExportModule : public ExportModule
{
protected:
  std::string filename;
  FILE* outputFile;
  std::vector<std::string> recs;

    void checkRecs();
    void printRecsHeader();
public:
    SolutionStatusExportModule(int n, EngngModel * e, FILE* out = nullptr);
    SolutionStatusExportModule(const SolutionStatusExportModule &) = delete;
    SolutionStatusExportModule &operator=(const SolutionStatusExportModule &) = delete;

    void initialize() override;
    void terminate() override;
    void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;

    const char *giveClassName() const override { return "SolutionStatusExportModule"; }
    const char *giveInputRecordName() const { return _IFT_SolutionStatusExportModule_Name; }
};
} // end namespace oofem
#endif // solutionstatusexportmodule_h_
