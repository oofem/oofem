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

#ifndef jsonexportmodule_h_
#define jsonexportmodule_h_

#include "outputexportmodule.h"
#include "json.h"

#define _IFT_JsonExportModule_Name "jsonoutput"

namespace oofem {

class OOFEM_EXPORT JsonExportModule : public OutputExportModule
{
    std::unique_ptr<JsonContext> ctx;
public:
    JsonExportModule(int n, EngngModel * e);
    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    void terminate() override;

    const char *giveClassName() const override { return "JsonExportModule"; }
    const char *giveInputRecordName() const { return _IFT_JsonExportModule_Name; }
};
} // end namespace oofem
#endif // outputexportmodule_h_
