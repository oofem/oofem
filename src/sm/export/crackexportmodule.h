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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef crackexportmodule_h
#define crackexportmodule_h

#include "exportmodule.h"
#include "floatarray.h"
#include "intarray.h"

///@name Input fields for MPSMaterial
//@{
#define _IFT_CrackExportModule_Name "crackvectorexport"
#define _IFT_CrackExportModule_cs "cs"
#define _IFT_CrackExportModule_threshold "threshold"
//#define _IFT_CrackExportModule_exporttype "exporttype"
//@}

namespace oofem {
/**
 * This one-purpose export module serves for estimation of the total water loss
 */
class CrackExportModule : public ExportModule
{
protected:
  IntArray crossSections;
  double threshold;

public:
    /// Constructor.
    CrackExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~CrackExportModule();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep, bool forcedOutput);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveClassName() const { return "CrackExportModule"; }
    virtual const char *giveInputRecordName() const { return _IFT_CrackExportModule_Name; }

    static void writeToOutputFile(const std :: string &iName, const std :: vector< FloatArray > &iPoints);
};

} // namespace oofem

#endif // crackexportmodule_h
