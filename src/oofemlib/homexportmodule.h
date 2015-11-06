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

#ifndef homexportmodule_h
#define homexportmodule_h

#include "exportmodule.h"
#include "floatarray.h"

///@name Input fields for Homogenization export module
//@{
#define _IFT_HOMExportModule_Name "hom"
#define _IFT_HOMExportModule_ISTs "ists" /// List of internal state types used for output
#define _IFT_HOMExportModule_scale "scale" ///[optional] Scales the output variables
//#define _IFT_HOMExportModule_matnum "matnum" ///[optional] If specified, only these materials are used
//@}

namespace oofem {
/**
 * Represents HOM (Homogenization) export module. It averages internal variables over the whole domain
 * and all elements in global coordinate system. Tensors are given in Voigt notation.
 * Thus various element types (beam, plane element, brick) can be combined and will give correct averaged values.
 *
 * @author Vit Smilauer
 * @author Mikael Öhman
 */
class OOFEM_EXPORT HOMExportModule : public ExportModule
{
protected:
    /// Scale of all homogenized values.
    double scale;
    /// Stream for file.
    FILE *stream;
    /// Material numbers over which averaging is performed. - replaced by 'regionsets'
    //IntArray matnum;
    /// Internal states to export
    IntArray ists;

public:
    /// Constructor. Creates empty Output Manager.
    HOMExportModule(int n, EngngModel * e);
    /// Destructor.
    virtual ~HOMExportModule();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveClassName() const { return "HOMExportModule"; }
    virtual const char *giveInputRecordName() const { return _IFT_HOMExportModule_Name; }
};
} // end namespace oofem

#endif
