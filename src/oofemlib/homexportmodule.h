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
#include "domain.h"
#include "engngm.h"
#include "intarray.h"


///@name Input fields for Homogenization export module
//@{
#define _IFT_HOMExportModule_Name "hom"
#define _IFT_HOMExportModule_scale "scale"
#define _IFT_HOMExportModule_matnum "matnum"
//@}

namespace oofem {
/**
 * Represents HOM (Homogenization) export module. It averages strain and stress tensors over the whole domain
 * and all elements in global coordinate system. The strain and stress tensor (reduced to six components) is used through whole procedure.
 * Appropriate element strain and stress components are placed using the mask of tensor indexes.
 * Thus various element types (beam, plane element, brick) can be combined and will give correct macroscopic response.
 *
 * @author Vit Smilauer
 */
class HOMExportModule : public ExportModule
{
protected:
    /// Scale of all homogenized values.
    double scale;
#ifdef RBR_SUPPORT
    enum omodeType { wdmode, rbrmode }; // WholeDomain or RegionByRegion output
    omodeType omode;
#endif

public:
    /// Constructor. Creates empty Output Manager.
    HOMExportModule(int n, EngngModel *e);
    /// Destructor.
    virtual ~HOMExportModule();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep, bool forcedOutput=false);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveClassName() const { return "HOMExportModule"; }
    virtual const char *giveInputRecordName() const { return _IFT_HOMExportModule_Name; }

protected:
    /// Stream for file.
    FILE *stream;
    /// Material numbers over which averaging is performed.
    IntArray matnum;
};
} // end namespace oofem

#endif
