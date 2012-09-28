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

#ifndef homexportmodule_h
#define homexportmodule_h

#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"
#include "intarray.h"
#include "linearelasticmaterial.h"

namespace oofem {
/**
 * Represents HOM (Homogenization) export module. It averages strain and stress tensors over the whole domain
 * and all elements in global coordinate system. The strain and stress tensor (reduced to six components) is used through whole procedure.
 * Appropriate element strain and stress components are placed using the mask of tensor indexes.
 * Thus various element types (beam, plane element, brick) can be combined and will give correct macroscopic response.
 *
 * @author Vit Smilauer
 */
class HOMExportModule : public ExportModule, public LinearElasticMaterial
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
    virtual void doOutput(TimeStep *tStep);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveClassName() const { return "HOMExportModule"; }

protected:
    /// Stream for file.
    FILE *stream;
    /// Material numbers over which averaging is performed.
    IntArray matnum;
};
} // end namespace oofem

#endif
