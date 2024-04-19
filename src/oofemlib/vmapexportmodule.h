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

#ifndef vmapexportmodule_h
#define vmapexportmodule_h

#include "exportmodule.h"
#include "exportregion.h"
#include "domain.h"
#include "intarray.h"
#include "VMAP.h"
#include "VMAPFile.h"
#include <set>


///@name Input fields for VTK XML export module
//@{
#define _IFT_VMAPExportModule_Name "vmapem"
#define _IFT_VMAPExportModule_cellvars "cellvars"
#define _IFT_VMAPExportModule_vars "vars"
#define _IFT_VMAPExportModule_primvars "primvars"
#define _IFT_VMAPExportModule_ipvars "ipvars"
#define _IFT_VMAPExportModule_regionstoskip "regionstoskip"
#define _IFT_VMAPExportModule_nvr "nvr"
#define _IFT_VMAPExportModule_vrmap "vrmap"
//@}



namespace oofem {
 /**
 * Represents VMAP (https://www.scai.fraunhofer.de/en/business-research-areas/multiphysics/projects/itea-vmap.html) export module. 
 * VMAP is a vendor-neutral standard for CAE data storage to enhance interoperability in virtual engineering workflows.
 * 
 * The export of data is done on Region By Region basis, taking care about (possibly) discontinuous character of
 * some internal variables at region boundaries.
 */
class OOFEM_EXPORT VMAPExportModule : public ExportModule
{
protected:
    /// List of InternalStateType values, identifying the selected vars for export.
    IntArray internalVarsToExport;
    /// List of primary unknowns to export.
    IntArray primaryVarsToExport;
    /// List of cell data to export.
    IntArray cellVarsToExport;
    /// List of internal variables to export directly in Integration Points (no smoothing to nodes)
    IntArray ipInternalVarsToExport;

    /// Map from Voigt to full tensor.
    IntArray redToFull;

    /// List of regions to skip.
    IntArray regionsToSkip;
    /// Number of virtual regions.
    int nvr;
    /// Real->virtual region map.
    IntArray vrmap;
    /// Scaling time in output, e.g. conversion from seconds to hours
    double timeScale;

    // VMAP fileName
    std::string vmapFileName;

    // exportRegions
    std::vector< ExportRegion > exportRegions;

    // vmap element type map (per domain)
    std::map<int,int> vmapelementtype ;


public:
    /// Constructor. Creates empty VMAP Export Manager with number n.
    VMAPExportModule(int n, EngngModel *e);
    virtual ~VMAPExportModule();
    /// Initializes receiver according to object description stored in input record.
    void initializeFrom(InputRecord &ir) override;
    /**
     * Writes the output in given solution step. 
     * @param tStep time step.
     * @param bool if true, no testTimeStepOutput should be done
     */
    void doOutput(TimeStep *tStep, bool forcedOutput=false) override;
    /**
     * Initializes receiver.
     * The init file messages should be printed.
     */
    virtual void initialize() override;
    /**
     * Terminates the receiver.
     * The terminating messages should be printed.
     * All the streams should be closed.
     */
    virtual void terminate() override;
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const override { return "VMAPExportModule"; }
protected:
    void writeGeometry(VMAP::VMAPFile &file, Domain* d);
    int initRegionNodeNumbering(ExportRegion& piece, Domain *domain, Set& region);
    void setupExportRegion(ExportRegion &exportRegion, Set &region);
    void setupElementTypes(VMAP::VMAPFile& vmapFile);
    int getElementTypeHash(Element* e);
    void getIntegrationPointMapping (IntArray& ipmap, integrationDomain id, int npoints);
     /**
     * Returns corresponding element VMAP cell_type.
     * Some common element types are supported, others can be supported via interface concept.
    */
    int giveElementType(Element *element);
    int giveElementType(int num) ;
    bool isElementComposite(Element *elem); /// Returns true if element geometry type is composite (not a single cell).
    /**
     * Returns the element cell connectivity in target VMAP format.
     */
    void giveElementConnectivity(IntArray &answer, Element *elem);
    /**
     * @brief writes variable data
     * 
     */
    void writePrimaryVariable (VMAP::VMAPFile& file, Domain*d, int region, UnknownType ut, TimeStep* tStep);
    void  getNodalVariableFromPrimaryField(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, int region);
    void writeInternalVariable(VMAP::VMAPFile& file, Domain*d, int region, InternalStateType ut, TimeStep* tStep);
    std::string getGroupName (Domain* d, int region);
};

} // end namespace oofem
#endif // vmapexportmodule_h
