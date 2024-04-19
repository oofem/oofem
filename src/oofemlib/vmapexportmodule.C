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

#include "vmapexportmodule.h"
#include "exportmodule.h"
#include "engngm.h"
#include "domain.h"
#include "node.h"
#include "dof.h"
#include "timestep.h"
#include "cltypes.h"
#include "classfactory.h"
#include "oofemcfg.h" // for oofem version
#include "VMAP.h"
#include "VMAPFile.h"
#include "VMAPIntegrationTypeFactory.h"
#include "VMAPElementTypeFactory.h"

namespace oofem {
REGISTER_ExportModule(VMAPExportModule)

VMAPExportModule :: VMAPExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{
   
}


VMAPExportModule :: ~VMAPExportModule()
{

}


void
VMAPExportModule :: initializeFrom(InputRecord &ir)
{
    ExportModule :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, _IFT_VMAPExportModule_cellvars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_VMAPExportModule_vars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_VMAPExportModule_primvars); // Macro - see unknowntype.h
    IR_GIVE_OPTIONAL_FIELD(ir, ipInternalVarsToExport, _IFT_VMAPExportModule_ipvars); // Macro - see internalstatetype.h

   /*
    // obsolete, replaced by user controlled regionSets, see exportmodule 
    regionsToSkip.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, regionsToSkip, _IFT_VMAPExportModule_regionstoskip); // Macro


    this->nvr = 0; // number of virtual regions
    IR_GIVE_OPTIONAL_FIELD(ir, nvr, _IFT_VMAPExportModule_nvr); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, vrmap, _IFT_VMAPExportModule_vrmap); // Macro
    */
}

void 
VMAPExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput( tStep ) || forcedOutput ) ) {
        return;
    }

    VMAP::VMAPFile file(vmapFileName, VMAP::VMAPFile::OPENREADWRITE);
    int ireg = 1;
    // export primary variables (as nodal values)
    for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
        this->writePrimaryVariable(file, this->emodel->giveDomain(1), ireg, 
                                   (UnknownType) primaryVarsToExport.at(i), tStep);
    }
    // export internal variables (as integration point values)
    for ( int j = 1; j <= internalVarsToExport.giveSize(); j++ ) {   
        this->writeInternalVariable(file, this->emodel->giveDomain(1), ireg, 
                                   (InternalStateType) internalVarsToExport.at(j), tStep);
    }

    file.closeFile(); 
}

 
void 
VMAPExportModule :: writeGeometry(VMAP::VMAPFile& vmapFile, Domain* d)
{
    fprintf (stderr, "VMAPExportModule :: writeGeometry\n");
    //std :: string outputFileName = this->emodel->giveOutputBaseFileName() + ".h5";
    //VMAP::VMAPFile vmapFile(outputFileName, VMAP::VMAPFile::OPENREADWRITE);
    int nreg = this->giveNumberOfRegions();
    this->exportRegions.clear();
    for ( int ireg = 1; ireg <= nreg; ireg++ ) {
        this->exportRegions.push_back(ExportRegion());
        ExportRegion &er = this->exportRegions[ireg-1];
        Set* regionSet = this->giveRegionSet(ireg);
        this->setupExportRegion(er, *regionSet);
 

        // create one geometry group per domain
        std::string groupName = "D" + std::to_string(d->giveNumber())+"R"+std::to_string(ireg);
        std::string geometryGroupPath = vmapFile.createGeometryGroup(d->giveNumber(),groupName);
        fprintf (stderr, "\tGeometry group:%s, Path:%s\n", groupName.c_str(), geometryGroupPath.c_str());

        // export domain nodes
        int ndofman = er.giveNumberOfNodes();
        IntArray& l2gdofmap = er.getMapL2G();
        VMAP::sPointsBlock points(ndofman);
        for(int i = 1; i <= ndofman; i++) {
                Node *n = d->giveNode(l2gdofmap.at(i));
                points.setPoint(i-1, d->giveNode(i)->giveLabel(), n->giveCoordinate(1), n->giveCoordinate(2), n->giveCoordinate(3));
        }
        // write points 
        vmapFile.writePointsBlock(geometryGroupPath,points);


        
        // write elements of the region
        int nelem = er.giveNumberOfCells();
        IntArray& re = er.getRegionCells();
        IntArray ec;
        VMAP::sElementBlock eblock =  VMAP::sElementBlock(nelem);
        for (int i=1; i<= nelem; i++) {
            Element* e = d->giveElement(re.at(i));
            VMAP::sElement se = VMAP::sElement();
            se.setIdentifier(e->giveLabel());
            /*
            Element type -> (VMAP::sElementType id) the one set up and reqistered. Need to register element type for each oofem element type used
            element type=geometry+interpolation+integration rule.
            The combination of gemetry and interpolation in oofem unique, problem is with integration (user controlled)
            Need mapping between ooofem egt type+interpolation+intagration -> elementTYpe (hash function?)
            another posibility mapping element_name -> elementType (assuming default integration rule)
            */
            se.setElementType(this->giveElementType(e));
            this->giveElementConnectivity(ec, e);
            se.setConnectivity(ec.getStdVector());
            se.setCoordinateSystem(1); // assume Global coordinate system with identifier 1
            fprintf(stdout, "Element %d: ", i);
            for (auto i : se.getConnectivity()) {
                fprintf(stdout, "%d ", i);
            }
            fprintf(stdout, "\n");
            ec.printYourself();
            eblock.setElement(i-1, se);
        }
        vmapFile.writeElementsBlock(geometryGroupPath, eblock);
    }

    //vmapFile.closeFile();
}

int
VMAPExportModule::getElementTypeHash(Element* e)
{
    Element_Geometry_Type egt= e->giveGeometryType();
    int nip = e->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(); 
    return egt*1000+nip;
}

void 
VMAPExportModule::setupElementTypes (VMAP::VMAPFile& vmapFile)
{
    // each element type (=geometry+interpolation+integration_rule) needs to be setup in vmap file
    Domain *d  = emodel->giveDomain(1);
    int nelem = d->giveNumberOfElements();
    std::vector<VMAP::sIntegrationType> itypes;
    std::vector<VMAP::sElementType> etypes;

    int etypeCounter = 0;
    int itypeCounter = 0;
    vmapelementtype.clear(); // clear hash map of element types 
    
    for (int i=1; i<= nelem; i++) { // loop over elements
        Element* e = d->giveElement(i);
        Element_Geometry_Type egt= e->giveGeometryType();
        int nip = e->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(); 
        // make a hash as unique combination of oofem element geon=metry type and integration
        int hash = this->getElementTypeHash(e);
        //if (!vmapelementtype.contains(hash)) {
        auto found = vmapelementtype.find(hash);
        if (found == vmapelementtype.end()) {
            // hash not present, set up and register a new vmap element type
            // first determine vmap eIntegrationType
            int eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::INVALID_TYPE;
            int ees = VMAP::sElementType::eElementShape::SHAPE_INVALID;
            int eed = VMAP::sElementType::eElementDimension::ELEM_INVALID;
            int eint = VMAP::sElementType::eFieldInterpolationType::CONSTANT;
            if (egt == EGT_line_1 || egt == EGT_line_2) {
                if (egt == EGT_line_1) {
                    ees= VMAP::sElementType::eElementShape::LINE_2;
                    eint = VMAP::sElementType::eFieldInterpolationType::LINEAR;
                } else {
                    ees= VMAP::sElementType::eElementShape::LINE_3;
                    eint = VMAP::sElementType::eFieldInterpolationType::QUADRATIC;
                }
                eed = VMAP::sElementType::eElementDimension::ELEM_3D;
                if (nip == 1) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_1;
                } else if (nip == 2) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_3;
                } else if (nip == 3) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_4;
                } else {
                    OOFEM_ERROR ("setupElementTypes: unsupported integration type on element %d", i);
                }
            }
            if (egt == EGT_triangle_1 || egt== EGT_triangle_2 ) {
                if (egt == EGT_triangle_1) {
                    ees= VMAP::sElementType::eElementShape::TRIANGLE_3;
                    eint = VMAP::sElementType::eFieldInterpolationType::LINEAR;
                } else {
                    ees= VMAP::sElementType::eElementShape::TRIANGLE_6;
                    eint = VMAP::sElementType::eFieldInterpolationType::QUADRATIC;
                }   
                eed = VMAP::sElementType::eElementDimension::ELEM_3D;
                if (nip == 1) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_TRIANGLE_1;
                } else if (nip == 3) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_TRIANGLE_3;
                } else if (nip == 4) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_TRIANGLE_4;
                } else if (nip == 6) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_TRIANGLE_6;
                } else {
                    OOFEM_ERROR ("setupElementTypes: unsupported integration type on element %d", i);
                }
            } else if (egt == EGT_quad_1 || egt== EGT_quad_2) {
                if (egt == EGT_quad_1 ) {
                    ees= VMAP::sElementType::eElementShape::QUAD_4;
                    eint = VMAP::sElementType::eFieldInterpolationType::BILINEAR;
                } else {
                    ees= VMAP::sElementType::eElementShape::QUAD_8;
                    eint = VMAP::sElementType::eFieldInterpolationType::BIQUADRATIC;
                }
                eed = VMAP::sElementType::eElementDimension::ELEM_3D;
                if (nip == 1) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_QUAD_1;
                } else if (nip == 4) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_QUAD_4;
                } else if (nip == 9) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_QUAD_9;
                } else {
                    OOFEM_ERROR ("setupElementTypes: unsupported integration type on element %d", i);
                }
            } else if (egt == EGT_tetra_1 || egt == EGT_tetra_2 ) {
                if (egt == EGT_tetra_1) {
                    ees= VMAP::sElementType::eElementShape::TETRAHEDRON_4;
                    eint = VMAP::sElementType::eFieldInterpolationType::LINEAR;
                } else {
                    ees= VMAP::sElementType::eElementShape::TETRAHEDRON_10;
                    eint = VMAP::sElementType::eFieldInterpolationType::QUADRATIC;
                }
                eed = VMAP::sElementType::eElementDimension::ELEM_3D;
                if (nip == 1) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_TETRAHEDRON_1;
                } else if (nip == 4) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_TETRAHEDRON_4;
                } else if (nip == 8) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_TETRAHEDRON_8;
                } else {
                    OOFEM_ERROR ("setupElementTypes: unsupported integration type on element %d", i);
                }
            } else if (egt == EGT_hexa_1 || egt == EGT_hexa_2 || egt == EGT_hexa_27) {
                if (egt == EGT_hexa_1) {
                    ees= VMAP::sElementType::eElementShape::HEXAHEDRON_8;
                    eint = VMAP::sElementType::eFieldInterpolationType::TRILINEAR;
                } else if (egt == EGT_hexa_2) {
                    ees= VMAP::sElementType::eElementShape::HEXAHEDRON_20;
                    eint = VMAP::sElementType::eFieldInterpolationType::TRIQUADRATIC;
                } else {
                    ees= VMAP::sElementType::eElementShape::HEXAHEDRON_27;
                    eint = VMAP::sElementType::eFieldInterpolationType::TRICUBIC;
                }
                eed = VMAP::sElementType::eElementDimension::ELEM_3D;
                if (nip == 1) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_HEXAHEDRON_1;
                } else if (nip == 8) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_HEXAHEDRON_8;
                } else if (nip == 27) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_HEXAHEDRON_27;
                } else {
                    OOFEM_ERROR ("setupElementTypes: unsupported integration type on element %d", i);
                }

            } else if (egt == EGT_wedge_1 || egt == EGT_wedge_2) {
                if (egt == EGT_wedge_1) {
                    ees= VMAP::sElementType::eElementShape::WEDGE_6;
                    eint = VMAP::sElementType::eFieldInterpolationType::LINEAR;
                } else {
                    ees= VMAP::sElementType::eElementShape::WEDGE_15;
                    eint = VMAP::sElementType::eFieldInterpolationType::QUADRATIC;
                }
                eed = VMAP::sElementType::eElementDimension::ELEM_3D;
                if (nip == 1) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_WEDGE_1;
                } else if (nip == 6) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_WEDGE_6;
                } else if (nip == 8) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_WEDGE_8;
                } else if (nip == 9) {
                    eit = VMAP::VMAPIntegrationTypeFactory::eIntegrationTypeId::GAUSS_WEDGE_9;
                } else {
                    OOFEM_ERROR ("setupElementTypes: unsupported integration type on element %d", i);
                }

            } else {
                OOFEM_ERROR ("setupElementTypes: unsupported element geometry type (element %d)", i);
            }

            // create vmap IntegrationType
            VMAP::sIntegrationType it = VMAP::VMAPIntegrationTypeFactory::createVMAPIntegrationType(eit);
            it.setIdentifier(itypeCounter);   
            itypes.push_back(it);

            VMAP::sElementType et = VMAP::VMAPElementTypeFactory::createVMAPElementType(eed, ees, eint, itypeCounter);
            et.setIdentifier(etypeCounter);
            etypes.push_back(et);
            // remember generated element type
            vmapelementtype[hash]=etypeCounter++;
            itypeCounter++;
        }
    }
    vmapFile.writeIntegrationTypes(itypes);
    vmapFile.writeElementTypes(etypes);

}

void
VMAPExportModule :: initialize()
{
    fprintf  (stderr, "VMAPExportModule :: initialize");
    ExportModule::initialize();
    VMAP::Initialize();

    // construct vmap container file name
    char fext[100];
	sprintf( fext, ".m%d.h5", this->number);
    this->vmapFileName = this->emodel->giveOutputBaseFileName() + fext;

    VMAP::VMAPFile file(vmapFileName, VMAP::VMAPFile::CREATEORREPLACE);
    // nitialize meta info
    VMAP::sMetaInformation info;
    info.setExporterName(PRG_VERSION);
    info.setAnalysisType(this->emodel->giveClassName());
    file.writeMetaInformation (info);
    // initialize unit system, use ISO default
    VMAP::sUnitSystem sus;
    file.writeUnitSystem(sus);
    // initialize default c.s.
    /*
    VMAP::sCoordinateSystem scs;
    scs.setType(VMAP::sCoordinateSystem::CARTESIAN_RIGHT_HAND);
    scs.setIdentifier(1);
    file.writeCoordinateSystem
    */
    // optionally define  sUnit structure to define the derived VMAP units (this needs to be done later, when one knows variables to export)
    // set up element types in VMAP
    this->setupElementTypes(file); 
    // write geometry
    writeGeometry(file, this->emodel->giveDomain(1));

    file.closeFile();
}


void
VMAPExportModule :: terminate()
{ }


void
VMAPExportModule:: writePrimaryVariable (VMAP::VMAPFile& file, Domain*d, int region, UnknownType ut, TimeStep* tStep)
{
    ExportRegion &er = this->exportRegions[region-1];
    IntArray& l2gmap = er.getMapL2G();

    VMAP::sStateVariable sv;
    sv.setIdentifier(int(ut)); // @todo: identifier should be unique across UnknownType and IntVarType
    sv.setVariableName(__UnknownTypeToString(ut));
    sv.setTimeValue(tStep->giveTargetTime());
    //sv.setUnit(int unit); 
    //sv.setCoordinateSystem(int csid) // @todo
    //sv.setDimension(int dim)
    InternalStateValueType isvt = giveInternalStateValueType(ut);
    int varsize = 1;
    VMAP::sStateVariable::eVariableDimension vd;
    if (isvt == ISVT_SCALAR) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_SCALAR;
    } else if (isvt == ISVT_VECTOR) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_VECTOR;
        varsize = 3;
    } else if (isvt == ISVT_TENSOR_S3 || isvt == ISVT_TENSOR_S3E) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_2ND_ORDER_TENSOR_SYMMETRIC;
        varsize = 6;
    } else if (isvt == ISVT_TENSOR_G) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_2ND_ORDER_TENSOR;
        varsize = 9;
    } else {
        OOFEM_ERROR ("VMAPExportModule::writePrimaryVariable: unsupported InternalStateValueType")
    }
    sv.setDimension(vd);
    sv.setMultiplicity(1);
    sv.setLocation(VMAP::sStateVariable::eVariableLocation::LOCATION_NODE); // primary values always in nodes
    sv.setEntity(VMAP::sStateVariable::eVariableEntity::ENTITY_REAL);
    // store nodal values
    std::vector<double> values;
    int nnodes = this->exportRegions[region-1].giveNumberOfNodes();
    values.reserve(nnodes*varsize);

    //    @todo
    FloatArray nodalvalues;
    for(int i = 1; i <= nnodes; i++) {
        Node *n = d->giveNode(l2gmap.at(i));
        this->getNodalVariableFromPrimaryField(nodalvalues, n, tStep, ut, region);
        for (int j =1; j<=nodalvalues.giveSize(); j++) {
            values.push_back(nodalvalues.at(j));
        }
    }
    printf ("Primary Values of %s:", __UnknownTypeToString(ut));
    for (auto const &i: values) {
        printf(" %e", i );
    }
    printf("\n");

    sv.setValues(values);
    // write data to VMAP file
    std::string groupname = this->getGroupName(d, region);
    int partid=d->giveNumber();
    std::string vargroup = file.createVariablesGroup(tStep->giveNumber(), partid);
    file.setVariableStateInformation(tStep->giveNumber(), d->getDomainTypeString(), tStep->giveTargetTime(), tStep->giveTargetTime(), tStep->giveNumber());
    file.writeVariable(vargroup, sv);
}

void 
VMAPExportModule::writeInternalVariable(VMAP::VMAPFile& file, Domain*d, int region, InternalStateType ut, TimeStep* tStep)
{
    std::vector<double> values;
    VMAP::sStateVariable sv;
    ExportRegion &er = this->exportRegions[region-1];
    IntArray &cells = er.getRegionCells();

    sv.setIdentifier(int(ut)); // @todo: identifier should be unique across UnknownType and IntVarType
    sv.setVariableName(__InternalStateTypeToString(ut));
    sv.setTimeValue(tStep->giveTargetTime());
    InternalStateValueType isvt = giveInternalStateValueType(ut);
    int varsize = 1;// variable size
    // determine variable dimension and oofem-vmap component mapping fro tensors
    IntArray compMap;
    VMAP::sStateVariable::eVariableDimension vd;
    if (isvt == ISVT_SCALAR) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_SCALAR;
    } else if (isvt == ISVT_VECTOR) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_VECTOR;
        varsize = 3;
    } else if (isvt == ISVT_TENSOR_S3 || isvt == ISVT_TENSOR_S3E) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_2ND_ORDER_TENSOR_SYMMETRIC;
        varsize = 6;
        compMap = {0,1,2,5,3,4};
    } else if (isvt == ISVT_TENSOR_G) {
        vd = VMAP::sStateVariable::eVariableDimension::DIMENSION_2ND_ORDER_TENSOR;
        varsize = 9;
        compMap = {0,4,8,1,5,2,3,7,6};
    } else {
        OOFEM_ERROR ("VMAPExportModule::writePrimaryVariable: unsupported InternalStateValueType")
    }
    sv.setDimension(vd);
    sv.setMultiplicity(1);
    sv.setLocation(VMAP::sStateVariable::eVariableLocation::LOCATION_INTEGRATION_POINT); // internal values defined in integration points
    sv.setEntity(VMAP::sStateVariable::eVariableEntity::ENTITY_REAL);
    
    int nelems = this->exportRegions[region-1].giveNumberOfCells();
    // another loop over elements and irules would be needed to prealocate ipvals.
    //ipval.reserve(nnodes*varsize);

    //    @todo
    FloatArray ipval;
    int size=0;
    // another loop over elements and irules would be needed to prealocate values.
    for(int i = 1; i <= nelems; i++) {
        Element *e = d->giveElement(cells.at(i));
        // loop over default integration rule of element
        IntegrationRule* irule = e->giveDefaultIntegrationRulePtr();
        size += irule->giveNumberOfIntegrationPoints() * varsize;
    }
    values.reserve(size);
    IntArray ipmap;

    for(int i = 1; i <= nelems; i++) {
        Element *e = d->giveElement(cells.at(i));
        // loop over default integration rule of element
        IntegrationRule* irule = e->giveDefaultIntegrationRulePtr();
        this->getIntegrationPointMapping (ipmap, irule->giveIntegrationDomain(), irule->giveNumberOfIntegrationPoints());
        // loop over integration points
        for (int j=0; j< irule->giveNumberOfIntegrationPoints(); j++) {
            // ip mapping between oofem->VMAP 
            GaussPoint* igp = irule->getIntegrationPoint(ipmap(j));
            if (e->giveIPValue(ipval, igp, ut, tStep) == 0) {
                ipval.resize(varsize);
                ipval.zero();
            }
            // now ve have a value, store it in correct VMAP component order
            if (vd == VMAP::sStateVariable::eVariableDimension::DIMENSION_SCALAR || vd == VMAP::sStateVariable::eVariableDimension::DIMENSION_VECTOR) {
                for (int k=0; k<varsize; k++) {
                    values.push_back(ipval[k]);
                }
            } else if (vd == VMAP::sStateVariable::eVariableDimension::DIMENSION_2ND_ORDER_TENSOR_SYMMETRIC || vd == VMAP::sStateVariable::eVariableDimension::DIMENSION_2ND_ORDER_TENSOR) {
                for (int k=0; k< varsize; k++) {
                    values.push_back(ipval.at(compMap[k]));
                }
            }
        }

    }
    printf ("Internal Values of %s:", __InternalStateTypeToString(ut));
    for (auto const &i: values) {
        printf(" %e", i );
    }
    printf("\n");

    sv.setValues(values);
    // write data to VMAP file
    std::string groupname = this->getGroupName(d, region);
    int partid=d->giveNumber();
    std::string vargroup = file.createVariablesGroup(tStep->giveNumber(), partid);
    file.setVariableStateInformation(tStep->giveNumber(), d->getDomainTypeString(), tStep->giveTargetTime(), tStep->giveTargetTime(), tStep->giveNumber());
    file.writeVariable(vargroup, sv);
}

void 
VMAPExportModule::getIntegrationPointMapping (IntArray& ipmap, integrationDomain id, int npoints)
{
    /* returns the mapping between VMAP integration point numbering and oofem integration point numbering */
    ipmap.resize(npoints);
    if (id == integrationDomain::_Line) {
        for (int i =0; i<npoints; i++) {
            ipmap(i)=i;
        }
    }
    else if (id == integrationDomain::_Triangle && npoints == 1) {
        ipmap={0};
    } else if (id == integrationDomain::_Triangle && npoints == 4) {
        ipmap={2,3,1,0};
    } else if (id == integrationDomain::_Triangle && npoints == 6) {
        ipmap = {4,5,3,0,1,2};
    } else if (id == integrationDomain::_Square && npoints == 1) {
        ipmap = {0};
    } else if (id == integrationDomain::_Square && npoints == 4) {
        ipmap={0,2,1,3};
    } else if (id == integrationDomain::_Square && npoints == 9) {
        ipmap={0, 3, 6, 4, 7, 2, 5, 8};
    } else if (id == integrationDomain::_Tetrahedra && npoints == 1) {
        ipmap={0};
    } else if (id == integrationDomain::_Tetrahedra && npoints== 4) {
        ipmap={3, 0, 1, 2};
    } else if (id == integrationDomain::_Wedge && npoints == 1) {
        ipmap={0};
    } else if (id == integrationDomain::_Cube && npoints == 1) {
        ipmap={0};
    } else if (id == integrationDomain::_Cube && npoints == 8) {
        ipmap={0,4,2,6,  1,5,3,7};
    } else {
        OOFEM_ERROR ("getIntegrationPointMapping: unsuported integrationDomain and number of integration points");
    }
}

int
VMAPExportModule::initRegionNodeNumbering(ExportRegion& piece, 
                                          Domain *domain, 
                                          //TimeStep *tStep, 
                                          Set& region)
{
    // regionG2LNodalNumbers is array with mapping from global numbering to local region numbering.
    // The i-th value contains the corresponding local region number (or zero, if global number is not in region).

    // regionL2GNodalNumbers is array with mapping from local to global numbering.
    // The i-th value contains the corresponding global node number.


    int nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int elementNode, node;
    int currOffset = 1;
    Element *element;

    IntArray &regionG2LNodalNumbers = piece.getMapG2L();
    IntArray &regionL2GNodalNumbers = piece.getMapL2G();

    regionG2LNodalNumbers.resize(nnodes);
    regionG2LNodalNumbers.zero();
    int regionDofMans = 0;
    int regionSingleCells = 0;

    const IntArray& elements = region.giveElementList();
    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        int ielem = elements.at(ie);
        element = domain->giveElement(ielem);

        if ( this->isElementComposite(element) ) {
            continue;                                    // composite cells exported individually
        }
/*
        if ( !element->isActivated(tStep) ) {                    //skip inactivated elements
            continue;
        }

        //skip materials with casting time > current time
        if ( !element->isCast(tStep) ) {
            continue;
        }
*/
        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        regionSingleCells++;
        elemNodes = element->giveNumberOfNodes();
        //  elemSides = element->giveNumberOfSides();

        // determine local region node numbering
        for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            node = element->giveNode(elementNode)->giveNumber();
            if ( regionG2LNodalNumbers.at(node) == 0 ) { // assign new number
                // mark for assignment. This is done later, as it allows to preserve
                // natural node numbering.
                //
                regionG2LNodalNumbers.at(node) = 1;
                regionDofMans++;
            }
        }
    }

    regionL2GNodalNumbers.resize(regionDofMans);

    for ( int i = 1; i <= nnodes; i++ ) {
        if ( regionG2LNodalNumbers.at(i) ) {
            regionG2LNodalNumbers.at(i) = currOffset++;
            regionL2GNodalNumbers.at(regionG2LNodalNumbers.at(i) ) = i;
        }
    }

    piece.setNumberOfCells(regionSingleCells);
    piece.setNumberOfNodes(regionDofMans);

    return 1;
}

void
VMAPExportModule::setupExportRegion(ExportRegion &exportRegion, 
                                    //TimeStep *tStep, 
                                    Set &region)
{
    Domain *d  = emodel->giveDomain(1);
    Element *elem;

    // output nodes Region By Region
   
    
    // Assemble local->global and global->local region map and get number of
    // single cells to process, the composite cells exported individually.
    this->initRegionNodeNumbering(exportRegion, d, region);
    const IntArray& mapG2L = exportRegion.getMapG2L();
    const IntArray& mapL2G = exportRegion.getMapL2G();
    const int numNodes = exportRegion.giveNumberOfNodes();
    const int numRegionEl = exportRegion.giveNumberOfCells();

    if ( numNodes > 0 && numRegionEl > 0 ) {
        // Setup region nodes
        exportRegion.setNumberOfNodes(numNodes);
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            const auto &coords = d->giveNode(mapL2G.at(inode) )->giveCoordinates();
            exportRegion.setNodeCoords(inode, coords);
        }


        //-------------------------------------------
        // Export all the cell data for the piece
        //-------------------------------------------
        IntArray cellNodes;
        exportRegion.setNumberOfCells(numRegionEl);

        int offset = 0;
        int cellNum = 0;
        IntArray elems = region.giveElementList();
        for ( int ei = 1; ei <= elems.giveSize(); ei++ ) {
            int elNum = elems.at(ei);
            elem = d->giveElement(elNum);

            // Skip elements that:
            // are inactivated or of composite type ( these are exported individually later)
            if ( this->isElementComposite(elem) /*|| !elem->isActivated(tStep)*/ ) {
                continue;
            }
/*
            //skip materials with casting time > current time
            if ( !elem->isCast(tStep) ) {
                continue;
            }
*/
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            exportRegion.getRegionCells().followedBy(elNum);

            cellNum++;

            // Set the connectivity
            this->giveElementConnectivity(cellNodes, elem);  // node numbering of the cell according to the VMAP format

            // Map from global to local node numbers for the current piece
            int numElNodes = cellNodes.giveSize();
            IntArray connectivity(numElNodes);
            for ( int i = 1; i <= numElNodes; i++ ) {
                connectivity.at(i) = mapG2L.at(cellNodes.at(i) );
            }

            exportRegion.setConnectivity(cellNum, connectivity);

            exportRegion.setCellType(cellNum, this->giveElementType(elem) );  // VMAP cell type

            offset += numElNodes;
            exportRegion.setOffset(cellNum, offset);
        }

       
    } // end of default piece for simple geometry elements
}


void
VMAPExportModule::giveElementConnectivity(IntArray &answer, Element *elem)
{
    // Gives the VMAP element connectivity for element

    Element_Geometry_Type elemGT = elem->giveGeometryType();
    IntArray nodeMapping(0);
    if ( ( elemGT == EGT_point ) ||
         ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
         ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
         ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
         ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
         ( elemGT == EGT_hexa_1 ) || ( elemGT == EGT_quad9_2 ) ||
         ( elemGT == EGT_wedge_1 ) ) 
    {

    } else if ( elemGT == EGT_hexa_27 ) {
        nodeMapping = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18, 23, 25, 26, 24, 22, 21, 27 // needs check
        };
    } else if ( elemGT == EGT_hexa_2 ) {
        nodeMapping = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18 // needs check
        };
    } else if ( elemGT == EGT_wedge_2 ) {
        nodeMapping = {
            4, 6, 5, 1, 3, 2, 12, 11, 10, 9, 8, 7, 13, 15, 14 // needs check
        };
    } else if ( elemGT == EGT_quad_1_interface ) {
        nodeMapping = {
            1, 2, 4, 3, 5, 6, 8, 7 // needs check
        };
    } else if ( elemGT == EGT_quad_21_interface ) {
        nodeMapping = {
            1, 2, 5, 4, 3, 6 // needs check
        };
    } else {
        OOFEM_ERROR("VTKXMLExportModule: unsupported element geometry type");
    }

    int nelemNodes = elem->giveNumberOfNodes();
    answer.resize(nelemNodes);
    if ( nodeMapping.giveSize() > 0 ) {
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(nodeMapping.at(i) )->giveNumber();
        }
    } else {
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber();
        }
    }
}


bool
VMAPExportModule::isElementComposite(Element *elem)
{
    return ( elem->giveGeometryType() == EGT_Composite );
}

/*
VMAP::sElementType::eElementShape
VMAPExportModule::giveElementType(Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    VMAP::sElementType::eElementShape et = VMAP::sElementType::eElementShape::SHAPE_INVALID;

    if ( elemGT == EGT_point ) {
        et = VMAP::sElementType::eElementShape::POINT;
    } else if ( elemGT == EGT_line_1 ) {
        et = VMAP::sElementType::eElementShape::LINE_2;
    } else if ( elemGT == EGT_line_2 ) {
        et = VMAP::sElementType::eElementShape::LINE_3;
    } else if ( elemGT == EGT_triangle_1 ) {
        et = VMAP::sElementType::eElementShape::TRIANGLE_3;
    } else if ( elemGT == EGT_triangle_2 ) {
        et = VMAP::sElementType::eElementShape::TRIANGLE_6;
    } else if ( elemGT == EGT_tetra_1 ) {
        et = VMAP::sElementType::eElementShape::TETRAHEDRON_4;
    } else if ( elemGT == EGT_tetra_2 ) {
        et = VMAP::sElementType::eElementShape::TETRAHEDRON_10;
    } else if ( elemGT == EGT_quad_1) {
        et = VMAP::sElementType::eElementShape::QUAD_4;
    // } else if ( elemGT == EGT_quad_21_interface ) {
    // vtkCellType = 30;
    } else if ( elemGT == EGT_quad_2 ) {
        et = VMAP::sElementType::eElementShape::QUAD_8;
    } else if ( elemGT == EGT_quad9_2 ) {
        et = VMAP::sElementType::eElementShape::QUAD_9;
    } else if ( elemGT == EGT_hexa_1  || elemGT == EGT_quad_1_interface ) {
        et = VMAP::sElementType::eElementShape::HEXAHEDRON_8;
    } else if ( elemGT == EGT_hexa_2 ) {
        et = VMAP::sElementType::eElementShape::HEXAHEDRON_20;
    } else if ( elemGT == EGT_hexa_27 ) {
        et = VMAP::sElementType::eElementShape::HEXAHEDRON_27;
    } else if ( elemGT == EGT_wedge_1 ) {
        et = VMAP::sElementType::eElementShape::WEDGE_6;
    } else if ( elemGT == EGT_wedge_2 ) {
        et = VMAP::sElementType::eElementShape::WEDGE_15;
    } else {
        OOFEM_ERROR("unsupported element geometry type on element %d", elem->giveNumber() );
    }

    return et;
}
*/
int
VMAPExportModule::giveElementType(Element *elem)
{
    // make a hash as unique combination of oofem element geon=metry type and integration
    int hash = this->getElementTypeHash(elem);
    //if (vmapelementtype.contains(hash)) {
    auto found = vmapelementtype.find(hash);
    if (found != vmapelementtype.end()) {
        return this->vmapelementtype[hash];
    } else {
        OOFEM_ERROR("VMAPExportModule::giveElementType: internal error (element type not initialized");
        return -1;
    }

}

int VMAPExportModule::giveElementType(int num) {
  return giveElementType(this->emodel->giveDomain(1)->giveElement(num));
}

std::string
VMAPExportModule::getGroupName (Domain* d, int region)  
{
    std::string groupName = "D" + std::to_string(d->giveNumber())+"R"+std::to_string(region);
    return groupName;
}

void
VMAPExportModule::getNodalVariableFromPrimaryField(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, int region)
{


    IntArray dofIDMask(3);
    dofIDMask.clear();
    
    if ( (type == DisplacementVector) || (type == ResidualForce) ) {
        dofIDMask = {
            ( int ) Undef, ( int ) Undef, ( int ) Undef
        };
        for ( Dof *dof : * dman ) {
            DofIDItem id = dof->giveDofID();
            if ( id == D_u ) {
                dofIDMask.at(1) = id;
            } else if ( id == D_v ) {
                dofIDMask.at(2) = id;
            } else if ( id == D_w ) {
                dofIDMask.at(3) = id;
            }
        }

        answer.resize(3);
    } else if ( type == VelocityVector ) {
        dofIDMask = {
            ( int ) Undef, ( int ) Undef, ( int ) Undef
        };
        for ( Dof *dof : * dman ) {
            DofIDItem id = dof->giveDofID();
            if ( id == V_u ) {
                dofIDMask.at(1) = id;
            } else if ( id == V_v ) {
                dofIDMask.at(2) = id;
            } else if ( id == V_w ) {
                dofIDMask.at(3) = id;
            }
        }

        answer.resize(3);
    } else if ( type == FluxVector || type == Humidity ) {
        dofIDMask.followedBy(C_1);
        answer.resize(1);
    } else if ( type == Temperature ) {
        dofIDMask.followedBy(T_f);
        answer.resize(1);
    } else if ( type == PressureVector ) {
        dofIDMask.followedBy(P_f);
        answer.resize(1);
    } else {
        OOFEM_ERROR("unsupported unknownType %s", __UnknownTypeToString(type) );
    }

    int size = dofIDMask.giveSize();
    answer.zero();

    for ( int j = 1; j <= size; j++ ) {
        DofIDItem id = ( DofIDItem ) dofIDMask.at(j);
        if ( id == Undef ) {
            answer.at(j) = 0.;
        } else if ( dman->hasDofID(id) ) {
            // primary variable available directly in DOF-manager
            answer.at(j) = dman->giveDofWithID(id)->giveUnknown(VM_Total, tStep);
            //mj - if incremental value needed: answer.at(j) = dman->giveDofWithID(id)->giveUnknown(VM_Incremental, tStep);
        }
    }

    InternalStateValueType valType = giveInternalStateValueType(type);
    //rotate back from nodal CS to global CS if applies
    if ( valType == ISVT_VECTOR ) { ///@todo in general, shouldn't this apply for 2nd order tensors as well? /JB
        Node *node = dynamic_cast< Node * >( dman );
        if ( node && node->hasLocalCS() ) {
            answer.rotatedWith(* node->giveLocalCoordinateTriplet(), 't');
        }
    }
}


} // end namespace oofem
