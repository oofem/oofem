#include "util.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "intarray.h"
#include "floatarray.h"


// Optional (only need the input fields defines)
#include "linearstatic.h"
#include "simplecrosssection.h"
#include "isolinearelasticmaterial.h"
#include "truss2d.h"
#include "generalboundarycondition.h"
#include "beam2d.h"
#include "constantedgeload.h"
#include "nodalload.h"
#include "structtemperatureload.h"
#include "peakfunction.h"
#include "node.h"
#include "outputmanager.h"
#include "boundarycondition.h"
#include "set.h"

using namespace oofem;

int main(int argc, char *argv[])
{
    DynamicDataReader myData;
    DynamicInputRecord *myInput;

    //Output File
    myData.setOutputFileName("beam01.out");

    //Description
    myData.setDescription("My custom beam problem");

    //Problem
    myInput = new DynamicInputRecord(_IFT_LinearStatic_Name);
    myInput->setField(0, _IFT_EngngModel_lstype);
    myInput->setField(0, _IFT_EngngModel_smtype);
    myInput->setField(3, _IFT_EngngModel_nsteps);
    myData.insertInputRecord(DataReader::IR_emodelRec, myInput);

    //Domain
    ///@todo Remove this.
    myInput = new DynamicInputRecord();
    std::string help = "2dbeam";
    // myInput->setRecordKeywordField("domain", 1);
    myInput->setField(help, _IFT_Domain_type);
    myData.insertInputRecord(DataReader::IR_domainRec, myInput);

    //Output
    ///@todo Make this a normal export module.
    myInput = new DynamicInputRecord();
    //myInput->setRecordKeywordField(_IFT_OutputManager_name, 1);
    myInput->setField(_IFT_OutputManager_Name);
    myInput->setField(_IFT_OutputManager_tstepall);
    myInput->setField(_IFT_OutputManager_dofmanall);
    myInput->setField(_IFT_OutputManager_elementall);
    myData.insertInputRecord(DataReader::IR_outManRec, myInput);

    //Components size record
    ///@todo This should probably have "domain 1" in the beginning. Makes sense.
    myInput = new DynamicInputRecord();
    myInput->setField(6, _IFT_Domain_ndofman);
    myInput->setField(5, _IFT_Domain_nelem);
    myInput->setField(1, _IFT_Domain_ncrosssect);
    myInput->setField(1, _IFT_Domain_nmat);
    myInput->setField(6, _IFT_Domain_nbc);
    myInput->setField(0, _IFT_Domain_nic);
    myInput->setField(3, _IFT_Domain_nfunct);
    myInput->setField(4, _IFT_Domain_nset);
    myData.insertInputRecord(DataReader::IR_domainCompRec, myInput);

    //Nodes
    FloatArray coord(3);
    IntArray bcs(3);
    IntArray loads(1);

    coord.setValues(3, 0.0, 0.0, 0.0);
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(1, _IFT_Node_Name, coord));

    coord.setValues(3, 2.4, 0.0, 0.0);
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(2, _IFT_Node_Name, coord));

    coord.setValues(3, 3.8, 0.0, 0.0);
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(3, _IFT_Node_Name, coord));

    coord.setValues(3, 5.8, 0.0, 1.5);
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(4, _IFT_Node_Name, coord));

    coord.setValues(3, 7.8, 0.0, 3.0);
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(5, _IFT_Node_Name, coord));

    coord.setValues(3, 2.4, 0.0, 3.0);
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(6, _IFT_Node_Name, coord));

    //Elements
    IntArray nodes(2);
    IntArray dtc(1);
    DynamicInputRecord* beam;

    nodes.setValues(2,1,2);
    beam = CreateElementIR(1, _IFT_Beam2d_Name, nodes);
    loads.setValues(2,3,1);
    beam->setField(loads, _IFT_Element_boundaryload);
    loads.setValues(1,5);
    beam->setField(loads, _IFT_Element_bodyload);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    nodes.setValues(2,2,3);
    beam = CreateElementIR(2, _IFT_Beam2d_Name, nodes);
    dtc.setValues(1,6);
    beam->setField(dtc, _IFT_Beam2d_dofstocondense);
    beam->setField(loads, _IFT_Element_bodyload);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    nodes.setValues(2,3,4);
    beam = CreateElementIR(3, _IFT_Beam2d_Name, nodes);
    dtc.setValues(1, 3);
    beam->setField(dtc, _IFT_Beam2d_dofstocondense);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    nodes.setValues(2,4,5);
    beam = CreateElementIR(4, _IFT_Beam2d_Name, nodes);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    nodes.setValues(2,6,2);
    beam = CreateElementIR(5, _IFT_Beam2d_Name, nodes);
    dtc.setValues(1, 6);
    beam->setField(dtc, _IFT_Beam2d_dofstocondense);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    //CrossSection
    myInput = new DynamicInputRecord(_IFT_SimpleCrossSection_Name, 1);
    myInput->setField(1, _IFT_CrossSection_SetNumber);
    myInput->setField(1, _IFT_SimpleCrossSection_MaterialNumber);
    myInput->setField(0.162, _IFT_SimpleCrossSection_area);
    myInput->setField(0.0039366, _IFT_SimpleCrossSection_iy);
    myInput->setField(1.e18, _IFT_SimpleCrossSection_shearcoeff);
    myInput->setField(0.54, _IFT_SimpleCrossSection_thick);
    myData.insertInputRecord(DataReader::IR_crosssectRec, myInput);

    //Material
    myInput = new DynamicInputRecord(_IFT_IsotropicLinearElasticMaterial_Name, 1);
    myInput->setField(1.0, _IFT_Material_density);
    myInput->setField(30.e6, _IFT_IsotropicLinearElasticMaterial_e);
    myInput->setField(0.2, _IFT_IsotropicLinearElasticMaterial_n);
    myInput->setField(1.2e-5, _IFT_IsotropicLinearElasticMaterial_talpha);
    myData.insertInputRecord(DataReader::IR_matRec, myInput);

    //Boundary Conditions
    FloatArray bcVals;
    IntArray dofs;

    myInput = new DynamicInputRecord(_IFT_BoundaryCondition_Name, 1);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    bcVals.setValues(1, 0.);
    dofs.setValues(1, D_w);
    myInput->setField(bcVals, _IFT_BoundaryCondition_values);
    myInput->setField(dofs, _IFT_GeneralBoundaryCondition_dofs);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_BoundaryCondition_Name, 2);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    bcVals.setValues(1, 0.);
    dofs.setValues(1, R_v);
    myInput->setField(bcVals, _IFT_BoundaryCondition_values);
    myInput->setField(dofs, _IFT_GeneralBoundaryCondition_dofs);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);
    
    myInput = new DynamicInputRecord(_IFT_BoundaryCondition_Name, 6);
    myInput->setField(2, _IFT_GeneralBoundaryCondition_timeFunct);
    bcVals.setValues(3, 0., 0., -0.006e-3);
    dofs.setValues(3, D_u, D_w, R_v);
    myInput->setField(bcVals, _IFT_BoundaryCondition_values);
    myInput->setField(dofs, _IFT_GeneralBoundaryCondition_dofs);
    myInput->setField(4, _IFT_GeneralBoundaryCondition_set);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);


    //Loads
    FloatArray comps(3);    

    myInput = new DynamicInputRecord(_IFT_ConstantEdgeLoad_Name, 3);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    comps.setValues(3, 0.0, 10.0, 0.0);
    myInput->setField(comps, _IFT_Load_components);
    myInput->setField(3, _IFT_BoundaryLoad_loadtype);
    myInput->setField(3, _IFT_BoundaryLoad_ndofs);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_NodalLoad_Name, 4);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    comps.setValues(3, -18.0, 24.0, 0.0);
    myInput->setField(comps, _IFT_Load_components);
    myInput->setField(3, _IFT_GeneralBoundaryCondition_set);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_StructuralTemperatureLoad_Name, 5);
    myInput->setField(3, _IFT_GeneralBoundaryCondition_timeFunct);
    comps.setValues(2, 30.0, -20.0);
    myInput->setField(comps, _IFT_Load_components);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    //Functions
    myInput = new DynamicInputRecord(_IFT_PeakFunction_Name, 1);
    myInput->setField(1.0, _IFT_PeakFunction_t);
    myInput->setField(1.0, _IFT_PeakFunction_ft);
    myData.insertInputRecord(DataReader::IR_funcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_PeakFunction_Name, 2);
    myInput->setField(2.0, _IFT_PeakFunction_t);
    myInput->setField(1.0, _IFT_PeakFunction_ft);
    myData.insertInputRecord(DataReader::IR_funcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_PeakFunction_Name, 3);
    myInput->setField(3.0, _IFT_PeakFunction_t);
    myInput->setField(1.0, _IFT_PeakFunction_ft);
    myData.insertInputRecord(DataReader::IR_funcRec, myInput);
    
    //Sets
    IntArray vals;

    myInput = new DynamicInputRecord(_IFT_Set_Name, 1);
    vals.setValues(5, 1, 2, 3, 4, 5);
    myInput->setField(vals, _IFT_Set_elements);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    myInput = new DynamicInputRecord(_IFT_Set_Name, 2);
    vals.setValues(2, 1, 5);
    myInput->setField(vals, _IFT_Set_nodes);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    myInput = new DynamicInputRecord(_IFT_Set_Name, 3);
    vals.setValues(1, 3);
    myInput->setField(vals, _IFT_Set_nodes);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    myInput = new DynamicInputRecord(_IFT_Set_Name, 4);
    vals.setValues(1, 6);
    myInput->setField(vals, _IFT_Set_nodes);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    // Writing to file (to verify, and for backups)
    myData.writeToFile("beam01.in");

    EngngModel *em = InstanciateProblem(&myData, _processor, 0);
    em->solveYourself();
    myData.finish();
}
