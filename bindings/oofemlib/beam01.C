#include "util.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "intarray.h"
#include "floatarray.h"


// Optional (only need the input fields defines)
#include "../sm/EngineeringModels/linearstatic.h"
#include "../sm/CrossSections/simplecrosssection.h"
#include "../sm/Materials/isolinearelasticmaterial.h"
#include "../sm/Elements/Bars/truss2d.h"
#include "generalboundarycondition.h"
#include "../sm/Elements/Beams/beam2d.h"
#include "constantedgeload.h"
#include "nodalload.h"
#include "../sm/Loads/structtemperatureload.h"
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
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(1, _IFT_Node_Name, {0.0, 0.0, 0.0}));
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(2, _IFT_Node_Name, {2.4, 0.0, 0.0}));
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(3, _IFT_Node_Name, {3.8, 0.0, 0.0}));
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(4, _IFT_Node_Name, {5.8, 0.0, 1.5}));
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(5, _IFT_Node_Name, {7.8, 0.0, 3.0}));
    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(6, _IFT_Node_Name, {2.4, 0.0, 3.0}));

    //Elements
    DynamicInputRecord* beam;

    beam = CreateElementIR(1, _IFT_Beam2d_Name, {1, 2});
    beam->setField(IntArray{3, 1}, _IFT_Element_boundaryload);
    beam->setField(IntArray{5}, _IFT_Element_bodyload);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    beam = CreateElementIR(2, _IFT_Beam2d_Name, {2, 3});
    beam->setField({6}, _IFT_Beam2d_dofstocondense);
    beam->setField(IntArray{5}, _IFT_Element_bodyload);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    beam = CreateElementIR(3, _IFT_Beam2d_Name, {3, 4});
    beam->setField(IntArray{3}, _IFT_Beam2d_dofstocondense);
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    beam = CreateElementIR(4, _IFT_Beam2d_Name, {4, 5});
    myData.insertInputRecord(DataReader::IR_elemRec, beam);

    beam = CreateElementIR(5, _IFT_Beam2d_Name, {6, 2});
    beam->setField(IntArray{6}, _IFT_Beam2d_dofstocondense);
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
    myInput = new DynamicInputRecord(_IFT_BoundaryCondition_Name, 1);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    myInput->setField(FloatArray{0.}, _IFT_BoundaryCondition_values);
    myInput->setField(IntArray{D_w}, _IFT_GeneralBoundaryCondition_dofs);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_BoundaryCondition_Name, 2);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    myInput->setField(FloatArray{0.}, _IFT_BoundaryCondition_values);
    myInput->setField(IntArray{R_v}, _IFT_GeneralBoundaryCondition_dofs);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_BoundaryCondition_Name, 6);
    myInput->setField(2, _IFT_GeneralBoundaryCondition_timeFunct);
    myInput->setField(FloatArray{0., 0., -0.006e-3}, _IFT_BoundaryCondition_values);
    myInput->setField(IntArray{D_u, D_w, R_v}, _IFT_GeneralBoundaryCondition_dofs);
    myInput->setField(4, _IFT_GeneralBoundaryCondition_set);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);


    //Loads
    myInput = new DynamicInputRecord(_IFT_ConstantEdgeLoad_Name, 3);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    myInput->setField(FloatArray{0.0, 10.0, 0.0}, _IFT_Load_components);
    myInput->setField(3, _IFT_BoundaryLoad_loadtype);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_NodalLoad_Name, 4);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    myInput->setField(FloatArray{-18.0, 24.0, 0.0}, _IFT_Load_components);
    myInput->setField(3, _IFT_GeneralBoundaryCondition_set);
    myData.insertInputRecord(DataReader::IR_bcRec, myInput);

    myInput = new DynamicInputRecord(_IFT_StructuralTemperatureLoad_Name, 5);
    myInput->setField(3, _IFT_GeneralBoundaryCondition_timeFunct);
    myInput->setField(FloatArray{30.0, -20.0}, _IFT_Load_components);
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
    myInput = new DynamicInputRecord(_IFT_Set_Name, 1);
    myInput->setField(IntArray{1, 2, 3, 4, 5}, _IFT_Set_elements);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    myInput = new DynamicInputRecord(_IFT_Set_Name, 2);
    myInput->setField(IntArray{1, 5}, _IFT_Set_nodes);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    myInput = new DynamicInputRecord(_IFT_Set_Name, 3);
    myInput->setField(IntArray{3}, _IFT_Set_nodes);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    myInput = new DynamicInputRecord(_IFT_Set_Name, 4);
    myInput->setField(IntArray{6}, _IFT_Set_nodes);
    myData.insertInputRecord(DataReader::IR_setRec, myInput);

    // Writing to file (to verify, and for backups)
    myData.writeToFile("beam01.in");

    EngngModel *em = InstanciateProblem(&myData, _processor, 0);
    em->solveYourself();
    myData.finish();
}
