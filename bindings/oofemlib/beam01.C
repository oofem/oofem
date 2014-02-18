#include "linearstatic.h"
#include "simplecrosssection.h"
#include "isolinearelasticmaterial.h"
#include "truss2d.h"
#include "generalboundarycondition.h"

#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "beam2d.h"
#include "constantedgeload.h"
#include "nodalload.h"
#include "structtemperatureload.h"
#include "peakfunction.h"

#include "node.h"
#include "datareader.h"
#include "outputmanager.h"
#include "boundarycondition.h"
#include "util.h"

using namespace oofem;

DynamicInputRecord* CreateNode(int i, const FloatArray& coord, const IntArray& bcs, const IntArray& loads = IntArray(0))
{
     DynamicInputRecord* result = new DynamicInputRecord();
     result->setRecordKeywordField(_IFT_Node_Name, i);
     result->setField(coord, _IFT_Node_coords);
     result->setField(bcs, _IFT_DofManager_bc);
     result->setField(loads, _IFT_DofManager_load);
     return result;
}

DynamicInputRecord* CreateBeam(int i, const IntArray& nodes, int mat, int cs)
{
     DynamicInputRecord* result = new DynamicInputRecord();
     result->setRecordKeywordField(_IFT_Beam2d_Name, i);
     result->setField(nodes, _IFT_Element_nodes);
     result->setField(mat, _IFT_Element_mat);
     result->setField(cs, _IFT_Element_crosssect);
     return result;
}


int main(int argc, char *argv[])
{
   DynamicDataReader myData;
    DynamicInputRecord *myInput;

    //Output File
    myData.setOutputFileName("outputfile.out");    

    //Description
    myData.setDescription("My custom problem");

    //Problem
    myInput = new DynamicInputRecord();
    myInput->setRecordKeywordField(_IFT_LinearStatic_Name, 1);
    myInput->setField(0, _IFT_EngngModel_lstype);
    myInput->setField(0, _IFT_EngngModel_smtype);
    myInput->setField(3, _IFT_EngngModel_nsteps);
    myData.insertInputRecord(DataReader::IR_emodelRec, myInput);

    //Domain
    myInput = new DynamicInputRecord();
    std::string help = "2dBeam";
    // myInput->setRecordKeywordField("domain", 1);
     myInput->setField(help, _IFT_Domain_type);
     myData.insertInputRecord(DataReader::IR_domainRec, myInput);

     //Output
     myInput = new DynamicInputRecord();
     //myInput->setRecordKeywordField(_IFT_OutputManager_name, 1);
     myInput->setField(_IFT_OutputManager_Name);
     myInput->setField(_IFT_OutputManager_tstepall);
     myInput->setField(_IFT_OutputManager_dofmanall);
     myInput->setField(_IFT_OutputManager_elementall);
     myData.insertInputRecord(DataReader::IR_outManRec, myInput);

     //Components size record
     myInput = new DynamicInputRecord();
     myInput->setField(6, _IFT_Domain_ndofman);
     myInput->setField(5, _IFT_Domain_nelem);
     myInput->setField(1, _IFT_Domain_ncrosssect);
     myInput->setField(1, _IFT_Domain_nmat);
     myInput->setField(5, _IFT_Domain_nbc);
     myInput->setField(0, _IFT_Domain_nic);
     myInput->setField(3, _IFT_Domain_nloadtimefunct);
     myData.insertInputRecord(DataReader::IR_domainCompRec, myInput);

     //Nodes
     FloatArray coord(3);
     IntArray bcs(3);
     IntArray loads(1);

     coord.setValues(3, 0.0, 0.0, 0.0);    
     bcs.setValues(3, 0, 1, 0);
     myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNode(1, coord, bcs));

     coord.setValues(3, 2.4, 0.0, 0.0);    
     bcs.setValues(3, 0, 0, 0);
     myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNode(2,coord, bcs));

     coord.setValues(3, 3.8, 0.0, 0.0);    
     bcs.setValues(3, 0, 0, 1);
     myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNode(3, coord, bcs));

     coord.setValues(3, 5.8, 0.0, 1.5);    
     bcs.setValues(3, 0, 0, 0);
     loads.setValues(1, 4);
     myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNode(4, coord, bcs, loads));

     coord.setValues(3, 7.8, 0.0, 3.0);    
     bcs.setValues(3, 0, 1, 0);
     myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNode(5, coord, bcs));

     coord.setValues(3, 2.4, 0.0, 3.0);    
     bcs.setValues(3, 1, 1, 2);
     myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNode(6, coord, bcs));

     //Elements
     IntArray nodes(2);    
     IntArray dtc(1);
     DynamicInputRecord* beam;

     nodes.setValues(2,1,2);
     beam = CreateBeam(1, nodes, 1, 1);    
     loads.setValues(2,3,1);
     beam->setField(loads, _IFT_Element_boundaryload);
     loads.setValues(1,5);
     beam->setField(loads, _IFT_Element_bodyload);
     myData.insertInputRecord(DataReader::IR_elemRec, beam);

     nodes.setValues(2,2,3);
     beam = CreateBeam(2, nodes, 1, 1);    
     dtc.setValues(1,6);
     beam->setField(dtc, _IFT_Beam2d_dofstocondense);    
     beam->setField(loads, _IFT_Element_bodyload);
     myData.insertInputRecord(DataReader::IR_elemRec, beam);

     nodes.setValues(2,3,4);
     beam = CreateBeam(3, nodes, 1, 1);    
     dtc.setValues(1, 3);
     beam->setField(dtc, _IFT_Beam2d_dofstocondense);
     myData.insertInputRecord(DataReader::IR_elemRec, beam);
    
     nodes.setValues(2,4,5);
     beam = CreateBeam(4, nodes, 1, 1);   
     myData.insertInputRecord(DataReader::IR_elemRec, beam);

     nodes.setValues(2,6,2);
     beam = CreateBeam(5, nodes, 1, 1);    
     dtc.setValues(1, 6);
     beam->setField(dtc, _IFT_Beam2d_dofstocondense);
     myData.insertInputRecord(DataReader::IR_elemRec, beam);

     //CrossSection
     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_SimpleCrossSection_Name, 1);
     myInput->setField(0.162, _IFT_SimpleCrossSection_area);
     myInput->setField(0.0039366, _IFT_SimpleCrossSection_iy);
     myInput->setField(1.e18, _IFT_SimpleCrossSection_shearcoeff);
     myInput->setField(0.54, _IFT_SimpleCrossSection_thick);
     myData.insertInputRecord(DataReader::IR_crosssectRec, myInput);

     //Material
     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_IsotropicLinearElasticMaterial_Name, 1);
     myInput->setField(1.0, _IFT_Material_density);
     myInput->setField(30.e6, _IFT_IsotropicLinearElasticMaterial_e);
     myInput->setField(0.2, _IFT_IsotropicLinearElasticMaterial_n);
     myInput->setField(1.2e-5, _IFT_IsotropicLinearElasticMaterial_talpha);
     myData.insertInputRecord(DataReader::IR_matRec, myInput);

     //Boundary Conditions
     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_BoundaryCondition_Name, 1);
     myInput->setField(1, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
     myInput->setField(0.0, _IFT_BoundaryCondition_PrescribedValue);
     myData.insertInputRecord(DataReader::IR_bcRec, myInput);

     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_BoundaryCondition_Name, 2);
     myInput->setField(2, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
     myInput->setField(-0.006e-3, _IFT_BoundaryCondition_PrescribedValue);
     myData.insertInputRecord(DataReader::IR_bcRec, myInput);

     //Loads
     FloatArray comps(3);    

     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_ConstantEdgeLoad_Name, 3);
     myInput->setField(1, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
     comps.setValues(3, 0.0, 10.0, 0.0);
     myInput->setField(comps, _IFT_Load_components);
     myInput->setField(3, _IFT_BoundaryLoad_loadtype);
     myInput->setField(3, _IFT_BoundaryLoad_ndofs);
     myData.insertInputRecord(DataReader::IR_bcRec, myInput);

     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_NodalLoad_Name, 4);
     myInput->setField(1, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
     comps.setValues(3, -18.0, 24.0, 0.0);
     myInput->setField(comps, _IFT_Load_components);
     myData.insertInputRecord(DataReader::IR_bcRec, myInput);

     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_StructuralTemperatureLoad_Name, 5);
     myInput->setField(3, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
     comps.setValues(2, 30.0, -20.0);
     myInput->setField(comps, _IFT_Load_components);
     myData.insertInputRecord(DataReader::IR_bcRec, myInput);
    
     //Functions
     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_PeakFunction_Name, 1);
     myInput->setField(1.0, _IFT_PeakFunction_t);
     myInput->setField(1.0, _IFT_PeakFunction_ft);
     myData.insertInputRecord(DataReader::IR_ltfRec, myInput);

     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_PeakFunction_Name, 2);
     myInput->setField(2.0, _IFT_PeakFunction_t);
     myInput->setField(1.0, _IFT_PeakFunction_ft);
     myData.insertInputRecord(DataReader::IR_ltfRec, myInput);

     myInput = new DynamicInputRecord();
     myInput->setRecordKeywordField(_IFT_PeakFunction_Name, 3);
     myInput->setField(3.0, _IFT_PeakFunction_t);
     myInput->setField(1.0, _IFT_PeakFunction_ft);
     myData.insertInputRecord(DataReader::IR_ltfRec, myInput);

    EngngModel *em = InstanciateProblem(&myData, _processor, 0);
    em->solveYourself();
    myData.finish();
}
