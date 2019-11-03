#include "util.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "intarray.h"
#include "floatarray.h"
#include "timer.h"

// Optional (only need the input fields defines)
#include "tm/stationarytransportproblem.h"
#include "nrsolver.h"
#include "tm/simpletransportcrosssection.h"
#include "tm/isoheatmat.h"
#include "tm/brick1_ht.h"
#include "tm/transportgradientdirichlet.h"
#include "tm/transportgradientneumann.h"
#include "tm/transportgradientperiodic.h"
#include "modulemanager.h"
#include "exportmodule.h"
#include "vtkxmlexportmodule.h"
#include "generalboundarycondition.h"
#include "constantfunction.h"
#include "node.h"
#include "outputmanager.h"
#include "boundarycondition.h"
#include "set.h"

#include "tm/transportgradientneumann.h"
#include "tm/transportgradientdirichlet.h"
#include "tm/transportgradientperiodic.h"

#include <random>
#include <fstream>

#include <H5Cpp.h>

#ifdef __PETSC_MODULE
#include <petsc.h>
#endif

using namespace oofem;
using namespace H5;


class OOFEM_EXPORT BasicInputRecord : public InputRecord
{
public:
    std::unique_ptr<InputRecord> clone() const override { return std::unique_ptr<InputRecord>(); }
    std :: string giveRecordAsString() const override { return ""; }
    void giveField(int &answer, InputFieldType id) override { }
    void giveField(double &answer, InputFieldType id) override { }
    void giveField(bool &answer, InputFieldType id) override { }
    void giveField(std :: string &answer, InputFieldType id) override { }
    void giveField(FloatArray &answer, InputFieldType id) override { }
    void giveField(IntArray &answer, InputFieldType id) override { }
    void giveField(FloatMatrix &answer, InputFieldType id) override { }
    void giveField(std :: vector< std :: string > &answer, InputFieldType id) override { }
    void giveField(Dictionary &answer, InputFieldType id) override { }
    void giveField(std :: list< Range > &answer, InputFieldType id) override { }
    void giveField(ScalarFunction &function, InputFieldType id) override { }
    void printYourself() override {}
    void finish(bool wrn = true) override {}
};

class OOFEM_EXPORT BasicNodeInputRecord : public BasicInputRecord
{
protected:
    int recordNumber;
    FloatArray coords;

public:
    BasicNodeInputRecord(int num, FloatArray coords): recordNumber(num), coords(std::move(coords)) { }
    virtual ~BasicNodeInputRecord() {}

    void giveRecordKeywordField(std :: string &answer, int &value) override { answer = "node"; value = recordNumber; }
    void giveRecordKeywordField(std :: string &answer) override { answer = "node";  }
    void giveField(FloatArray &answer, InputFieldType id) override {
        if (std::string(id) == _IFT_Node_coords) answer = coords;
        else throw MissingKeywordInputException(this*, id, recordNumber);
    }

    bool hasField(InputFieldType id) override { return std::string(id) == _IFT_Node_coords; }
};


class OOFEM_EXPORT BasicElementInputRecord : public BasicInputRecord
{
protected:
    int recordNumber;
    IntArray enodes;

public:
    BasicElementInputRecord(int num, IntArray enodes): recordNumber(num), enodes(std::move(enodes)) { }
    virtual ~BasicElementInputRecord() {}

    void giveRecordKeywordField(std :: string &answer, int &value) override { answer = _IFT_Brick1_ht_Name; value = recordNumber; }
    void giveRecordKeywordField(std :: string &answer) override { answer = _IFT_Brick1_ht_Name; }
    void giveField(IntArray &answer, InputFieldType id) override { 
        if (std::string(id) == _IFT_Element_nodes) answer = enodes;
        else throw MissingKeywordInputException(this*, id, recordNumber);
    }

    bool hasField(InputFieldType id) override { return std::string(id) == _IFT_Element_nodes; }
};


class OOFEM_EXPORT DeactivatedElementInputRecord : public BasicElementInputRecord
{
public:
    DeactivatedElementInputRecord(int num, IntArray enodes): BasicElementInputRecord(num, enodes) { }
    virtual ~DeactivatedElementInputRecord() {}

    virtual void giveField(IntArray &answer, InputFieldType id) {
        if (std::string(id) == _IFT_Element_nodes) answer = enodes;
        else throw MissingKeywordInputException(this*, id, recordNumber);
    }
    virtual void giveField(int &answer, InputFieldType id) { 
        if (std::string(id) == _IFT_Element_activityTimeFunction) answer = 2;
        else if (std::string(id) == _IFT_Element_nip) answer = 0;
        else throw MissingKeywordInputException(this*, id, recordNumber);
    }
    virtual bool hasField(InputFieldType id) { return std::string(id) == _IFT_Element_activityTimeFunction || 
        std::string(id) == _IFT_Element_nip || std::string(id) == _IFT_Element_nodes; }
};



inline int nC(int nX, int nY, int nZ, const IntArray &n)
{
    return nX + nY*n[0] + nZ*n[1]*n[0];
}

FloatArray read_float_dataset(const DataSet &d)
{
    auto dspace = d.getSpace();
    hsize_t dims[2];
    auto rank = dspace.getSimpleExtentDims(dims, nullptr);
    if (rank != 1) { OOFEM_ERROR("Dataset is not 1D"); }
    FloatArray x((int)dims[0]);
    d.read(x.givePointer(), PredType::NATIVE_DOUBLE, dspace);
    return x;
}


IntArray read_int_dataset(const DataSet &d)
{
    auto dspace = d.getSpace();
    hsize_t dims[2];
    auto rank = dspace.getSimpleExtentDims(dims, nullptr);
    if (rank != 1) { OOFEM_ERROR("Dataset is not 1D"); }
    IntArray x((int)dims[0]);
    d.read(x.givePointer(), PredType::NATIVE_INT, dspace);
    return x;
}


int main(int argc, char *argv[])
{
#ifdef __PARALLEL_MODE
 #ifdef __USE_MPI
    int rank;
    MPI_Init(& argc, & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
 #endif
#endif
#ifdef __PETSC_MODULE
    PetscInitialize(& argc, & argv, PETSC_NULL, PETSC_NULL);
#endif

    std :: string filename = argv[1];
    std :: string name = argv[2];
    std :: string bc = argv[3];
    int tangentProblem = atoi(argv[4]);

    printf("Input: %s\n", filename.c_str());
    printf("Ouput basename: %s\n", name.c_str());
    printf("BC type: %s\n", bc.c_str());
    printf("tangent problem: %d\n", tangentProblem);
    // elname, w = _Hex21Stokes_Name, 3
    // elname, w = _Q27Space_Name, 3
    // elname, w = _IFT_LSpace_Name, 2
    // elname, w = _IFT_QBrick1_ht_Name, 3
    std :: string elname = _IFT_Brick1_ht_Name;

    std :: unique_ptr< EngngModel > em;

    {
        Timer timer;
        timer.startTimer();
        DynamicDataReader myData("dream3d_analysis");

        std::unique_ptr<DynamicInputRecord> myInput;
        myData.setOutputFileName(name + ".out");
        myData.setDescription("Internally generated hex grid from Dream3D file");

        // Read the file with all inclusions:
        H5File file(filename, H5F_ACC_RDONLY );
        Group voxel_data = file.openGroup("DataContainers/VoxelDataContainer");
        Group simpl_geometry = voxel_data.openGroup("_SIMPL_GEOMETRY");
        Group cell_data = voxel_data.openGroup("CellData");
        IntArray nelem = read_int_dataset(simpl_geometry.openDataSet("DIMENSIONS"));
        IntArray n = nelem; n.add(1);
        FloatArray spacing = read_float_dataset(simpl_geometry.openDataSet("SPACING"));
        IntArray phases = read_int_dataset(cell_data.openDataSet("Phases"));
        //IntArray grain_ids = read_int_data(cell_data.openDataSet("GrainIds"));
        //std::vector< FloatMatrix > euler_angles = read_matrix_data(cell_data.openDataSet("EulerAngles"));
        file.close();

        FloatArray rveSize = {spacing[0]*nelem[0], spacing[1]*nelem[1], spacing[2]*nelem[2]}; 

        printf("%s: Number of voxels = %d\n", name.c_str(), nelem[0] * nelem[1] * nelem[2]);
        rveSize.printYourself("rve size");

        //Problem
        myInput = std::make_unique<DynamicInputRecord>(_IFT_StationaryTransportProblem_Name);
        myInput->setField(1, _IFT_EngngModel_nsteps);
        myInput->setField(1e-6, _IFT_NRSolver_rtolf);
        myInput->setField(3, _IFT_EngngModel_lstype);
        myInput->setField(7, _IFT_EngngModel_smtype);
        myInput->setField(1, _IFT_ModuleManager_nmodules);
        myInput->setField(_IFT_EngngModel_suppressOutput);
        myInput->setField(_IFT_StationaryTransportProblem_keepTangent);
        myData.insertInputRecord(DataReader::IR_emodelRec, std::move(myInput));

        // VTKXML tstep_all domain_all primvars 1 6 cellvars 3 103 56 41'
        myInput = std::make_unique<DynamicInputRecord>(_IFT_VTKXMLExportModule_Name);
        myInput->setField(_IFT_ExportModule_tstepall);
        myInput->setField(_IFT_ExportModule_domainall);
        myInput->setField(IntArray{6}, _IFT_VTKXMLExportModule_primvars);
        myInput->setField(IntArray{103, 56, 41}, _IFT_VTKXMLExportModule_cellvars);
        myData.insertInputRecord(DataReader::DataReader::IR_expModuleRec, std::move(myInput));

        //Domain
        ///@todo Remove this.
        myInput = std::make_unique<DynamicInputRecord>();
        std::string help = "3d";
        // myInput->setRecordKeywordField("domain", 1);
        myInput->setField(help, _IFT_Domain_type);
        myData.insertInputRecord(DataReader::IR_domainRec, std::move(myInput));

        //Output
        ///@todo Remove this.
        myInput = std::make_unique<DynamicInputRecord>();
        //myInput->setRecordKeywordField(_IFT_OutputManager_name, 1);
        myInput->setField(_IFT_OutputManager_Name);
        myData.insertInputRecord(DataReader::IR_outManRec, std::move(myInput));

        //Components size record
        myInput = std::make_unique<DynamicInputRecord>();
        myInput->setField(n[0]*n[1]*n[2], _IFT_Domain_ndofman);
        myInput->setField(nelem[0]*nelem[1]*nelem[2], _IFT_Domain_nelem);
        myInput->setField(2, _IFT_Domain_ncrosssect);
        myInput->setField(2, _IFT_Domain_nmat);
        myInput->setField(2, _IFT_Domain_nbc);
        myInput->setField(0, _IFT_Domain_nic);
        myInput->setField(2, _IFT_Domain_nfunct);
        myInput->setField(12, _IFT_Domain_nset);
        myInput->setField(3, _IFT_Domain_numberOfSpatialDimensions);
        myData.insertInputRecord(DataReader::IR_domainCompRec, std::move(myInput));

        //Nodes
        for (int nz = 0; nz < n[2]; ++nz) {
            for (int ny = 0; ny < n[1]; ++ny) {
                for (int nx = 0; nx < n[0]; ++nx) {
                    int node = nC(nx, ny, nz, n) + 1;
#if 1
                    myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(node, _IFT_Node_Name, 
                        {nx * spacing[0] - rveSize[0] * 0.5,
                        ny * spacing[1] - rveSize[1] * 0.5,
                        nz * spacing[2] - rveSize[2] * 0.5}));
#else
                    auto x = std::make_unique<BasicNodeInputRecord>(node,
                        {nx * spacing[0] - rveSize[0] * 0.5,
                        ny * spacing[1] - rveSize[1] * 0.5,
                        nz * spacing[2] - rveSize[2] * 0.5});

                    myData.insertInputRecord(DataReader::IR_dofmanRec, std::move(x));
#endif
                }
            }
        }

        IntArray co_nodes(n[0]*n[1]*n[2]);

        // It would be exceedingly painful to renumber everything based on the WC elements, so lets just deactive them.
        // This way, no special considerations are needed in the BC code either.
        //Elements
        IntArray enode_base = {n[0]*n[1]+n[2], n[0]*n[1]+n[0]+1, n[0]*n[1] + 1, n[1]*n[0], n[0], n[0]+1, 1, 0};
        for (int eZ = 0; eZ < nelem[2]; ++eZ) {
            for (int eY = 0; eY < nelem[1]; ++eY) {
                for (int eX = 0; eX < nelem[0]; ++eX) {
                    int e = nC(eX, eY, eZ, nelem) + 1;
                    //printf("e = %d\n", e);
                    IntArray enodes = enode_base;
                    enodes.add(nC(eX, eY, eZ, n) + 1);
#if 1
                    auto eir = CreateElementIR(e, elname.c_str(), enodes);
                    if ( phases.at(e) != 1 ) { // Not Co phase, so lets deactivate
                        eir->setField(2, _IFT_Element_activityTimeFunction);
                    } else {
                        for ( auto node : enodes )
                            co_nodes.at(node) = 1;
                    }
                    myData.insertInputRecord(DataReader::IR_elemRec, std::move(eir));
#else
                    if ( phases.at(e) != 1 ) { // Not Co phase, so lets deactivate
                        myData.insertInputRecord(DataReader::IR_elemRec, std::make_unique<DeactivatedElementInputRecord>(e, enodes));
                        //eir->setField(2, _IFT_Element_activityTimeFunction);
                    } else {
                        myData.insertInputRecord(DataReader::IR_elemRec, std::make_unique<BasicElementInputRecord>(e, enodes));
                        for ( auto node : enodes )
                            co_nodes.at(node) = 1;
                    }
#endif
                }
            }
        }

        // Special concern for b.c. If Dirichlet, then we must allow the boundary nodes to be free.
        // For periodic b.c. then we must mark the master surface (+z,+y,+x side) to Co nodes, as well as free up the slaves (they follow the master anyway).
        // Here, only periodic is supported:
        // X-Y plane:
        for (int ny = 0; ny < n[1]; ++ny)
            for (int nx = 0; nx < n[0]; ++nx)
                if ( co_nodes[nC(nx, ny, n[2]-1, n)] != 0 )
                    co_nodes[nC(nx, ny, 0, n)] = 1;
        // X-Z plane:
        for (int nz = 0; nz < n[2]; ++nz)
            for (int nx = 0; nx < n[0]; ++nx)
                if ( co_nodes[nC(nx, n[1]-1, nz, n)] != 0 )
                    co_nodes[nC(nx, 0, nz, n)] = 1;
        // Y-Z plane:
        for (int nz = 0; nz < n[2]; ++nz)
            for (int ny = 0; ny < n[1]; ++ny)
                if ( co_nodes[nC(n[0]-1, ny, nz, n)] != 0 )
                    co_nodes[nC(0, ny, nz, n)] = 1;

        // Free all slaves:
        for (int ny = 0; ny < n[1]; ++ny)
            for (int nx = 0; nx < n[0]; ++nx)
                co_nodes[nC(nx, ny, n[2]-1, n)] = 1;
        for (int nz = 0; nz < n[2]; ++nz)
            for (int nx = 0; nx < n[0]; ++nx)
                co_nodes[nC(nx, n[1]-1, nz, n)] = 1;
        for (int nz = 0; nz < n[2]; ++nz)
            for (int ny = 0; ny < n[1]; ++ny)
                co_nodes[nC(n[0]-1, ny, nz, n)] = 1;

        //Sets
        IntArray xp;  // x+
        xp.preallocate(2*nelem[2]*nelem[1]);
        for (int ez = 0; ez < nelem[2]; ++ez)
            for (int ey = 0; ey < nelem[1]; ++ey)
                xp.followedBy({nC(nelem[0]-1, ey, ez, nelem)+1, 4});

        IntArray xm;  // x-
        xm.preallocate(2*nelem[2]*nelem[1]);
        for (int ez = 0; ez < nelem[2]; ++ez)
            for (int ey = 0; ey < nelem[1]; ++ey)
                xm.followedBy({nC(0, ey, ez, nelem)+1, 6});

        IntArray yp;  // y+
        yp.preallocate(2*nelem[2]*nelem[0]);
        for (int ez = 0; ez < nelem[2]; ++ez)
            for (int ex = 0; ex < nelem[0]; ++ex)
                yp.followedBy({nC(ex, nelem[1]-1, ez, nelem)+1, 3});

        IntArray ym;  // y-
        ym.preallocate(2*nelem[2]*nelem[0]);
        for (int ez = 0; ez < nelem[2]; ++ez)
            for (int ex = 0; ex < nelem[0]; ++ex)
                ym.followedBy({nC(ex, 0, ez, nelem)+1, 5});

        IntArray zp;  // z+
        zp.preallocate(2*nelem[1]*nelem[0]);
        for (int ey = 0; ey < nelem[1]; ++ey)
            for (int ex = 0; ex < nelem[0]; ++ex)
                zp.followedBy({nC(ex, ey, nelem[2]-1, nelem)+1, 1});

        IntArray zm;  // z-
        zm.preallocate(2*nelem[1]*nelem[0]);
        for (int ey = 0; ey < nelem[1]; ++ey)
            for (int ex = 0; ex < nelem[0]; ++ex)
                zm.followedBy({nC(ex, ey, 0, nelem)+1, 2});

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 1);
        myInput->setField(xm, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 2);
        myInput->setField(ym, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 3);
        myInput->setField(zm, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 4);
        myInput->setField(xp, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 5);
        myInput->setField(yp, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 6);
        myInput->setField(zp, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        // Set of WC only nodes which gets prescribed as to not cause trouble!
        IntArray wc_nodes_index;
        co_nodes.add(-1); wc_nodes_index.findNonzeros(co_nodes);
        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 7);
        myInput->setField(wc_nodes_index, _IFT_Set_nodes);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        IntArray totm, totp, tot;
        totm.followedBy(xm);
        totm.followedBy(ym);
        totm.followedBy(zm);

        totp.followedBy(xp);
        totp.followedBy(yp);
        totp.followedBy(zp);

        tot.followedBy(totm);
        tot.followedBy(totp);

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 8);
        myInput->setField(tot, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 9);
        myInput->setField(totp, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 10);
        myInput->setField(totm, _IFT_Set_elementBoundaries);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        IntArray tmp, co_elems, wc_elems;
        tmp = phases; tmp.add(-1); wc_elems.findNonzeros(tmp);
        tmp = phases; tmp.add(-2); co_elems.findNonzeros(tmp);

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 11);
        myInput->setField(co_elems, _IFT_Set_elements);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 12);
        myInput->setField(wc_elems, _IFT_Set_elements);
        myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_SimpleTransportCrossSection_Name, 1);
        myInput->setField(1, _IFT_SimpleTransportCrossSection_material);
        myInput->setField(11, _IFT_CrossSection_SetNumber);
        myData.insertInputRecord(DataReader::IR_crosssectRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_SimpleTransportCrossSection_Name, 2);
        myInput->setField(2, _IFT_SimpleTransportCrossSection_material);
        myInput->setField(12, _IFT_CrossSection_SetNumber);
        myData.insertInputRecord(DataReader::IR_crosssectRec, std::move(myInput));

        //Material
        myInput = std::make_unique<DynamicInputRecord>(_IFT_IsotropicHeatTransferMaterial_Name, 1);
        myInput->setField(ScalarFunction(1.0), _IFT_IsotropicHeatTransferMaterial_k);
        myInput->setField(ScalarFunction(1.0), _IFT_IsotropicHeatTransferMaterial_c);
        myInput->setField(1.0, _IFT_Material_density);
        myData.insertInputRecord(DataReader::IR_matRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_IsotropicHeatTransferMaterial_Name, 2);
        myInput->setField(ScalarFunction(0.0), _IFT_IsotropicHeatTransferMaterial_k);
        myInput->setField(ScalarFunction(1.0), _IFT_IsotropicHeatTransferMaterial_c);
        myInput->setField(1.0, _IFT_Material_density);
        myData.insertInputRecord(DataReader::IR_matRec, std::move(myInput));

        //Boundary Conditions
        if ( bc == "d" or bc == "md" ) {
            myInput = std::make_unique<DynamicInputRecord>(_IFT_TransportGradientDirichlet_Name, 1);
            myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
            myInput->setField(IntArray{10}, _IFT_GeneralBoundaryCondition_dofs);
            myInput->setField(FloatArray{0., 0., 0.}, _IFT_TransportGradientDirichlet_centerCoords);
            myInput->setField(FloatArray{0., 0., 1.}, _IFT_TransportGradientDirichlet_gradient);
            myInput->setField(8, _IFT_GeneralBoundaryCondition_set);
            if ( bc == "md" ) {
                myInput->setField(IntArray{1, 2, 3, 4, 5, 6}, _IFT_TransportGradientDirichlet_surfSets);
                myInput->setField(_IFT_TransportGradientDirichlet_tractionControl);
            }
            myData.insertInputRecord(DataReader::IR_bcRec, std::move(myInput));
        } else if ( bc == "n" or bc == "mn" ) {
            myInput = std::make_unique<DynamicInputRecord>(_IFT_TransportGradientNeumann_Name, 1);
            myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
            myInput->setField(IntArray{10}, _IFT_GeneralBoundaryCondition_dofs);
            myInput->setField(FloatArray{0., 0., 0.}, _IFT_TransportGradientNeumann_centerCoords);
            myInput->setField(FloatArray{0., 0., 1.}, _IFT_TransportGradientNeumann_gradient);

            myInput->setField(IntArray{1, 2, 3, 4, 5, 6}, _IFT_TransportGradientNeumann_surfSets);
            if ( bc == "mn" ) {
                myInput->setField(_IFT_TransportGradientNeumann_dispControl);
            }
            myData.insertInputRecord(DataReader::IR_bcRec, std::move(myInput));
        } else if ( bc == "p" ) {
            myInput = std::make_unique<DynamicInputRecord>(_IFT_TransportGradientPeriodic_Name, 1);
            myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
            myInput->setField(IntArray{10}, _IFT_GeneralBoundaryCondition_dofs);
            myInput->setField(FloatArray{0., 0., 0.}, _IFT_TransportGradientPeriodic_centerCoords);
            myInput->setField(FloatArray{0., 0., 1.}, _IFT_TransportGradientPeriodic_gradient);

            myInput->setField(rveSize, _IFT_TransportGradientPeriodic_jump);
            myInput->setField(9, _IFT_GeneralBoundaryCondition_set);
            myInput->setField(10, _IFT_TransportGradientPeriodic_masterSet);

            myData.insertInputRecord(DataReader::IR_bcRec, std::move(myInput));
        } else {
            OOFEM_ERROR("Unrecognized bc type, must be d, md, n, mn, or p");
        }

        // Fixing WC-only nodes
        myInput = std::make_unique<DynamicInputRecord>(_IFT_BoundaryCondition_Name, 2);
        myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
        myInput->setField(FloatArray{0.}, _IFT_BoundaryCondition_values);
        myInput->setField(IntArray{T_f}, _IFT_GeneralBoundaryCondition_dofs);
        myInput->setField(7, _IFT_GeneralBoundaryCondition_set);
        myData.insertInputRecord(DataReader::IR_bcRec, std::move(myInput));

        //Functions
        myInput = std::make_unique<DynamicInputRecord>(_IFT_ConstantFunction_Name, 1);
        myInput->setField(1.0, _IFT_ConstantFunction_f);
        myData.insertInputRecord(DataReader::IR_funcRec, std::move(myInput));

        myInput = std::make_unique<DynamicInputRecord>(_IFT_ConstantFunction_Name, 2);
        myInput->setField(0.0, _IFT_ConstantFunction_f);
        myData.insertInputRecord(DataReader::IR_funcRec, std::move(myInput));

        timer.stopTimer();
        printf("Mesh generation time %.3f s\n", timer.getWtime());
        // Writing to file (to verify, and for backups)
        //myData.writeToFile((name + "_debug.in").c_str());

        printf("Initializing problem\n");
        timer.startTimer();
        em = InstanciateProblem(myData, _processor, 0);

        myData.finish();
        timer.stopTimer();
        printf("Instanciation time %.3f s\n", timer.getWtime());
    }

    //std::cin.get();

    Timer timer;
    timer.startTimer();
    printf("Starting analysis\n");
    if ( !tangentProblem ) {
        em->solveYourself();
    } else {
        printf("Solving tangent problem\n");
        TimeStep *tStep = em->giveNextStep();
        FloatMatrix tangent;
        if ( dynamic_cast< TransportGradientNeumann* >( em->giveDomain(1)->giveBc(1) ) ) {
            dynamic_cast< TransportGradientNeumann* >( em->giveDomain(1)->giveBc(1) )->computeTangent(tangent, tStep);
        } else if ( dynamic_cast< TransportGradientDirichlet* >( em->giveDomain(1)->giveBc(1) ) ) {
            dynamic_cast< TransportGradientDirichlet* >( em->giveDomain(1)->giveBc(1) )->computeTangent(tangent, tStep);
        } else if ( dynamic_cast< TransportGradientPeriodic* >( em->giveDomain(1)->giveBc(1) ) ) {
            dynamic_cast< TransportGradientPeriodic* >( em->giveDomain(1)->giveBc(1) )->computeTangent(tangent, tStep);
        }
        tangent.printYourself("tangent");

        std :: ofstream fout(em->giveOutputBaseFileName() + ".data", std :: ios :: out);
        fout.setf(std::ios::scientific);
        fout.precision(6);
        for ( int i = 0; i < tangent.giveNumberOfRows(); ++i ) {
            for ( int j = 0; j < tangent.giveNumberOfColumns(); ++j ) {
                fout << tangent(i, j) << " ";
            }
            fout << "\n";
        }
        fout.close();
    }

    timer.stopTimer();
    printf("Compute time %.3f s\n", timer.getWtime());

#ifdef __PETSC_MODULE
    PetscFinalize();
#endif
#ifdef __USE_MPI
    MPI_Finalize();
#endif
}
