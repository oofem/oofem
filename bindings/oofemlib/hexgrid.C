#include "util.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "intarray.h"
#include "floatarray.h"
#include "timer.h"

// Optional (only need the input fields defines)
#include "tm/EngineeringModels/stationarytransportproblem.h"
#include "nrsolver.h"
#include "tm/simpletransportcrosssection.h"
#include "tm/Materials/isoheatmat.h"
#include "tm/Elements/brick1_ht.h"
#include "tm/BoundaryCondition/transportgradientdirichlet.h"
#include "tm/BoundaryCondition/transportgradientneumann.h"
#include "tm/BoundaryCondition/transportgradientperiodic.h"
#include "modulemanager.h"
#include "exportmodule.h"
#include "vtkxmlexportmodule.h"
#include "generalboundarycondition.h"
#include "constantfunction.h"
#include "node.h"
#include "outputmanager.h"
#include "boundarycondition.h"
#include "set.h"

#include "tm/BoundaryCondition/transportgradientneumann.h"
#include "tm/BoundaryCondition/transportgradientdirichlet.h"
#include "tm/BoundaryCondition/transportgradientperiodic.h"

#include <random>
#include <fstream>

#ifdef __PETSC_MODULE
 #include <petsc.h>
#endif

using namespace oofem;


typedef std :: pair< FloatArray, double > inclusion;

std :: vector< inclusion > getInclusionsInBox(FloatArray corner, double rveSize, std :: vector< inclusion > &inclusions)
{
    std :: vector< inclusion > intersectedInclusions;
    for ( auto &inc : inclusions ) {
        FloatArray &ic = inc.first;
        double ir = inc.second;
        // Check if sphere and cube intersects
        double dmin = 0;
        for ( int i = 0; i < corner.giveSize(); ++i ) {
            if ( ic[i] < corner[i] )
                dmin += (ic[i] - corner[i])*(ic[i] - corner[i]);
            else if ( ic[i] > (corner[i] + rveSize) )
                dmin += (ic[i] - corner[i] - rveSize)*(ic[i] - corner[i] - rveSize);
        }
        if ( dmin <= ir*ir )
            intersectedInclusions.push_back(inc);
    }
    return intersectedInclusions;
}


int main(int argc, char *argv[])
{
#ifdef __MPI_PARALLEL_MODE
 #ifdef __USE_MPI
    int rank;
    MPI_Init(& argc, & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
    oofem_logger.setComm(MPI_COMM_WORLD);
 #endif
#endif
#ifdef __PETSC_MODULE
    PetscInitialize(& argc, & argv, PETSC_NULL, PETSC_NULL);
#endif

    Timer timer;
    std :: string inclusion_file = argv[1];
    std :: string name = argv[2];
    std :: string bc = argv[3];
    double k = atof(argv[4]);
    double rveSize = atof(argv[5]);
    int nelem = atoi(argv[6]) * rveSize;
    int sample = atoi(argv[7]);
    int tangentProblem = atoi(argv[8]);
    FloatArray rvePosition;
    
    DynamicDataReader myData("hexgrid");
    std::unique_ptr<DynamicInputRecord> myInput;

    // Read the file with all inclusions:
    std :: ifstream datafile(inclusion_file, std :: ios :: binary);
    double boxSize;
    FloatArray coord(3);
    double radius;
    int num_inclusions;
    std :: vector< inclusion > inclusions;
    datafile.read( reinterpret_cast< char* >( &boxSize ), sizeof(double));
    datafile.read( reinterpret_cast< char* >( &num_inclusions ), sizeof(int));
    inclusions.reserve(num_inclusions);
    for ( int i = 0; i <= num_inclusions; ++i ) {
        datafile.read( reinterpret_cast< char* >( coord.givePointer() ), 3 * sizeof(double));
        datafile.read( reinterpret_cast< char* >( &radius ), sizeof(double));
        //coord.printYourself("coord");
        //printf("radius = %e\n", radius);
        inclusions.push_back(std :: make_pair(coord, radius));
    }
    printf("inclusions.data:  boxSize = %e, %d inclusions\n", boxSize, num_inclusions);

    std :: default_random_engine rd(rveSize * (sample+1));
    std :: uniform_real_distribution<> dis(0, boxSize - rveSize);
    rvePosition = {dis(rd), dis(rd), dis(rd)};
    
    printf("%s: bc = %s, k = %.1e, sample = %d, rvePosition = [%.3e, %.3e, %.3e], rveSize = %.2e, nelem = %d\n", 
           name.c_str(), bc.c_str(), k, sample, rvePosition[0], rvePosition[1], rvePosition[2], rveSize, nelem);

    // elname, w = _Hex21Stokes_Name, 3
    // elname, w = _Q27Space_Name, 3
    // elname, w = _IFT_LSpace_Name, 2
    // elname, w = _IFT_QBrick1_ht_Name, 3
    std :: string elname = _IFT_Brick1_ht_Name;
    int w = 2;
    int n = nelem*(w-1) + 1;
    
    timer.startTimer();
    //Output File
    myData.setOutputFileName(name + "." + std :: to_string(sample) + ".out");

    //Description
    myData.setDescription("Internally generated hex grid");

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
    myInput = std::make_unique<DynamicInputRecord>();
    //myInput->setRecordKeywordField(_IFT_OutputManager_name, 1);
    myInput->setField(_IFT_OutputManager_Name);
    myData.insertInputRecord(DataReader::IR_outManRec, std::move(myInput));

    //Components size record
    myInput = std::make_unique<DynamicInputRecord>();
    myInput->setField(n*n*n, _IFT_Domain_ndofman);
    myInput->setField(nelem*nelem*nelem, _IFT_Domain_nelem);
    myInput->setField(2, _IFT_Domain_ncrosssect);
    myInput->setField(2, _IFT_Domain_nmat);
    myInput->setField(2, _IFT_Domain_nbc);
    myInput->setField(0, _IFT_Domain_nic);
    myInput->setField(1, _IFT_Domain_nfunct);
    myInput->setField(12, _IFT_Domain_nset);
    myInput->setField(3, _IFT_Domain_numberOfSpatialDimensions);
    myData.insertInputRecord(DataReader::IR_domainCompRec, std::move(myInput));
    
    //Nodes
    for (int nz = 0; nz < n; ++nz) {
        for (int ny = 0; ny < n; ++ny) {
            for (int nx = 0; nx < n; ++nx) {
                int node = nx + ny * n + nz * n * n + 1;
                myData.insertInputRecord(DataReader::IR_dofmanRec, CreateNodeIR(node, _IFT_Node_Name, 
                    {nx * rveSize / nelem - rveSize * 0.5,
                     ny * rveSize / nelem - rveSize * 0.5,
                     nz * rveSize / nelem - rveSize * 0.5}));
            }
        }
    }

    //Elements
    #define nC(nX, nY, nZ) (nX) + (nY)*n + (nZ)*n*n + 1
    for (int eZ = 0; eZ < nelem; ++eZ) {
        for (int eY = 0; eY < nelem; ++eY) {
            for (int eX = 0; eX < nelem; ++eX) {
                int e = eX + eY*nelem + eZ*nelem*nelem + 1;
                IntArray enodes;
                if ( w == 2 ) {
                    enodes = {
                        nC(0+eX, 0+eY, 1+eZ),
                        nC(0+eX, 1+eY, 1+eZ),
                        nC(1+eX, 1+eY, 1+eZ),
                        nC(1+eX, 0+eY, 1+eZ),
                        nC(0+eX, 0+eY, 0+eZ),
                        nC(0+eX, 1+eY, 0+eZ),
                        nC(1+eX, 1+eY, 0+eZ),
                        nC(1+eX, 0+eY, 0+eZ)};
                } else {
                    enodes = {
                        nC(0+eX*2, 0+eY*2, 2+eZ*2),
                        nC(0+eX*2, 2+eY*2, 2+eZ*2),
                        nC(2+eX*2, 2+eY*2, 2+eZ*2),
                        nC(2+eX*2, 0+eY*2, 2+eZ*2),

                        nC(0+eX*2, 0+eY*2, 0+eZ*2),
                        nC(0+eX*2, 2+eY*2, 0+eZ*2),
                        nC(2+eX*2, 2+eY*2, 0+eZ*2),
                        nC(2+eX*2, 0+eY*2, 0+eZ*2),

                        nC(0+eX*2, 1+eY*2, 2+eZ*2),
                        nC(1+eX*2, 2+eY*2, 2+eZ*2),
                        nC(2+eX*2, 1+eY*2, 2+eZ*2),
                        nC(1+eX*2, 0+eY*2, 2+eZ*2),

                        nC(0+eX*2, 1+eY*2, 0+eZ*2),
                        nC(1+eX*2, 2+eY*2, 0+eZ*2),
                        nC(2+eX*2, 1+eY*2, 0+eZ*2),
                        nC(1+eX*2, 0+eY*2, 0+eZ*2),

                        nC(0+eX*2, 0+eY*2, 1+eZ*2),
                        nC(0+eX*2, 2+eY*2, 1+eZ*2),
                        nC(2+eX*2, 2+eY*2, 1+eZ*2),
                        nC(2+eX*2, 0+eY*2, 1+eZ*2),

                        nC(1+eX*2, 1+eY*2, 2+eZ*2),
                        nC(1+eX*2, 1+eY*2, 0+eZ*2),
                        nC(0+eX*2, 1+eY*2, 1+eZ*2),
                        nC(1+eX*2, 2+eY*2, 1+eZ*2),
                        nC(2+eX*2, 1+eY*2, 1+eZ*2),
                        nC(1+eX*2, 0+eY*2, 1+eZ*2),

                        nC(1+eX*2, 1+eY*2, 1+eZ*2)};
                }

                myData.insertInputRecord(DataReader::IR_elemRec, CreateElementIR(e, elname.c_str(), enodes));
            }
        }
    }

    //Sets
    IntArray xp;  // x+
    xp.preallocate(2*nelem*nelem);
    for (int ez = 0; ez < nelem; ++ez)
        for (int ey = 0; ey < nelem; ++ey) {
            int ex = nelem-1;
            int e = ex + ey*nelem + ez*nelem*nelem + 1;
            xp.followedBy({e, 5});
        }

    IntArray xm;  // x-
    xm.preallocate(2*nelem*nelem);
    for (int ez = 0; ez < nelem; ++ez)
        for (int ey = 0; ey < nelem; ++ey) {
            int ex = 0;
            int e = ex + ey*nelem + ez*nelem*nelem + 1;
            xm.followedBy({e, 3});
        }

    IntArray yp;  // y+
    yp.preallocate(2*nelem*nelem);
    for (int ez = 0; ez < nelem; ++ez)
        for (int ex = 0; ex < nelem; ++ex) {
            int ey = nelem-1;
            int e = ex + ey*nelem + ez*nelem*nelem + 1;
            yp.followedBy({e, 4});
        }

    IntArray ym;  // y-
    ym.preallocate(2*nelem*nelem);
    for (int ez = 0; ez < nelem; ++ez)
        for (int ex = 0; ex < nelem; ++ex) {
            int ey = 0;
            int e = ex + ey*nelem + ez*nelem*nelem + 1;
            ym.followedBy({e, 6});
        }

    IntArray zp;  // z+
    zp.preallocate(2*nelem*nelem);
    for (int ey = 0; ey < nelem; ++ey)
        for (int ex = 0; ex < nelem; ++ex) {
            int ez = nelem-1;
            int e = ex + ey*nelem + ez*nelem*nelem + 1;
            zp.followedBy({e, 1});
        }

    IntArray zm;  // z-
    zm.preallocate(2*nelem*nelem);
    for (int ey = 0; ey < nelem; ++ey)
        for (int ex = 0; ex < nelem; ++ex) {
            int ez = 0;
            int e = ex + ey*nelem + ez*nelem*nelem + 1;
            zm.followedBy({e, 2});
        }

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

    myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 7);
    myInput->setField(IntArray{(n*n*n + n*n + n) / 2}, _IFT_Set_nodes);
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

    // Check for inclusions:
    auto rveInclusions = getInclusionsInBox(rvePosition, rveSize, inclusions);
    IntArray emat1, emat2;
    double eh = rveSize/nelem; // Element size
    for (int eZ = 0; eZ < nelem; ++eZ)
        for (int eY = 0; eY < nelem; ++eY)
            for (int eX = 0; eX < nelem; ++eX) {
                int e = eX + eY*nelem + eZ*nelem*nelem + 1;
                FloatArray c = rvePosition;
                c.add({eh * (0.5 + eX), eh * (0.5 + eY), eh * (0.5 + eZ)});
                
                bool found = false;
                for ( auto &inc : rveInclusions )
                    if ( distance(c, inc.first) <= inc.second ) {
                        found = true;
                        break;
                    }
                if ( found )
                    emat2.followedBy(e);
                else
                    emat1.followedBy(e);
            }

    printf("Final inclusion fraction: %.3f\n", (double)emat2.giveSize() / (nelem * nelem * nelem));
    myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 11);
    myInput->setField(emat1, _IFT_Set_elements);
    myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));

    myInput = std::make_unique<DynamicInputRecord>(_IFT_Set_Name, 12);
    myInput->setField(emat2, _IFT_Set_elements);
    myData.insertInputRecord(DataReader::IR_setRec, std::move(myInput));
    
    //CrossSection
    for ( int i = 1; i <= 2; ++i ) {
        myInput = std::make_unique<DynamicInputRecord>(_IFT_SimpleTransportCrossSection_Name, i);
        myInput->setField(i, _IFT_SimpleTransportCrossSection_material);
        myInput->setField(10 + i, _IFT_CrossSection_SetNumber);
        myData.insertInputRecord(DataReader::IR_crosssectRec, std::move(myInput));
    }

    //Material
    myInput = std::make_unique<DynamicInputRecord>(_IFT_IsotropicHeatTransferMaterial_Name, 1);
    myInput->setField(1.0, _IFT_IsotropicHeatTransferMaterial_k);
    myInput->setField(1.0, _IFT_IsotropicHeatTransferMaterial_c);
    myInput->setField(1.0, _IFT_Material_density);
    myData.insertInputRecord(DataReader::IR_matRec, std::move(myInput));

    myInput = std::make_unique<DynamicInputRecord>(_IFT_IsotropicHeatTransferMaterial_Name, 2);
    myInput->setField(k, _IFT_IsotropicHeatTransferMaterial_k);
    myInput->setField(1.0, _IFT_IsotropicHeatTransferMaterial_c);
    myInput->setField(1.0, _IFT_Material_density);
    myData.insertInputRecord(DataReader::IR_matRec, std::move(myInput));

    //Boundary Conditions
    if ( bc == "d" || bc == "md" ) {
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
    } else if ( bc == "n" || bc == "mn" ) {
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
        
        myInput->setField(FloatArray{rveSize, rveSize, rveSize}, _IFT_TransportGradientPeriodic_jump);
        myInput->setField(9, _IFT_GeneralBoundaryCondition_set);
        myInput->setField(10, _IFT_TransportGradientPeriodic_masterSet);

        myData.insertInputRecord(DataReader::IR_bcRec, std::move(myInput));
    }
    myInput = std::make_unique<DynamicInputRecord>(_IFT_BoundaryCondition_Name, 2);
    myInput->setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
    myInput->setField(FloatArray{0.}, _IFT_BoundaryCondition_values);
    myInput->setField(IntArray{T_f}, _IFT_GeneralBoundaryCondition_dofs);
    /// @note If the mesh doesn't have a "center" node, then the centerCoords need to be changed to reflect this fixed point!
    /// Fixing the center point should only be done for Neumann b.c.s, though not actually needed with KSP-solvers.
    myInput->setField(0, _IFT_GeneralBoundaryCondition_set);
    myData.insertInputRecord(DataReader::IR_bcRec, std::move(myInput));

    //Functions
    myInput = std::make_unique<DynamicInputRecord>(_IFT_ConstantFunction_Name, 1);
    myInput->setField(1.0, _IFT_ConstantFunction_f);
    myData.insertInputRecord(DataReader::IR_funcRec, std::move(myInput));

    timer.stopTimer();
    printf("Mesh generation time %.3f s\n", timer.getUtime());
    // Writing to file (to verify, and for backups)
    //myData.writeToFile((name + "." + std :: to_string(sample) + ".in").c_str());

    printf("Initializing problem\n");
    timer.startTimer();
    auto em = InstanciateProblem(myData, _processor, 0);
    timer.stopTimer();
    printf("Instanciation time %.3f s\n", timer.getUtime());
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
        //FILE *file = fopen((this->dataOutputFileName + ".data").c_str(), "w");
        for ( int i = 0; i < tangent.giveNumberOfRows(); ++i ) {
            for ( int j = 0; j < tangent.giveNumberOfColumns(); ++j ) {
                //fprintf(file, "%.9e ", tangent(i, j) );
                fout << tangent(i, j) << " ";
            }
            fout << "\n";
            //fprintf(file, "\n");
        }
        fout.close();
        //fclose(file);
    }
    
    myData.finish();

#ifdef __PETSC_MODULE
    PetscFinalize();
#endif
#ifdef __USE_MPI
    MPI_Finalize();
#endif
}
