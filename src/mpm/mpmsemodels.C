#include "mpmsemodels.h"
#include "classfactory.h"

namespace oofem {
    REGISTER_EngngModel(StationaryMPMSProblem);

    int
    MPMSymbolicProblem:: instanciateMPM (DataReader &dr, InputRecord &ir) {
        std::shared_ptr<InputRecord> irPtr(ir.ptr());
        std::string name;
        int num=-1;
        DataReader::RecordGuard scope(dr,irPtr.get());
        for(auto& mir: dr.giveGroupRecords(irPtr,"nvariables","MPMVariables",DataReader::IR_mpmVarRec,/*optional*/true)){
            IR_GIVE_FIELD(mir, name, "name");
            std::unique_ptr< Variable > var = std :: make_unique< Variable >();
            var->initializeFrom(mir);
            variableMap[name] = std::move(var);
        }
        //if(variableMap.empty()) OOFEM_ERROR("No MPM Variables defined.");
        for(auto& mir: dr.giveGroupRecords(irPtr,"nterms","MPMTerms",DataReader::IR_mpmTermRec,/*optional*/true)){
            IR_GIVE_RECORD_KEYWORD_FIELD(mir, name, num);
            std::unique_ptr< Term > term = classFactory.createTerm(name.c_str());
            term->initializeFrom(mir, this);
            termList.push_back(std::move(term));
        }
        //if(termList.empty()) OOFEM_ERROR("No MPM Terms defined.");
        for(auto& mir: dr.giveGroupRecords(irPtr,"nintegrals","MPMIntegrals",DataReader::IR_mpmIntegralRec,/*optional*/true)){
            std::unique_ptr< Integral > integral = std :: make_unique< Integral >(nullptr, &dummySet, nullptr);
            integral->initializeFrom(mir, this);
            this->addIntegral(std::move(integral));
        }
        //if(integralList.empty()) OOFEM_ERROR("No MPM Integrals defined.");
        return 1;
    }

} // end namespace oofem
