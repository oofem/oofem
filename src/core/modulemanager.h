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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef modulemanager_h
#define modulemanager_h

#include "datareader.h"
#include "error.h"
#include "exportmodule.h"


#include <memory>
#include <vector>

///@name Input fields for module managers
//@{
#define _IFT_ModuleManager_nmodules "nmodules"
//@}

namespace oofem {
class EngngModel;

/**
 * Class representing and implementing ModuleManager. It is attribute of EngngModel.
 * It manages the modules of given type.
 */
template< class M >
class OOFEM_EXPORT ModuleManager
{
protected:
    /// Module list.
    std::vector< std::unique_ptr< M > > moduleList;
    /// Number of modules.
    int numberOfModules = 0;
    /// Associated Engineering model.
    EngngModel *emodel = nullptr;

public:
    ModuleManager(EngngModel * emodel) {
        this->emodel = emodel;
        numberOfModules = 0;
    }
    virtual ~ModuleManager() {
    }
    ModuleManager(const ModuleManager &) = delete;
    ModuleManager & operator=(const ModuleManager &) = delete;
    /**
     * Creates new instance of module.
     * @param name Name of module.
     * @param n Number associated with module.
     * @param emodel Engineering model which receiver belongs to.
     */
    virtual std::unique_ptr<M> CreateModule(const char *name, int n, EngngModel *emodel) = 0;
    /**
     *
     * Reads receiver description from input stream and creates corresponding modules components accordingly.
     * It scans input file, each line is assumed to be single record describing particular module.
     * The record line is converted to lowercase letters.
     * After new output module object is created, its initializeForm member function is
     * called with its record as parameter.
     * @param dr Data reader for input records.
     * @param ir Record for receiver.
     * @return Nonzero if o.k.
     */
    virtual int instanciateYourself(DataReader &dr, const std::shared_ptr<InputRecord>& irPtr, InputFieldType ift, const std::string& name, DataReader::InputRecordType irType)
    {
        // read modules
        DataReader::GroupRecords modRecs=dr.giveGroupRecords(irPtr,ift,name,irType,/*optional*/true);
        moduleList.reserve(modRecs.size());
        int modIndex0=0;
        for (auto& mir: modRecs){
            std::string modName;
            mir.giveRecordKeywordField(modName);

            // read type of module
            std :: unique_ptr< M > module = this->CreateModule(modName.c_str(), modIndex0, emodel);
            if ( !module ) {
                OOFEM_ERROR("unknown module (%s)", modName.c_str());
            }

            module->initializeFrom(mir);
            registerModule(module);
            modIndex0++;
        }

        #  ifdef VERBOSE
            VERBOSE_PRINT0("Instanciated %d modules", numberOfModules);
        #  endif
        return 1;
    }

    /**
     * Stores a module in moduleList. Useful when 
     * adding modules externally, e.g. from Python
     */
    
    virtual void registerModule (std::unique_ptr< M > &module){
        moduleList.push_back(std::move(module));
    }
    
    
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const = 0;

    /**
     * Returns the required module.
     * @param num Module number.
     */
    M *giveModule(int num) {
        M *elem = NULL;

        if ( num >= 1 && num <= (int)moduleList.size() ) {
            elem = moduleList[num-1].get();
        } else {
            OOFEM_ERROR("No module no. %d defined", num);
        }

        return elem;
    }

    int giveNumberOfModules() const { return (int)moduleList.size(); }
};

template class OOFEM_EXPORT ModuleManager<ExportModule>;

} // end namespace oofem
#endif // modulemanager_h
