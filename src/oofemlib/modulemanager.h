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

#ifndef modulemanager_h
#define modulemanager_h

#include "datareader.h"
#include "error.h"

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
    std :: vector< std :: unique_ptr< M > > moduleList;
    /// Number of modules.
    int numberOfModules;
    /// Associated Engineering model.
    EngngModel *emodel;

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
    virtual M *CreateModule(const char *name, int n, EngngModel *emodel) = 0;
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
    virtual int instanciateYourself(DataReader *dr, InputRecord *ir)
    {
        IRResultType result;                   // Required by IR_GIVE_FIELD macro

        std :: string name;

        // read modules
        moduleList.reserve(numberOfModules);
        for ( int i = 0; i < numberOfModules; i++ ) {
            InputRecord *mir = dr->giveInputRecord(DataReader :: IR_expModuleRec, i + 1);
            result = mir->giveRecordKeywordField(name);
            if ( result != IRRT_OK ) {
                IR_IOERR("", mir, result);
            }

            // read type of module
            std :: unique_ptr< M > module( this->CreateModule(name.c_str(), i, emodel) );
            if ( !module ) {
                OOFEM_ERROR("unknown module (%s)", name.c_str());
            }

            module->initializeFrom(mir);
            moduleList.push_back(std :: move(module));
        }

#  ifdef VERBOSE
        VERBOSE_PRINT0("Instanciated output modules ", numberOfModules)
#  endif
        return 1;
    }

    /**
     * Instanciates the receiver from input record. Called from instanciateYourself to initialize yourself
     * from corresponding record. Should be called before instanciateYourself.
     */
    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
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
} // end namespace oofem
#endif // modulemanager_h
