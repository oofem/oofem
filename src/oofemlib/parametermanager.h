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

#ifndef parametermanager_h
#define parametermanager_h

#include <unordered_map>
#include <string>
#include <tuple>
#include <vector>
#include <mutex>
#include <shared_mutex>
#include <variant>
#include <optional>

#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include "error.h"

namespace oofem {


#define PM_UPDATE_PARAMETER(_val, _pm, _ir, _componentnum, _paramkey, _prio) \
    { \
        std::size_t _indx=_paramkey.getIndex(); \
        const char* _kwd = _paramkey.getName().c_str(); \
        if ((_prio > _pm.getPriority(_componentnum-1, _indx)) && (_ir.hasField(_kwd))) { \
            _ir.giveField(_val, _kwd); \
            _pm.setPriority(_componentnum-1, _indx, _prio); \
        } \
    }

#define PM_UPDATE_PARAMETER_AND_REPORT(_val, _pm, _ir, _componentnum, _paramkey, _prio, _flag) \
    { \
        std::size_t _indx=_paramkey.getIndex(); \
        const char* _kwd = _paramkey.getName().c_str(); \
        if ((_prio > _pm.getPriority(_componentnum-1, _indx)) && (_ir.hasField(_kwd))) { \
            _ir.giveField(_val, _kwd); \
            _pm.setPriority(_componentnum-1, _indx, _prio); \
            _flag=true; \
        } else { \
            _flag=false; \
        } \
    }

#define PM_UPDATE_TEMP_PARAMETER(_type, _pm, _ir, _componentnum, _paramkey, _prio) \
    { \
        std::size_t _indx=_paramkey.getIndex(); \
        const char* _kwd = _paramkey.getName().c_str(); \
        if ((_prio > _pm.getPriority(_componentnum-1, _indx)) && (_ir.hasField(_kwd))) { \
            _type _val; \
            _ir.giveField(_val, _kwd); \
            _pm.setPriority(_componentnum-1, _indx, _prio); \
            _pm.setTemParam(_componentnum-1, _indx, _val); \
        } \
    }

#define PM_CHECK_FLAG_AND_REPORT(_pm, _ir, _componentnum, _paramkey, _prio, _flag) \
    { \
        std::size_t _indx=_paramkey.getIndex(); \
        const char* _kwd = _paramkey.getName().c_str(); \
        if ((_prio > _pm.getPriority(_componentnum-1, _indx)) && (_ir.hasField(_kwd))) { \
            _pm.setPriority(_componentnum-1, _indx, _prio); \
            _flag=true; \
        } else { \
            _flag=false; \
        } \
    }

#define PM_ELEMENT_ERROR_IFNOTSET(_pm, _componentnum, _paramkey) \
    { \
        if (!_pm.checkIfSet(_componentnum-1, _paramkey.getIndex())) { \
            OOFEM_ERROR("Element %d: Parameter %s not set", _componentnum, _paramkey.getName());\
        }\
    }    

#define PM_DOFMAN_ERROR_IFNOTSET(_pm, _componentnum, _paramkey) \
    { \
        if (!_pm.checkIfSet(_componentnum-1, _paramkey.getIndex())) { \
            OOFEM_ERROR("DofManager %d: Parameter %s not set", _componentnum, _paramkey.getName());\
        }\
    }    

 /**
  * ParameterPriorityManager class manages the actual priority for parameters that can be set from multiple sources, typically from component record or via set.
  * It allows setting and getting parameters with a priority system, where higher priority values override lower priority ones.
  * The manager manages parameters for multiple components, each identified by a unique index. 
  * The parameters are stored in a shared map, where the key is a string representing the parameter name and the value is a the parameter index.
  * It uses a shared mutex for thread-safe access to the parameters.
  * The class is designed to be used in a multi-threaded environment where multiple threads may attempt to set or get parameters simultaneously.
  * The class provides methods to register parameters, set their priorities, and retrieve their priorities for a collection of components.
  * 
  * It also allows to manage temporary parameters for each component which can be set and retrieved.
  * The temporary params are handy for pattern with multi-source initialization and single finalization.
  */
// ParameterManager class to track parameter priorities
class ParameterManager {
public:
    using paramValue = std::variant<int, double, std::string, bool, IntArray, FloatArray, FloatMatrix>;
    void setPriority(size_t componentIndex, size_t paramIndex, int priority) {
        std::unique_lock lock(mtx);
        if (componentIndex >= priorities.size()) {
            priorities.resize(componentIndex + 1);
        }
        priorities[componentIndex][paramIndex] = priority;
    }

    int getPriority(size_t componentIndex, size_t paramIndex) const {
        std::shared_lock lock(mtx);
        if (componentIndex < priorities.size() && priorities[componentIndex].find(paramIndex) != priorities[componentIndex].end()) {
            return priorities[componentIndex].at(paramIndex);
        }
        return -1; // Return -1 if priority is not found
    }

    void clear() {
        std::unique_lock lock(mtx);
        priorities.clear();
        tempParams.clear();
    }

    bool checkIfSet(size_t componentIndex, size_t paramIndex) {       
        return priorities[componentIndex].find(paramIndex) != priorities[componentIndex].end();
    }

    void setTemParam(size_t componentIndex, size_t paramIndex, const paramValue &value) {
        std::unique_lock lock(mtx);
        tempParams[paramIndex] = value;
    }
    std::optional<paramValue> getTempParam(size_t componentIndex, size_t paramIndex) const {
        std::shared_lock lock(mtx);
        auto it = tempParams.find(paramIndex);
        if (it != tempParams.end()) {
            return it->second;
        }
        return std::nullopt; // Return nullopt if parameter is not found
    }
    bool hasTempParam(size_t componentIndex, size_t paramIndex) const {
        std::shared_lock lock(mtx);
        return tempParams.find(paramIndex) != tempParams.end();
    }

private:
    std::vector<std::unordered_map<size_t, int>> priorities;
    mutable std::shared_mutex mtx;

    std::unordered_map<size_t, paramValue> tempParams;
};




} // end namespace oofem
 

#endif // parametermanager_h