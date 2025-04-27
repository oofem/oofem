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
#include <variant>
#include <vector>

#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "dictionary.h"
#include "scalarfunction.h"
#include "range.h"


 

namespace oofem {

 // ParameterManager class using std::variant with custom types
 /**
  * ParameterManager class manages the parameters for set of components, each with unique index (number).
  * It allows setting and getting parameters with different types (int, double, bool, string, IntArray, FloatArray, FloatMatrix, vector<string>, Dictionary, list<Range>, ScalarFunction).
  * The parameters are stored in a vector of unordered maps, where each map corresponds to a component index and contains the parameter name as the key and a tuple of value and priority as the value.
  * The priority is used to determine which value to keep when multiple values are set for the same parameter.
  * The class provides methods to set and get parameters, as well as to clear and resize the parameter storage.
  * The class is designed to be used in a context where parameters can be set and retrieved dynamically, such as in a simulation or modeling environment.
  * The class is not thread-safe and should be used in a single-threaded context or with appropriate synchronization mechanisms.
  * It is designed to be extensible, allowing for the addition of new parameter types in the future.
  */
 class ParameterManager {
 public:
     using ParamValue = std::variant<int, double, bool, std::string, IntArray, FloatArray, FloatMatrix, std :: vector< std :: string >, Dictionary, std :: list< Range >, ScalarFunction>;
 
     ParameterManager(int size = 0) : params(size) {}
     void clear() {params.clear(); }    
     void resize(size_t newSize) { params.resize(newSize); }

     void setParam(size_t componentIndex, const std::string& paramName, const ParamValue& value, int priority) {
         if (componentIndex >= params.size()) {
             params.resize(componentIndex + 1);
         }
         if (params[componentIndex].find(paramName) == params[componentIndex].end() || std::get<1>(params[componentIndex][paramName]) < priority) {
             params[componentIndex][paramName] = std::make_tuple(value, priority);
         }
     }
 
     ParamValue getParam(size_t componentIndex, const std::string& paramName) {
         if (componentIndex < params.size() && params[componentIndex].find(paramName) != params[componentIndex].end()) {
             return std::get<0>(params[componentIndex][paramName]);
         }
         return {};
     }
 
 private:
     std::vector<std::unordered_map<std::string, std::tuple<ParamValue, int>>> params;
 };
 




} // end namespace oofem
 

#endif // parametermanager_h