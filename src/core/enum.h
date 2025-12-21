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

#ifndef enum_h
#define enum_h


#include<optional>
#include<initializer_list>
#include<utility>
#include<cstring>
#include<iostream>
#include<vector>
#include<map>

namespace oofem {
   /* including enum-impl.h will define the enum SomeEnum and specialize EnumData<SomeEnum>; EnumTraits<SomeEnum> the accesses EnumData<SomeEnum> to work with the enum itself */
   template<typename T> struct EnumData;

   template<typename Enum> struct EnumTraits{
      typedef EnumData<Enum> D;
      // just forward
      static constexpr auto enum_name = D::enum_name;
      // helper function
      static const bool _starts_with_prefix(const char* s){
         return D::opt_prefix_len>0 && std::strlen(s)>D::opt_prefix_len && std::strncmp(s,D::opt_prefix,D::opt_prefix_len)==0;
      }
      // convert enum value to name
      static std::optional<const char*> name(typename D::Enum v){
         for(const auto& vn: D::value_to_name){ if(vn.value==v) return vn.name; }
         return {};
      };
      // convert integer to value
      static std::optional<typename D::Enum> value(int i){
         for(const auto& vn: D::value_to_name){ if((int)vn.value==i) return vn.value; }
         return {};
      };
      // convert name to value; considers optional prefix (opt_prefix), so e.g. VM_Total with opt_prefix "VM_" will also match "Total"
      static std::optional<typename D::Enum> value(const char* n){
         for(const auto& vn: D::value_to_name){
            if(std::strcmp(vn.name,n)==0) return vn.value;
            if(_starts_with_prefix(vn.name) && std::strcmp(vn.name+D::opt_prefix_len,n)==0) return vn.value;
         }
         return {};
      };
      // return all permissible names (with and without opt_prefix, if used)
      static std::map<int,std::vector<std::string>> all_values_to_names(){
         std::map<int,std::vector<std::string>> ret;
         for(const auto& vn: D::value_to_name){
            auto [I,_]=ret.insert({(int)vn.value,{vn.name,}});
            // add any aliases for the value here
            if(_starts_with_prefix(vn.name)) I->second.push_back(vn.name+D::opt_prefix_len);
         }
         return ret;
      }
   };

}
#endif // enum_h
