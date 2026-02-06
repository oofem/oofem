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

#ifndef json_h_
#define json_h_

#ifndef _USE_JSON
    #error This file may be only included with _USE_JSON
#endif

#include<nlohmann/json.hpp>

namespace oofem{

   using json = nlohmann::json;
   using ordered_json = nlohmann::ordered_json;

   struct JsonContext{
      FILE* f=nullptr;
      ordered_json ctx;
      JsonContext(FILE* f_, const ordered_json& ctx_): f(f_), ctx(ctx_){ };
      JsonContext prepend(const ordered_json& head) const {
          ordered_json ctx2(head); ctx2.update(ctx);
          return JsonContext(f,ctx2);
      }
      JsonContext append(const ordered_json& tail) const {
         ordered_json ctx2(ctx); ctx2.update(tail);
         return JsonContext(f,ctx2);
      }
      void print(const ordered_json& j) const {
         ordered_json j2(ctx);
         j2.update(j);
         fprintf(f,"%s\n",j2.dump(-1).c_str());
      }
};

}

#endif /* json_h_ */
