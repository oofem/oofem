#pragma once

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


