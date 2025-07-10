#include"oofemcfg.h"
// macro values are defiend in CMakeLists.txt

OOFEM_EXPORT const char* PRG_VERSION =  "OOFEM version " __OOFEM_MAJOR_VERSION "." __OOFEM_MINOR_VERSION;
OOFEM_EXPORT const char* OOFEM_VERSION = __OOFEM_VERSION;
OOFEM_EXPORT const char* OOFEM_COPYRIGHT = __OOFEM_COPYRIGHT;
OOFEM_EXPORT const char* OOFEM_GIT_HASH = __OOFEM_GIT_HASH;
OOFEM_EXPORT const char* OOFEM_GIT_REPOURL = __OOFEM_GIT_REPOURL;
OOFEM_EXPORT const char* OOFEM_GIT_BRANCH = __OOFEM_GIT_BRANCH;
#


OOFEM_EXPORT const char* PRG_HEADER_SM = ""
"____________________________________________________\n"
"  ____  ____  __________  ___ \n"
" / __ \\/ __ \\/ __/ __/  |/  /___  _______\n"
"/ /_/ / /_/ / _// _// /|_/ // _ \\/ __/ _ \\\n"
"\\____/\\____/_/ /___/_/  /_(_)___/_/  \\_, / \n" 
__OOFEM_COPYRIGHT " /__/\n"
"____________________________________________________\n";

OOFEM_EXPORT const char* PRG_HEADER = ""
"############################################################## www.oofem.org ###\n"
"  ____  ____  __________  ___\n"
" / __ \\/ __ \\/ __/ __/  |/  /___  _______\n"
"/ /_/ / /_/ / _// _// /|_/ // _ \\/ __/ _ \\\n"
"\\____/\\____/_/ /___/_/  /_(_)___/_/  \\_, /  OOFEM ver. " __OOFEM_VERSION "\n"
"                                    /___/   " __OOFEM_COPYRIGHT "\n"
"################################################################################\n"
"Git URL:" __OOFEM_GIT_REPOURL ",\n    branch:" __OOFEM_GIT_BRANCH ", hash:" __OOFEM_GIT_HASH "\n\n";

                               
OOFEM_EXPORT const char* HOST_TYPE = __HOST_TYPE;
OOFEM_EXPORT const char* HOST_NAME = __HOST_NAME;
OOFEM_EXPORT const char* MODULE_LIST = __MODULE_LIST;

