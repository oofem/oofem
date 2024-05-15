#include "oofem_version.h"
#include "oofemcfg.h"

const char* PRG_VERSION =  "OOFEM version" __OOFEM_MAJOR_VERSION "." __OOFEM_MINOR_VERSION;
const char* OOFEM_VERSION = __OOFEM_VERSION;
const char* OOFEM_COPYRIGHT = __OOFEM_COPYRIGHT;
const char* OOFEM_GIT_HASH = __OOFEM_GIT_HASH;
const char* OOFEM_GIT_REPOURL = __OOFEM_GIT_REPOURL;
const char* OOFEM_GIT_BRANCH = __OOFEM_GIT_BRANCH;
#


const char* PRG_HEADER_SM = ""
"____________________________________________________\n"
"  ____  ____  __________  ___ \n"
" / __ \\/ __ \\/ __/ __/  |/  /___  _______\n"
"/ /_/ / /_/ / _// _// /|_/ // _ \\/ __/ _ `\n"
"\\____/\\____/_/ /___/_/  /_(_)___/_/  \\_, / \n" 
__OOFEM_COPYRIGHT " /__/\n"
"____________________________________________________\n";

const char* PRG_HEADER = ""
"############################################################## www.oofem.org ###\n"
"  ____  ____  __________  ___\n"
" / __ \\/ __ \\/ __/ __/  |/  /___  _______\n"
"/ /_/ / /_/ / _// _// /|_/ // _ \\/ __/ _ `\n"
"\\____/\\____/_/ /___/_/  /_(_)___/_/  \\_, /  OOFEM ver. " __OOFEM_VERSION "\n"
"                                    /___/   " __OOFEM_COPYRIGHT "\n"
"################################################################################\n"
"Git URL:" __OOFEM_GIT_REPOURL ",\n    branch:" __OOFEM_GIT_BRANCH ", hash:" __OOFEM_GIT_HASH "\n\n";

                               
const char* HOST_TYPE = __HOST_TYPE;
const char* HOST_NAME = __HOST_NAME;
const char* MODULE_LIST = __MODULE_LIST;

