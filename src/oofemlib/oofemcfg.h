#ifndef oofemcfg_h
#define oofemcfg_h

#ifdef _USE_SHARED
    #include "oofem_export.h"
#else
    #define OOFEM_EXPORT
    #define OOFEM_NO_EXPORT
#endif

OOFEM_EXPORT extern const char* PRG_VERSION;
OOFEM_EXPORT extern const char* OOFEG_VERSION;
OOFEM_EXPORT extern const char* OOFEM_COPYRIGHT;
OOFEM_EXPORT extern const char* PRG_HEADER;
OOFEM_EXPORT extern const char* PRG_HEADER_SM;
OOFEM_EXPORT extern const char* HOST_TYPE;
OOFEM_EXPORT extern const char* HOST_NAME;
OOFEM_EXPORT extern const char* MODULE_LIST;
OOFEM_EXPORT extern const char* OOFEM_GIT_HASH;
OOFEM_EXPORT extern const char* OOFEM_GIT_REPOURL;
OOFEM_EXPORT extern const char* OOFEM_GIT_BRANCH;
#endif /* oofemcfg.h */
