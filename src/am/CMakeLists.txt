set (am_sources
    additivemanufacturingproblem.C
    freeconstantsurfaceload.C
    )

set (tm 
    ${am_sources}
    )

add_library (am
    ${LIB_TYPE}
    ${tm}
    )

add_dependencies(am version)