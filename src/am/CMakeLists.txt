set (am_sources
    additivemanufacturingproblem.C
    freeconstantsurfaceload.C
    depositedheatsource.C
    )

set (tm 
    ${am_sources}
    )

add_library (am
    ${LIB_TYPE}
    ${tm}
    )

