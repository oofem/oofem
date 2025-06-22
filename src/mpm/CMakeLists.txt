set (mpm_core
    mpm.C
    integral.C
    )


set (mpm_element
    termlibrary.C
    termlibrary2.C
    termlibrary3.C
    prototype.C
    prototype2.C
    up.C
    tm.C
    sadg.C
    termlibrary4.C
    )

    set (mpm_problems
    mpmproblem.C
    dg.C
    mpmsemodels.C
    )

set (mpm
    ${mpm_core}    
    ${mpm_element}
    ${mpm_problems}
    )

add_library (mpm
    ${LIB_TYPE}
    ${mpm}
    )

