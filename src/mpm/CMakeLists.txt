set (mpm_element
    termlibrary.C
    termlibrary2.C
    termlibrary3.C
    prototype.C
    prototype2.C
    up.C
    tm.C
    )

    set (mpm_problems
    mpmproblem.C
    )

set (mpm
    ${mpm_element}
    ${mpm_problems}
    )

add_library (mpm
    ${LIB_TYPE}
    ${mpm}
    )

add_dependencies(mpm version)
