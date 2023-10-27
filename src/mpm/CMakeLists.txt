set (mpm_element
    termlibrary.C
    prototype.C
    up.C
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
