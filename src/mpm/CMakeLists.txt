set (mpm_element
    termlibrary.C
    prototype.C
    )

set (mpm
    ${mpm_element}
    )

add_library (mpm
    ${LIB_TYPE}
    ${mpm}
    )
