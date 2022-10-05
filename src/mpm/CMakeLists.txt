set (mpm_element
    prototype.C
    )

set (mpm
    ${mpm_element}
    )

add_library (mpm
    ${LIB_TYPE}
    ${mpm}
    )
