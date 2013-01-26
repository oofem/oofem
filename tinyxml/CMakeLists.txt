set (tinyxml
    tinystr.cpp
    tinyxml.cpp
    tinyxmlerror.cpp
    tinyxmlparser.cpp
    )

add_library (tinyxml
    ${LIB_TYPE}
    ${tinyxml}
    )
