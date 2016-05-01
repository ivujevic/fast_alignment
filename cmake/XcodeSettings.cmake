macro(adjust_xcode_settings _TARGET_NAME)
    if(CMAKE_GENERATOR STREQUAL Xcode)
        set_property(TARGET ${_TARGET_NAME} PROPERTY
                XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT "dwarf-with-dsym")
        set_property(TARGET ${_TARGET_NAME} PROPERTY
                XCODE_ATTRIBUTE_GCC_GENERATE_DEBUGGING_SYMBOLS[variant=MinSizeRel] YES)
        set_property(TARGET ${_TARGET_NAME} PROPERTY
                XCODE_ATTRIBUTE_GCC_GENERATE_DEBUGGING_SYMBOLS[variant=Release] YES)
    endif()
endmacro()