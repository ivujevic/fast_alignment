if(CODE_COVERAGE)
    set(CMAKE_BUILD_TYPE Debug)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
    endif()
endif()

macro(code_coverage _TARGET_NAME)
    set(_coverage_list "")
    foreach(_orig_file ${ARGN})
        get_filename_component(_full_path ${_orig_file} ABSOLUTE)
        set(_coverage_list "${_coverage_list};${_full_path}")
    endforeach()
    add_custom_target(${_TARGET_NAME}_coverage
            COMMAND ${_TARGET_NAME}
            COMMAND xcrun llvm-profdata merge
            -o ${_TARGET_NAME}.profdata default.profraw
            COMMAND xcrun llvm-cov
            show $<TARGET_FILE:${_TARGET_NAME}>
            -instr-profile=${_TARGET_NAME}.profdata
            ${_coverage_list})
endmacro()