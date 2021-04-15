enable_language(CXX)
set(LIB_NAME fitpackpp)

set(SRC
    src/BSplineCurve.cpp
)

add_library(${LIB_NAME} STATIC ${SRC})
target_include_directories(${LIB_NAME}
    PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(${LIB_NAME}
    PRIVATE
    fitpack
)

install(TARGETS ${LIB_NAME}
    EXPORT fitpackpp-targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include DESTINATION .)

