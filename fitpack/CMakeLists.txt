enable_language(Fortran)
set(LIB_NAME fitpack)

set(SRC
  src/splev.f
  src/curfit.f
  src/fpbspl.f
  src/fpchec.f
  src/fpcurf.f
  src/fprota.f
  src/fpback.f
  src/fpknot.f
  src/fpdisc.f
  src/fpgivs.f
  src/fprati.f
)

add_library(${LIB_NAME} STATIC ${SRC})

include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE ""
                         SYMBOLS curfit splev)
FortranCInterface_VERIFY(CXX)

target_include_directories(${LIB_NAME}
    PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

install(TARGETS ${LIB_NAME}
    EXPORT fitpackpp-targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
