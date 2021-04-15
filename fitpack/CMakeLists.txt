enable_language(Fortran)

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

add_library(fitpack STATIC ${SRC})

include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE ""
                         SYMBOLS curfit splev)
FortranCInterface_VERIFY(CXX)

target_include_directories(fitpack
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

