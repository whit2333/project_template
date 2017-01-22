set(aname DERP)
set(needs_libs )
set(libname "${aname}")
set(dictname "${libname}Dict")
set(lib_LINKDEF "${PROJECT_SOURCE_DIR}/src/${aname}/include/LinkDef.h")
set(lib_DICTIONARY_SRC "${libname}Dict.cxx")
set(lib_PCM_FILE "${PROJECT_BINARY_DIR}/src/${aname}/${libname}Dict_rdict.pcm")

include_directories(${PROJECT_SOURCE_DIR}/src/${aname}/include)

# add src file names here (without the extention or src/ prefix)
set(lib_files
  Derp
  #qcdnum2c
)

set(lib_Fortran_files
   )

set(lib_SRCS)
set(lib_HEADERS)
foreach(infileName ${lib_files})
   SET(lib_SRCS ${lib_SRCS} ${PROJECT_SOURCE_DIR}/src/${aname}/src/${infileName}.cxx)
   SET(lib_HEADERS ${lib_HEADERS} ${PROJECT_SOURCE_DIR}/src/${aname}/include/${infileName}.h)
endforeach(infileName)

set(lib_Fortran_SRCs)
foreach(infileName ${lib_Fortran_files})
   SET(lib_Fortran_SRCs ${lib_Fortran_SRCs} "${PROJECT_SOURCE_DIR}/src/${aname}/src/${infileName}.f")
endforeach(infileName)


#MESSAGE("(${dictname} ${lib_HEADERS} LINKDEF ${lib_LINKDEF}) ")
ROOT_GENERATE_DICTIONARY(${dictname} ${lib_HEADERS} LINKDEF ${lib_LINKDEF} OPTIONS -p)

SET(lib_HEADERS ${lib_HEADERS})
SET(lib_SRCS ${lib_Fortran_SRCs} ${lib_SRCS} ${lib_DICTIONARY_SRC})

SET(lib_VERSION "${${PROJECT_NAME}_VERSION}")
SET(lib_MAJOR_VERSION "${${PROJECT_NAME}_MAJOR_VERSION}")
SET(lib_LIBRARY_PROPERTIES 
    VERSION "${lib_VERSION}"
    SOVERSION "${lib_MAJOR_VERSION}"
    SUFFIX ".so")

ADD_CUSTOM_TARGET(${aname}_ROOTDICTS DEPENDS ${lib_SRCS} ${lib_HEADERS} ${lib_DICTIONARY_SRC} ${lib_DICTIONARY_HEADER})

ADD_LIBRARY(${libname} SHARED ${lib_SRCS})
TARGET_LINK_LIBRARIES(${libname} ${LINK_LIBRARIES} ${needs_libs})
SET_TARGET_PROPERTIES(${libname} PROPERTIES ${lib_LIBRARY_PROPERTIES})
ADD_DEPENDENCIES(${libname} ${needs_libs} ${aname}_ROOTDICTS)

INSTALL(
   TARGETS ${libname} 
   EXPORT ${PROJECT_NAME}Targets
   DESTINATION lib )

INSTALL(
   FILES ${lib_PCM_FILE} 
   DESTINATION lib )

file(GLOB template_headers ${PROJECT_SOURCE_DIR}/src/${aname}/include/*.hxx)
SET(lib_HEADERS ${lib_HEADERS} ${template_headers})
install(
   FILES ${lib_HEADERS} 
   DESTINATION include/${PROJECT_NAME} )

#install(TARGETS foo
#  # IMPORTANT: Add the foo library to the "export-set"
#  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
#  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
#  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/foo"
#    COMPONENT dev)

