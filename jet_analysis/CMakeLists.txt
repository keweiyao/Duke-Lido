# compile everything except the main source file into a statis lib to be linked
# to the main executable
add_library(${LIBRARY_JETFIND} STATIC
jet_finding.cpp
../src/lorentz.cpp
../src/simpleLogger.cpp
../src/hcubature.cpp
../src/TableBase.cpp
)

#set_target_properties(${LIBRARY_JETFIND} PROPERTIES PREFIX "")
# install the library headers
#file(GLOB HEADERS *.h)
#install(FILES ${HEADERS} DESTINATION include/${LIBRARY_NAME})
