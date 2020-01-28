# compile everything except the main source file into a statis lib to be linked
# to the main executable
add_library(${LIBRARY_JETFIND} STATIC
jet_finding.cpp
)

set_target_properties(${LIBRARY_JETFIND} PROPERTIES PREFIX "")
