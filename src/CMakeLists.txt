# compile everything except the main source file into a statis lib to be linked
# to the main executable
add_library(${LIBRARY_NAME} STATIC
JetFunctions.cpp
MediumLoader.cpp
)

set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")

# compile the actual executable
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} ${Boost_INCLUDE_DIRS} ${GSL_LIBRARIES} ${GSLCALAS_LIBRARIES} ${HDF5_LIBRARIES} -pthread -lpthread)

# install executable
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
# also install the generated c++ shared library
install(TARGETS ${LIBRARY_NAME} DESTINATION lib/${LIBRARY_NAME})
# install the library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/${LIBRARY_NAME})
