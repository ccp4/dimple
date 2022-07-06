cmake_minimum_required(VERSION 3.15)
project(dimple C CXX)

find_package(CCP4 REQUIRED mmdb2 ccp4c clipper-core clipper-ccp4)
include_directories(${CCP4_INCLUDE_DIRS})
add_definitions(-DBUILD_EXE)

find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter)
if(WIN32)
    string(REPLACE "\\" "/" Python_SITELIB "${Python_SITELIB}")
endif(WIN32)

add_executable(find-blobs cxx/find-blobs.cc)
target_link_libraries(find-blobs ${CCP4_LIBRARIES})

install(TARGETS find-blobs DESTINATION bin)
install(PROGRAMS dimple DESTINATION bin)
if (WIN32)
  install(PROGRAMS dimple.bat DESTINATION bin)
endif()
install(FILES ccp4i-dimple.py DESTINATION ${Python_SITELIB})
install(FILES __init__.py __main__.py coots.py pdb.py utils.py
              cell.py main.py mtz.py workflow.py
        DESTINATION ${Python_SITELIB}/dimple)
install(FILES contaminants/__init__.py contaminants/data.py
              contaminants/prepare.py
        DESTINATION ${Python_SITELIB}/dimple/contaminants)
