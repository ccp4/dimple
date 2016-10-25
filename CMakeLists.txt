cmake_minimum_required(VERSION 2.8)
project(dimple C CXX)

find_package(CCP4 REQUIRED mmdb2 ccp4c clipper-core clipper-ccp4)
include_directories(${CCP4_INCLUDE_DIRS})
add_definitions(-DBUILD_EXE)

add_executable(find-blobs cxx/find-blobs.cc)
target_link_libraries(find-blobs ${CCP4_LIBRARIES})

install(TARGETS find-blobs DESTINATION bin)
install(PROGRAMS dimple DESTINATION bin)
if (WIN32)
  install(PROGRAMS dimple.bat DESTINATION bin)
endif()
install(FILES ccp4i-dimple.py DESTINATION ${CCP4_PY2_DIR})
install(FILES __init__.py __main__.py coots.py pdb.py utils.py
              cell.py main.py mtz.py workflow.py
        DESTINATION ${CCP4_PY2_DIR}/dimple)
install(FILES contaminants/__init__.py contaminants/data.py
        DESTINATION ${CCP4_PY2_DIR}/dimple/contaminants)
