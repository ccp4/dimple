cmake_minimum_required(VERSION 2.8)
project(dimple C CXX)

find_package(CCP4 REQUIRED mmdb ccp4c clipper-core clipper-ccp4)
include_directories(${CCP4_INCLUDE_DIRS})
add_definitions(-DBUILD_EXE)

add_executable(find-blobs c4/find-blobs.cc)
target_link_libraries(find-blobs ${CCP4_LIBRARIES})

if (WIN32)
    install(PROGRAMS dimple.bat DESTINATION bin)
endif()
install(TARGETS find-blobs DESTINATION bin)
install(PROGRAMS dimple DESTINATION bin)
install(FILES dimple.py ccp4i-dimple.py DESTINATION ${CCP4_PY2_DIR})
install(DIRECTORY c4 DESTINATION ${CCP4_PY2_DIR}
        FILES_MATCHING PATTERN "*.py")
