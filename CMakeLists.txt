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
install(FILES dimple.py DESTINATION lib/python2.7/site-packages)
install(DIRECTORY c4 DESTINATION lib/python2.7/site-packages
        FILES_MATCHING PATTERN "*.py")
