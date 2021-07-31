add_compile_options(-std=c++11)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

include_directories(
  ../include
)

add_executable(sample1
  Sample1.cpp
)

install(TARGETS sample1
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )
