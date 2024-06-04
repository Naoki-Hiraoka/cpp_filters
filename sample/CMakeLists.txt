set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

include_directories(
  ../include
)

add_executable(sample1
  Sample1.cpp
  )

target_link_libraries(sample1
  ${PROJECT_NAME}
  )

install(TARGETS sample1
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )
