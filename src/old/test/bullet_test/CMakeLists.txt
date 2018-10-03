cmake_minimum_required(VERSION 2.8)
project(bullet_test)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-g -Wall -Werror -std=c++11")

find_package (Bullet REQUIRED)
find_package (Eigen3 REQUIRED)

include_directories(
  ${BULLET_INCLUDE_DIRS} 
  ${EIGEN3_INCLUDE_DIR} 
)

message(${BULLET_INCLUDE_DIRS})

# add_executable(bt_test bullet_collision_checker.cpp include/types.h)
add_executable(bt_test hello_bullet.cpp)
target_link_libraries(bt_test ${BULLET_LIBRARIES})
