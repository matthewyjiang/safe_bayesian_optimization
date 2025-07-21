cmake_minimum_required(VERSION 3.8)
project(safe_bayesian_optimization)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(eigen3_cmake_module REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(trusses_custom_interfaces REQUIRED)
find_package(geometry_msgs REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add executables
add_executable(safe_bayesian_optimization_node src/safe_bayesian_optimization_node.cpp)
add_executable(goal_point_publisher src/goal_point_publisher.cpp)

# Add dependencies

ament_target_dependencies(
    safe_bayesian_optimization_node rclcpp trusses_custom_interfaces Eigen3 OpenCV geometry_msgs
)

ament_target_dependencies(
    goal_point_publisher rclcpp geometry_msgs
)

# Install executables
install(TARGETS
  safe_bayesian_optimization_node
  goal_point_publisher
  DESTINATION lib/${PROJECT_NAME}
)

# Install launch files
install(DIRECTORY
  launch
  DESTINATION share/${PROJECT_NAME}/
)

# Install config files
install(DIRECTORY
  config
  DESTINATION share/${PROJECT_NAME}/
)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
