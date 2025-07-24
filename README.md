# Safe Bayesian Optimization

A ROS 2 package for safe path planning using Bayesian optimization with uncertainty estimation.

## Overview

This package implements a safe Bayesian optimization algorithm for robotic navigation, using terrain mapping with uncertainty to select optimal subgoals while avoiding unsafe areas.

## Build Command

```bash
colcon build --cmake-args -DBUILD_EXAMPLES=OFF
```

## Features

- **Bayesian Optimization**: Uses acquisition functions to select next best points for exploration
- **Safety Constraints**: Incorporates uncertainty estimates to avoid unsafe terrain
- **Reactive Navigation**: Diffeomorphism-based obstacle avoidance for real-time path planning
- **Obstacle Processing**: Polygon-based obstacle representation and processing
- **Terrain Mapping Integration**: Connects with terrain mapping services for real-time updates
- **Debug Visualization**: Integrated debug image publishing for development and testing
- **ROS 2 Native**: Built for ROS 2 with standard message interfaces

## Nodes

- `safe_bayesian_optimization_node`: Main optimizer node that computes safe subgoals
- `goal_point_publisher`: Publishes test goal points for navigation
- `reactive_navigation_node`: Reactive navigation node for obstacle avoidance using diffeomorphism-based path planning

## Dependencies

### System Requirements
- ROS 2 (Humble)
- C++14 compiler
- CMake >= 3.8

### Core Libraries
- **Eigen3**: Linear algebra operations
- **OpenCV**: Computer vision and image processing
- **Boost**: System utilities and geometry operations
- **CGAL**: Computational geometry algorithms
- **PCL 1.12**: Point cloud processing (common, io, surface components)
- **Qt5**: Core and Widgets components
- **Qhull**: Convex hull computation

### ROS 2 Packages
- **Standard ROS 2 Packages**:
  - `rclcpp` - ROS 2 C++ client library
  - `rclpy` - ROS 2 Python client library
  - `geometry_msgs` - Geometric primitive messages
  - `sensor_msgs` - Sensor data messages
  - `nav_msgs` - Navigation messages
  - `visualization_msgs` - Visualization markers
  - `tf2`, `tf2_ros`, `tf2_geometry_msgs` - Transform libraries
  - `cv_bridge` - OpenCV-ROS bridge
  - `turtlesim` - Turtle simulator for testing

### Custom Interfaces
- **`trusses_custom_interfaces`**: Custom service definitions for spatial data and terrain mapping

### External Dependencies (Fetched)
- **MyGAL**: Voronoi diagram generation library (fetched from GitHub)

## Usage

```bash
ros2 launch safe_bayesian_optimization safe_bayesian_optimization.launch.py
```

## Configuration

Parameters can be configured in:

### `config/safe_bayesian_optimization.yaml`:
- `opt.beta`: Exploration parameter for acquisition function
- `opt.f_min`: Minimum function value threshold
- `terrain_map.*`: Terrain map dimensions and resolution

### `config/reactive_planner.yaml`:
- `reactive_planner.p`: Diffeomorphism parameter
- `reactive_planner.epsilon`: Smoothness parameter
- `reactive_planner.varepsilon`: Small perturbation parameter
- `reactive_planner.mu_1`, `reactive_planner.mu_2`: Scaling parameters
- `reactive_planner.robot_radius`: Robot collision radius
- `reactive_planner.linear_gain`: Linear velocity gain
- `reactive_planner.angular_gain`: Angular velocity gain
- `reactive_planner.linear_cmd_limit`: Maximum linear velocity command
- `reactive_planner.angular_cmd_limit`: Maximum angular velocity command
- `reactive_planner.goal_tolerance`: Distance tolerance to stop near subgoal

## Topic Interfaces

### Published Topics

#### safe_bayesian_optimization_node
- `/current_subgoal` (`geometry_msgs/Point`) - Computed safe subgoals for reactive navigation
- `/subgoal_marker` (`visualization_msgs/Marker`) - Visualization marker for current subgoal
- `/polygon_array` (`safe_bayesian_optimization/PolygonArray`) - Obstacle polygons from unsafe terrain
- `/envelope_polygon` (`geometry_msgs/Polygon`) - Workspace envelope boundary
- `/concave_markers` (`visualization_msgs/MarkerArray`) - Safe region boundary visualization
- `/debug_polygons_image` (`sensor_msgs/Image`) - Debug visualization image (optional)

#### reactive_navigation_node
- `/turtle1/cmd_vel` (`geometry_msgs/Twist`) - Velocity commands for robot control
- `/freespace_markers` (`visualization_msgs/MarkerArray`) - Local free space visualization
- `/envelope_markers` (`visualization_msgs/MarkerArray`) - Workspace envelope visualization

#### goal_point_publisher
- `/goal_point` (`geometry_msgs/Point`) - Test goal points for navigation
- `/goal_marker` (`visualization_msgs/Marker`) - Goal point visualization marker

### Subscribed Topics

#### safe_bayesian_optimization_node
- `/goal_point` (`geometry_msgs/Point`) - Target goal points for navigation planning

#### reactive_navigation_node
- `/polygon_array` (`safe_bayesian_optimization/PolygonArray`) - Obstacle polygons from optimizer
- `/envelope_polygon` (`geometry_msgs/Polygon`) - Workspace envelope boundary
- `/turtle1/pose` (`turtlesim/Pose`) - Current robot pose from turtlesim
- `/current_subgoal` (`geometry_msgs/Point`) - Navigation subgoals from optimizer

### Service Interfaces

#### safe_bayesian_optimization_node (Client)
- `get_spatial_data` (`trusses_custom_interfaces/SpatialData`) - Requests spatial data for terrain analysis
- `get_terrain_map_with_uncertainty` (`trusses_custom_interfaces/GetTerrainMapWithUncertainty`) - Requests terrain map with uncertainty

### Custom Messages
- `safe_bayesian_optimization/PolygonArray` - Array of geometry_msgs/Polygon for obstacle representation

## License

Apache 2.0
