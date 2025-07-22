# Safe Bayesian Optimization

A ROS 2 package for safe path planning using Bayesian optimization with uncertainty estimation.

## Overview

This package implements a safe Bayesian optimization algorithm for robotic navigation, using terrain mapping with uncertainty to select optimal subgoals while avoiding unsafe areas.

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

- ROS 2 (Humble)
- Eigen3
- OpenCV
- Boost Geometry
- Custom interfaces: `trusses_custom_interfaces`
- Custom mapping nodes and visualization: `SpiritHighLevel`

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

## License

Apache 2.0
