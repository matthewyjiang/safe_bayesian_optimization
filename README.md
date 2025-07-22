# Safe Bayesian Optimization

A ROS 2 package for safe path planning using Bayesian optimization with uncertainty estimation.

## Overview

This package implements a safe Bayesian optimization algorithm for robotic navigation, using terrain mapping with uncertainty to select optimal subgoals while avoiding unsafe areas.

## Features

- **Bayesian Optimization**: Uses acquisition functions to select next best points for exploration
- **Safety Constraints**: Incorporates uncertainty estimates to avoid unsafe terrain
- **Terrain Mapping Integration**: Connects with terrain mapping services for real-time updates
- **ROS 2 Native**: Built for ROS 2 with standard message interfaces

## Nodes

- `safe_bayesian_optimization_node`: Main optimizer node that computes safe subgoals
- `goal_point_publisher`: Publishes goal points for navigation

## Dependencies

- ROS 2 (Humble/Iron)
- Eigen3
- OpenCV
- Custom interfaces: `trusses_custom_interfaces`

## Usage

```bash
ros2 launch safe_bayesian_optimization safe_bayesian_optimization.launch.py
```

## Configuration

Parameters can be configured in `config/safe_bayesian_optimization.yaml`:
- `opt.beta`: Exploration parameter for acquisition function
- `opt.f_min`: Minimum function value threshold
- `terrain_map.*`: Terrain map dimensions and resolution

## License

Apache 2.0
