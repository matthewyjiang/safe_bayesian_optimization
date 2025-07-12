from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration, PathJoinSubstitution
from launch_ros.substitutions import FindPackageShare


def generate_launch_description():
    config_file_arg = DeclareLaunchArgument(
        'config_file',
        default_value=PathJoinSubstitution([
            FindPackageShare('safe_bayesian_optimization'),
            'config',
            'safe_bayesian_optimization.yaml'
        ]),
        description='Path to the configuration file'
    )

    safe_bayesian_optimization_node = Node(
        package='safe_bayesian_optimization',
        executable='safe_bayesian_optimization_node',
        name='safe_bayesian_optimization_node',
        parameters=[LaunchConfiguration('config_file')],
        output='screen'
    )

    return LaunchDescription([
        config_file_arg,
        safe_bayesian_optimization_node
    ])