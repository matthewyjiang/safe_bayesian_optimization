from launch import LaunchDescription
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare
from launch.actions import IncludeLaunchDescription, DeclareLaunchArgument
from launch.launch_description_sources import PythonLaunchDescriptionSource
from ament_index_python.packages import get_package_share_directory
import os
from launch.launch_description_sources import FrontendLaunchDescriptionSource
from launch import LaunchDescription
from launch.substitutions import LaunchConfiguration, PathJoinSubstitution



import pathlib
import yaml


def generate_launch_description():
    config_file_arg = DeclareLaunchArgument(
            'config_file',
            default_value=PathJoinSubstitution([
                FindPackageShare('safe_legged_scouting'),
                'config',
                'safe_legged_scouting.yaml'
                ]),
            description='Path to the configuration file'
            )
    
    reactive_planner_config_arg = DeclareLaunchArgument(
            'reactive_planner_config',
            default_value=PathJoinSubstitution([
                FindPackageShare('safe_legged_scouting'),
                'config',
                'reactive_planner.yaml'
                ]),
            description='Path to the reactive planner configuration file'
            )

    yaml_dir = get_package_share_directory("safe_legged_scouting")
    config_file = os.path.join(yaml_dir, 'config/lpsc.yaml')
    print(config_file)
    # LaunchConfiguration('ros_control_config').perform(context)

    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    visualizer_params = config.get('visualizer', {}).get('ros__parameters', {})
    print(visualizer_params)
    mapping_params = config.get('mapping_node', {}).get('ros__parameters', {})
    data_collector_params = config.get('data_collector', {}).get('ros__parameters', {})

    foxglove_bridge_launch = IncludeLaunchDescription(
            FrontendLaunchDescriptionSource([
                os.path.join(
                    get_package_share_directory('foxglove_bridge'),
                    'launch',
                    'foxglove_bridge_launch.xml'
                    )
                ])
            )

    return LaunchDescription([
        config_file_arg,
        reactive_planner_config_arg, 
         Node(
            package='safe_legged_scouting',
             executable='safe_legged_scouting_node',
             name='safe_legged_scouting_node',
             parameters=[LaunchConfiguration('config_file'), LaunchConfiguration('reactive_planner_config')],
             output='screen'
             ),
         Node(
             package='safe_legged_scouting',
             executable='goal_point_publisher',
             name='goal_point_publisher',
             output='screen'
             ),
         Node(
            package='safe_legged_scouting',
            executable='reactive_navigation_node',
            name='reactive_navigation_node',
            parameters=[LaunchConfiguration('reactive_planner_config')],
            output='screen'
            ),
                Node(
            package='foxglove_visualization',  # Replace with the package where Foxglove is defined
            executable='visualizer',  # Replace with the executable name of Foxglove
            name='visualizer',
            output='screen',
            parameters=[visualizer_params]
        ),
        Node(
            package='foxglove_visualization',  # Replace with the package where FakeDataPublisher is defined
            executable='leg_measurements_publisher',  # Replace with the executable name of FakeDataPublisher
            name='leg_measurements_publisher',
            output='screen'
        ),
        Node(
            package='foxglove_visualization',  # Replace with the package where FakeDataPublisher is defined
            executable='fake_data_publisher',  # Replace with the executable name of FakeDataPublisher
            name='fake_data_publisher',
            output='screen'
        ),
        Node(
            package='mapping_collector',  # Replace with the package where FakeDataPublisher is defined
            executable='data_collector',  # Replace with the executable name of FakeDataPublisher
            name='data_collector',
            output='screen',
            parameters=[data_collector_params]
        ),

        Node(
            package='mapping_package',  # Replace with the package where FakeDataPublisher is defined
            executable='terrain_mapping_node',  # Replace with the executable name of FakeDataPublisher
            name='mapping_node',
            output='screen',
            parameters=[mapping_params]
        ),
        Node(
            package='turtlesim',
            executable='turtlesim_node',
            name='turtlesim_node',
            output='screen',
            parameters=[{'background_r': 255, 'background_g': 255, 'background_b': 255}]
            ),
        Node(
            package='safe_legged_scouting',
            executable='turtlesim_spatial_publisher.py',
            name='spirit_spatial_publisher',
            output='screen'
            ),

        foxglove_bridge_launch,


    ]) 
