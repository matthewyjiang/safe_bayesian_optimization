#!/usr/bin/env python3

import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Twist
from turtlesim.msg import Pose
from trusses_custom_interfaces.msg import SpatialMeasurement
import math
import time
import numpy as np
import random
import os


class TurtlesimSpatialPublisher(Node):
    
    def __init__(self):
        super().__init__('turtlesim_spatial_publisher')
        
        # Load terrain data
        self.load_terrain_data()
        
        # Publisher for spatial measurements
        self.spatial_pub = self.create_publisher(
            SpatialMeasurement, 
            '/spatial_measurements', 
            10
        )
        
        # Subscriber to turtlesim pose
        self.pose_sub = self.create_subscription(
            Pose,
            '/turtle1/pose',
            self.pose_callback,
            10
        )
        
        # Store current pose
        self.current_pose = None
        self.initial_points_published = False
        
        # Timer for publishing spatial measurements
        self.timer = self.create_timer(1, self.publish_spatial_measurement)
        
        self.get_logger().info('Turtlesim Spatial Publisher initialized')
    
    def load_terrain_data(self):
        """Load terrain data from CSV file"""
        try:
            # Use ament_index to find package share directory
            from ament_index_python.packages import get_package_share_directory
            package_share_dir = get_package_share_directory('safe_bayesian_optimization')
            terrain_file = os.path.join(package_share_dir, 'data', 'terrain.csv')
            
            # Load terrain data using numpy
            self.terrain_data = np.loadtxt(terrain_file, delimiter=',', skiprows=1)
            self.get_logger().info(f'Loaded terrain data with shape: {self.terrain_data.shape}')
            
            # Define terrain grid bounds (assuming turtlesim world is 11x11, terrain data covers this)
            self.terrain_rows, self.terrain_cols = self.terrain_data.shape
            self.x_min, self.x_max = 0.0, 11.0
            self.y_min, self.y_max = 0.0, 11.0
            
        except Exception as e:
            self.get_logger().error(f'Failed to load terrain data: {e}')
            self.terrain_data = None
    
    def get_terrain_value(self, x, y):
        """Get terrain value at given x, y coordinates"""
        if self.terrain_data is None:
            return 100.0  # Default value
        
        # Convert world coordinates to array indices
        col = int((x - self.x_min) / (self.x_max - self.x_min) * (self.terrain_cols - 1))
        row = int((y - self.y_min) / (self.y_max - self.y_min) * (self.terrain_rows - 1))
        
        # Clamp to valid indices
        col = max(0, min(col, self.terrain_cols - 1))
        row = max(0, min(row, self.terrain_rows - 1))
        
        return float(self.terrain_data[row, col])
    
    def pose_callback(self, msg):
        """Store the current turtle pose"""
        self.current_pose = msg
        
        # Publish initial random points around robot on first pose received
        if not self.initial_points_published:
            self.publish_initial_random_points()
            self.initial_points_published = True
    
    def publish_initial_random_points(self):
        """Publish random spatial measurements around the robot's initial position"""
        if self.current_pose is None:
            return
        
        num_points = 10  # Number of random points to publish
        max_distance = 2.0  # Maximum distance from robot
        
        for i in range(num_points):
            # Generate random point in vicinity of robot
            angle = random.uniform(0, 2 * math.pi)
            distance = random.uniform(0.5, max_distance)
            
            # Calculate random point coordinates
            x = self.current_pose.x + distance * math.cos(angle)
            y = self.current_pose.y + distance * math.sin(angle)
            
            # Ensure point is within turtlesim bounds (0-11)
            x = max(0.1, min(10.9, x))
            y = max(0.1, min(10.9, y))
            
            # Create and publish spatial measurement
            spatial_msg = SpatialMeasurement()
            spatial_msg.position.x = x
            spatial_msg.position.y = y
            spatial_msg.position.z = 0.0
            spatial_msg.value = self.get_terrain_value(x, y)
            spatial_msg.unit = 'm'
            spatial_msg.source_name = 'turtlesim_init'
            spatial_msg.uncertainty = 3.0  # Lower uncertainty for initial points
            spatial_msg.time = self.get_clock().now().to_msg()
            
            self.spatial_pub.publish(spatial_msg)
            
        self.get_logger().info(f'Published {num_points} initial random points around robot')
    
    def publish_spatial_measurement(self):
        """Publish spatial measurement based on turtle position"""
        if self.current_pose is None:
            return
        
        # Create spatial measurement message
        spatial_msg = SpatialMeasurement()
        
        # Set position from turtle pose
        spatial_msg.position.x = self.current_pose.x
        spatial_msg.position.y = self.current_pose.y
        spatial_msg.position.z = 0.0  # Turtlesim is 2D
        
        # Get terrain value from loaded data
        spatial_msg.value = self.get_terrain_value(self.current_pose.x, self.current_pose.y)
        
        spatial_msg.unit = 'm'
        spatial_msg.source_name = 'turtlesim'
        # Add some uncertainty based on velocity and distance from known points
        velocity_magnitude = math.sqrt(self.current_pose.linear_velocity**2 + self.current_pose.angular_velocity**2)
        base_uncertainty = 5.0
        velocity_uncertainty = velocity_magnitude * 2.0
        spatial_msg.uncertainty = base_uncertainty + velocity_uncertainty
        
        
        # Set timestamp
        spatial_msg.time = self.get_clock().now().to_msg()
        
        # Publish the measurement
        self.spatial_pub.publish(spatial_msg)
        
        self.get_logger().info(
            f'Published spatial measurement: pos=({self.current_pose.x:.2f}, {self.current_pose.y:.2f}), '
            f'value={spatial_msg.value:.2f}, uncertainty={spatial_msg.uncertainty:.2f}'
        )


def main(args=None):
    rclpy.init(args=args)
    
    node = TurtlesimSpatialPublisher()
    
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    
    node.destroy_node()
    rclpy.shutdown()


if __name__ == '__main__':
    main()
