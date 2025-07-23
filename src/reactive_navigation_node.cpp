#include "polygeom_lib.h"
#include "reactive_planner_lib.h"
#include "safe_bayesian_optimization/msg/polygon_array.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/strategies/buffer.hpp>
#include <geometry_msgs/msg/polygon.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <rclcpp/rclcpp.hpp>
#include <tf2/LinearMath/Matrix3x3.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include <turtlesim/msg/pose.hpp>
#include <voro++/cell.hh>
#include <voro++/container.hh>
#include <voro++/voro++.hh>

namespace bg = boost::geometry;

class ReactiveNavigationNode : public rclcpp::Node {
public:
  ReactiveNavigationNode()
      : Node("reactive_navigation_node"), tf_buffer_(this->get_clock()),
        tf_listener_(tf_buffer_) {
    // Declare reactive planner parameters
    this->declare_parameter("reactive_planner.p", 2.0);
    this->declare_parameter("reactive_planner.epsilon", 0.1);
    this->declare_parameter("reactive_planner.varepsilon", 0.05);
    this->declare_parameter("reactive_planner.mu_1", 1.0);
    this->declare_parameter("reactive_planner.mu_2", 1.0);
    this->declare_parameter("reactive_planner.robot_radius", 0.3);
    this->declare_parameter("reactive_planner.linear_gain", 1.0);
    this->declare_parameter("reactive_planner.angular_gain", 1.0);
    this->declare_parameter("reactive_planner.linear_cmd_limit", 0.5);
    this->declare_parameter("reactive_planner.angular_cmd_limit", 1.0);

    // Read parameters and configure diffeomorphism parameters
    double p = this->get_parameter("reactive_planner.p").as_double();
    double epsilon =
        this->get_parameter("reactive_planner.epsilon").as_double();
    double varepsilon =
        this->get_parameter("reactive_planner.varepsilon").as_double();
    double mu_1 = this->get_parameter("reactive_planner.mu_1").as_double();
    double mu_2 = this->get_parameter("reactive_planner.mu_2").as_double();
    robot_radius_ =
        this->get_parameter("reactive_planner.robot_radius").as_double();
    linear_gain_ =
        this->get_parameter("reactive_planner.linear_gain").as_double();
    angular_gain_ =
        this->get_parameter("reactive_planner.angular_gain").as_double();
    linear_cmd_limit_ =
        this->get_parameter("reactive_planner.linear_cmd_limit").as_double();
    angular_cmd_limit_ =
        this->get_parameter("reactive_planner.angular_cmd_limit").as_double();

    // Set default workspace (will be updated when envelope is received)
    std::vector<std::vector<double>> default_workspace = {
        {-10.0, -10.0}, {10.0, -10.0}, {10.0, 10.0}, {-10.0, 10.0}};

    // Configure diffeomorphism parameters
    diffeo_params_.set_all_params(p, epsilon, varepsilon, mu_1, mu_2,
                                  default_workspace);

    RCLCPP_INFO(this->get_logger(),
                "Configured reactive planner with p=%.2f, epsilon=%.3f, "
                "varepsilon=%.3f, mu_1=%.2f, mu_2=%.2f, robot_radius=%.3f, "
                "linear_gain=%.2f, angular_gain=%.2f, linear_limit=%.2f, "
                "angular_limit=%.2f",
                p, epsilon, varepsilon, mu_1, mu_2, robot_radius_, linear_gain_,
                angular_gain_, linear_cmd_limit_, angular_cmd_limit_);

    // Subscribe to obstacle polygons
    obstacles_sub_ = this->create_subscription<
        safe_bayesian_optimization::msg::PolygonArray>(
        "polygon_array", 10,
        std::bind(&ReactiveNavigationNode::obstacles_callback, this,
                  std::placeholders::_1));

    // Subscribe to envelope polygon
    envelope_sub_ = this->create_subscription<geometry_msgs::msg::Polygon>(
        "envelope_polygon", 10,
        std::bind(&ReactiveNavigationNode::envelope_callback, this,
                  std::placeholders::_1));

    // Subscribe to turtlesim pose
    pose_sub_ = this->create_subscription<turtlesim::msg::Pose>(
        "/turtle1/pose", 10,
        std::bind(&ReactiveNavigationNode::pose_callback, this,
                  std::placeholders::_1));

    // Subscribe to current subgoal
    subgoal_sub_ = this->create_subscription<geometry_msgs::msg::Point>(
        "/current_subgoal", 10,
        std::bind(&ReactiveNavigationNode::subgoal_callback, this,
                  std::placeholders::_1));

    // Create control command publisher
    cmd_vel_pub_ = this->create_publisher<geometry_msgs::msg::Twist>(
        "/turtle1/cmd_vel", 10);

    // Create timer for control loop
    control_timer_ = this->create_wall_timer(
        std::chrono::milliseconds(500),
        std::bind(&ReactiveNavigationNode::control_callback, this));

    RCLCPP_INFO(this->get_logger(), "Reactive Navigation Node initialized");
  }

private:
  rclcpp::Subscription<safe_bayesian_optimization::msg::PolygonArray>::SharedPtr
      obstacles_sub_;
  rclcpp::Subscription<geometry_msgs::msg::Polygon>::SharedPtr envelope_sub_;
  rclcpp::Subscription<turtlesim::msg::Pose>::SharedPtr pose_sub_;
  rclcpp::Subscription<geometry_msgs::msg::Point>::SharedPtr subgoal_sub_;
  rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel_pub_;
  rclcpp::TimerBase::SharedPtr control_timer_;

  tf2_ros::Buffer tf_buffer_;
  tf2_ros::TransformListener tf_listener_;

  DiffeoParamsClass diffeo_params_;
  double robot_radius_;
  double linear_gain_;
  double angular_gain_;
  double linear_cmd_limit_;
  double angular_cmd_limit_;
  std::vector<polygon> obstacle_polygons_;
  std::vector<std::vector<PolygonClass>> diffeo_tree_array_;

  // Robot state
  geometry_msgs::msg::Point current_position_;
  double current_yaw_;
  bool has_robot_pose_ = false;

  // Subgoal state
  geometry_msgs::msg::Point current_subgoal_;
  bool has_subgoal_ = false;
  int env_x_min_ = 0;
  int env_x_max_ = 0;
  int env_y_min_ = 0;
  int env_y_max_ = 0;

  std::vector<polygon> get_merged_dilated_polygons() {
    if (obstacle_polygons_.empty()) {
      return {};
    }

    std::vector<polygon> polygon_list;
    polygon_list.reserve(obstacle_polygons_.size());

    // Dilate each polygon by robot radius using full buffer strategies
    bg::strategy::buffer::distance_symmetric<double> distance_strategy(
        robot_radius_);
    bg::strategy::buffer::join_round join_strategy;
    bg::strategy::buffer::end_round end_strategy;
    bg::strategy::buffer::point_circle point_strategy;
    bg::strategy::buffer::side_straight side_strategy;

    for (const auto &poly : obstacle_polygons_) {
      multi_polygon dilated_result;
      boost::geometry::buffer(poly, dilated_result, distance_strategy,
                              side_strategy, join_strategy, end_strategy,
                              point_strategy);

      if (!dilated_result.empty()) {
        polygon_list.push_back(dilated_result.front());
      }
    }

    RCLCPP_INFO(this->get_logger(), "Dilated %zu polygons",
                polygon_list.size());

    if (polygon_list.empty()) {
      return {};
    }

    // Merge all overlapping polygons using union operation
    multi_polygon output_union;
    if (polygon_list.size() >= 1) {
      output_union.push_back(polygon_list.back());
      polygon_list.pop_back();

      while (!polygon_list.empty()) {
        polygon next_polygon = polygon_list.back();
        polygon_list.pop_back();
        multi_polygon temp_result;
        bg::union_(output_union, next_polygon, temp_result);
        output_union = temp_result;
      }
    }

    // Simplify and create final merged polygon list
    std::vector<polygon> polygon_list_merged;
    for (size_t i = 0; i < output_union.size(); i++) {
      polygon simplified_component;
      bg::simplify(output_union[i], simplified_component, 0.2);
      polygon_list_merged.push_back(simplified_component);
    }

    RCLCPP_INFO(this->get_logger(),
                "Merged %zu dilated polygons into %zu merged polygons",
                obstacle_polygons_.size(), polygon_list_merged.size());

    return polygon_list_merged;
  }

  void obstacles_callback(
      const safe_bayesian_optimization::msg::PolygonArray::SharedPtr msg) {
    RCLCPP_INFO(this->get_logger(), "Received %zu obstacle polygons",
                msg->polygons.size());

    // Clear previous obstacles
    obstacle_polygons_.clear();
    obstacle_polygons_.reserve(msg->polygons.size());

    // Convert ROS polygons to boost geometry polygons
    for (const auto &ros_polygon : msg->polygons) {
      if (ros_polygon.points.size() < 3) {
        RCLCPP_WARN(this->get_logger(),
                    "Skipping polygon with less than 3 points");
        continue;
      }

      polygon boost_poly;

      // Add points to the polygon
      for (const auto &point : ros_polygon.points) {
        boost::geometry::append(
            boost_poly.outer(),
            bg::model::point<double, 2, bg::cs::cartesian>(
                static_cast<double>(point.x), static_cast<double>(point.y)));
      }

      // Close the polygon by adding the first point again if needed
      if (!boost::geometry::equals(boost_poly.outer().front(),
                                   boost_poly.outer().back())) {
        boost::geometry::append(boost_poly.outer(), boost_poly.outer().front());
      }

      // Correct the polygon (ensure proper orientation and closure)
      boost::geometry::correct(boost_poly);

      obstacle_polygons_.push_back(boost_poly);
    }
    RCLCPP_INFO(this->get_logger(),
                "Converted %zu obstacle polygons to boost geometry format",
                obstacle_polygons_.size());

    auto merged_polygons = get_merged_dilated_polygons();

    for (const auto &merged_poly : merged_polygons) {
      std::vector<PolygonClass> tree;
      diffeoTreeConvex(BoostPointToStd(BoostPolyToBoostPoint(merged_poly)),
                       diffeo_params_, &tree);
      diffeo_tree_array_.push_back(tree);
    }
    RCLCPP_INFO(this->get_logger(),
                "Converted %zu merged polygons to diffeomorphism trees",
                diffeo_tree_array_.size());
  }

  void envelope_callback(const geometry_msgs::msg::Polygon::SharedPtr msg) {
    // Convert ROS polygon to workspace format
    std::vector<std::vector<double>> workspace;
    workspace.reserve(msg->points.size());

    for (const auto &point : msg->points) {
      std::vector<double> vertex = {static_cast<double>(point.x),
                                    static_cast<double>(point.y)};
      workspace.push_back(vertex);
    }

    // Update workspace in diffeomorphism parameters
    diffeo_params_.set_workspace(workspace);
    env_x_min_ = workspace[0][0];
    env_x_max_ = workspace[2][0];
    env_y_min_ = workspace[0][1];
    env_y_max_ = workspace[2][1];
  }

  void subgoal_callback(const geometry_msgs::msg::Point::SharedPtr msg) {
    current_subgoal_ = *msg;
    has_subgoal_ = true;

    RCLCPP_INFO(this->get_logger(), "Received subgoal: x=%.3f, y=%.3f, z=%.3f",
                current_subgoal_.x, current_subgoal_.y, current_subgoal_.z);
  }

  void pose_callback(const turtlesim::msg::Pose::SharedPtr msg) {
    // Update current position directly from turtlesim pose
    current_position_.x = msg->x;
    current_position_.y = msg->y;
    current_position_.z = 0.0; // Turtlesim is 2D

    // Use theta directly as yaw
    current_yaw_ = msg->theta;

    has_robot_pose_ = true;
  }

  void control_callback() {
    if (!has_robot_pose_) {
      RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 1000,
                           "No robot pose available for control");
      return;
    }

    if (obstacle_polygons_.empty()) {
      RCLCPP_INFO(this->get_logger(),
                  "No obstacles detected, skipping control");
      return;
    }

    geometry_msgs::msg::Twist cmd_vel;
    cmd_vel.linear.x = 0.0;
    cmd_vel.linear.y = 0.0;
    cmd_vel.angular.z = 0.0;

    std::vector<double> robot_position = {current_position_.x,
                                          current_position_.y};

    DiffeoTransformResult transform_result = computeDiffeoTransform(
        robot_position, current_yaw_, diffeo_tree_array_, diffeo_params_);

    std::vector<point> model_obstacle_centers;
    std::vector<double> model_obstacle_radii;

    for (const auto &tree : diffeo_tree_array_) {
      if (tree.empty()) {
        RCLCPP_WARN(this->get_logger(), "Diffeo tree is empty, skipping");
        continue;
      }

      point root_center = tree.back().get_center();
      double root_radius = tree.back().get_radius();

      model_obstacle_centers.push_back(root_center);
      model_obstacle_radii.push_back(root_radius);
    }
    RCLCPP_INFO(this->get_logger(), "computing local workspace polygon");

    polygon local_workspace_polygon = compute_local_workspace_polygon(
        transform_result.transformed_position, model_obstacle_centers,
        model_obstacle_radii);

    polygon local_free_space_polygon;
    if (local_workspace_polygon.outer().size() < 3) {
      RCLCPP_WARN(this->get_logger(),
                  "Local workspace polygon has less than 3 points, skipping "
                  "control");
      return;
    }
    bg::strategy::buffer::distance_symmetric<double> distance_strategy(
        robot_radius_);
    bg::strategy::buffer::join_round join_strategy;
    bg::strategy::buffer::end_round end_strategy;
    bg::strategy::buffer::point_circle point_strategy;
    bg::strategy::buffer::side_straight side_strategy;

    // Dilate the local workspace polygon to create free space
    multi_polygon local_free_space_multi;
    bg::buffer(local_workspace_polygon, local_free_space_multi,
               distance_strategy, side_strategy, join_strategy, end_strategy,
               point_strategy);

    // Extract the first polygon if buffer operation succeeded
    if (!local_free_space_multi.empty()) {
      local_free_space_polygon = local_free_space_multi[0];
    } else {
      RCLCPP_WARN(this->get_logger(),
                  "Buffer operation failed, using original polygon");
      local_free_space_polygon = local_workspace_polygon;
    }
    point robot_position_point(transform_result.transformed_position[0],
                               transform_result.transformed_position[1]);

    if (!has_subgoal_) {
      RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 1000,
                           "No subgoal available for control");
      return;
    }

    point subgoal_point(current_subgoal_.x, current_subgoal_.y);
    point local_goal_linear = compute_local_goal_linear(
        robot_position_point, transform_result.transformed_orientation,
        local_free_space_polygon, subgoal_point);

    point local_goal_angular =
        compute_local_goal(local_free_space_polygon, subgoal_point);

    // Compute the basis for the virtual control inputs
    double tV = (local_goal_linear.get<0>() -
                 transform_result.transformed_position[0]) *
                    cos(transform_result.transformed_orientation) +
                (local_goal_linear.get<1>() -
                 transform_result.transformed_position[1]) *
                    sin(transform_result.transformed_orientation);
    double tW1 = (local_goal_angular.get<0>() -
                  transform_result.transformed_position[0]) *
                     cos(transform_result.transformed_orientation) +
                 (local_goal_angular.get<1>() -
                  transform_result.transformed_position[1]) *
                     sin(transform_result.transformed_orientation);
    double tW2 = -(local_goal_angular.get<0>() -
                   transform_result.transformed_position[0]) *
                     sin(transform_result.transformed_orientation) +
                 (local_goal_angular.get<1>() -
                  transform_result.transformed_position[1]) *
                     cos(transform_result.transformed_orientation);

    // Compute the basis for transforming to actual control inputs

    double alpha1 = transform_result.alpha1;
    double alpha2 = transform_result.alpha2;
    double beta1 = transform_result.beta1;
    double beta2 = transform_result.beta2;

    double e_norm = sqrt(pow(alpha1, 2) + pow(alpha2, 2));
    double dksi_dpsi =
        MatrixDeterminant(transform_result.transformed_jacobian) /
        pow(e_norm, 2);
    double DksiCosSin = (alpha1 * beta1 + alpha2 * beta2) / pow(e_norm, 2);

    double linear_ctl_gain, angular_ctl_gain;
    std::vector<double> limit_check_vector_linear = {
        linear_gain_, linear_cmd_limit_ * e_norm / std::abs(tV),
        0.4 * angular_cmd_limit_ * dksi_dpsi * e_norm /
            std::abs(tV * DksiCosSin)};
    linear_ctl_gain = *std::min_element(limit_check_vector_linear.begin(),
                                        limit_check_vector_linear.end());
    std::vector<double> limit_check_vector_angular = {
        angular_gain_,
        0.6 * angular_cmd_limit_ * dksi_dpsi / std::abs(std::atan2(tW2, tW1))};

    angular_ctl_gain = *std::min_element(limit_check_vector_angular.begin(),
                                         limit_check_vector_angular.end());

    // Compute virtual and actual inputs
    double dV_virtual = linear_ctl_gain * tV;
    double linear_cmd = dV_virtual / e_norm;
    double dW_virtual = angular_ctl_gain * std::atan2(tW2, tW1);
    double angular_cmd = (dW_virtual - linear_cmd * DksiCosSin) / dksi_dpsi;

    // Publish control commands
    cmd_vel.linear.x =
        std::max(-linear_cmd_limit_, std::min(linear_cmd, linear_cmd_limit_));
    cmd_vel.angular.z = std::max(-angular_cmd_limit_,
                                 std::min(angular_cmd, angular_cmd_limit_));
    cmd_vel_pub_->publish(cmd_vel);

    RCLCPP_INFO(this->get_logger(),
                "Published cmd_vel: linear.x=%.3f, angular.z=%.3f",
                cmd_vel.linear.x, cmd_vel.angular.z);
  }

  polygon compute_local_workspace_polygon(
      std::vector<double> robot_position_transformed,
      std::vector<point> &model_obstacle_centers,
      std::vector<double> &model_obstacle_radii) {

    int num_blocks_x = env_x_max_ - env_x_min_;
    int num_blocks_y = env_y_max_ - env_y_min_;

    voro::container_poly con(env_x_min_, env_x_max_, env_y_min_, env_y_max_,
                             0.5, 0.5, num_blocks_x, num_blocks_y, 1, false,
                             false, false, 4);
    con.put(0, robot_position_transformed.at(0),
            robot_position_transformed.at(1), 0.0, robot_radius_);
    for (size_t i = 0; i < model_obstacle_centers.size(); ++i) {
      const point &center = model_obstacle_centers[i];
      double radius = model_obstacle_radii[i];
      con.put(i + 1, center.get<0>(), center.get<1>(), 0.0, radius);
    }

    voro::voronoicell_neighbor c;
    polygon local_workspace_polygon;

    if (con.compute_ghost_cell(c, robot_position_transformed.at(0),
                               robot_position_transformed.at(1), 0,
                               robot_radius_)) {
      // Extract the vertices of the ghost cell
      std::vector<double> vertices;
      c.vertices(vertices);
      // find all vertices that have z = 0.5 (they are either -0.5 or 0.5)
      for (size_t i = 0; i < vertices.size(); i += 3) {
        if (vertices[i + 2] > 0) {
          point p(vertices[i], vertices[i + 1]);
          local_workspace_polygon.outer().push_back(p);
        }
      }
    }
    return local_workspace_polygon;
  }

  line compute_linear_free_space(point robot_position, double robot_orientation,
                                 polygon local_free_space_polygon) {

    line free_space_line;
    if (bg::area(local_free_space_polygon) < 0.01) {
      free_space_line.push_back(robot_position);
      free_space_line.push_back(robot_position);
      return free_space_line;
    } else {
      free_space_line.push_back(robot_position);
      free_space_line.push_back(
          polyxray(local_free_space_polygon, robot_position,
                   point(cos(robot_orientation), sin(robot_orientation))));
      return free_space_line;
    }
  }

  point compute_local_goal_linear(point robot_position,
                                  double robot_orientation,
                                  polygon local_free_space_polygon,
                                  point goal) {
    line lfl = compute_linear_free_space(robot_position, robot_orientation,
                                         local_free_space_polygon);
    if (bg::is_empty(lfl)) {
      return robot_position;
    } else {
      ProjectionResultStruct projection_result = linedist(lfl, goal);
      return projection_result.projected_point;
    }
  }

  point compute_local_goal(polygon lfs, point goal) {
    if (bg::area(lfs) < 0.01) {
      return goal;
    } else {
      if (bg::within(goal, lfs)) {
        return goal; // Goal is already within local free space
      } else {
        ProjectionResultStruct projection_result = polydist(lfs, goal);
        return projection_result.projected_point;
      }
    }
  }
};

int main(int argc, char **argv) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<ReactiveNavigationNode>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
