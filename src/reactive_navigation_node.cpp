#include "polygeom_lib.h"
#include "reactive_planner_lib.h"
#include "safe_bayesian_optimization/msg/polygon_array.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/strategies/buffer.hpp>
#include <geometry_msgs/msg/polygon.hpp>
#include <rclcpp/rclcpp.hpp>

namespace bg = boost::geometry;

class ReactiveNavigationNode : public rclcpp::Node {
public:
  ReactiveNavigationNode() : Node("reactive_navigation_node") {
    // Declare reactive planner parameters
    this->declare_parameter("reactive_planner.p", 2.0);
    this->declare_parameter("reactive_planner.epsilon", 0.1);
    this->declare_parameter("reactive_planner.varepsilon", 0.05);
    this->declare_parameter("reactive_planner.mu_1", 1.0);
    this->declare_parameter("reactive_planner.mu_2", 1.0);
    this->declare_parameter("reactive_planner.robot_radius", 0.3);

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

    // Set default workspace (will be updated when envelope is received)
    std::vector<std::vector<double>> default_workspace = {
        {-10.0, -10.0}, {10.0, -10.0}, {10.0, 10.0}, {-10.0, 10.0}};

    // Configure diffeomorphism parameters
    diffeo_params_.set_all_params(p, epsilon, varepsilon, mu_1, mu_2,
                                  default_workspace);

    RCLCPP_INFO(this->get_logger(),
                "Configured reactive planner with p=%.2f, epsilon=%.3f, "
                "varepsilon=%.3f, mu_1=%.2f, mu_2=%.2f, robot_radius=%.3f",
                p, epsilon, varepsilon, mu_1, mu_2, robot_radius_);

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

    RCLCPP_INFO(this->get_logger(), "Reactive Navigation Node initialized");
  }

private:
  rclcpp::Subscription<safe_bayesian_optimization::msg::PolygonArray>::SharedPtr
      obstacles_sub_;
  rclcpp::Subscription<geometry_msgs::msg::Polygon>::SharedPtr envelope_sub_;
  DiffeoParamsClass diffeo_params_;
  double robot_radius_;
  std::vector<polygon> obstacle_polygons_;
  std::vector<std::vector<PolygonClass>> diffeo_tree_array_;

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

      auto merged_polygons = get_merged_dilated_polygons();

      for (const auto &merged_poly : merged_polygons) {
        std::vector<PolygonClass> tree;
        diffeoTreeConvex(BoostPointToStd(BoostPolyToBoostPoint(merged_poly)),
                         diffeo_params_, &tree);
        diffeo_tree_array_.push_back(tree);
      }

      RCLCPP_INFO(this->get_logger(),
                  "Converted %zu obstacle polygons to boost geometry format",
                  obstacle_polygons_.size());
    }
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
  }
};

int main(int argc, char **argv) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<ReactiveNavigationNode>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
