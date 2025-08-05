#include "polygeom_lib.h"
#include "safe_bayesian_optimization/msg/polygon_array.hpp"
#include "trusses_custom_interfaces/srv/get_terrain_map_with_uncertainty.hpp"
#include "trusses_custom_interfaces/srv/spatial_data.hpp"
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/algorithms/detail/envelope/interface.hpp>
#include <boost/geometry/algorithms/difference.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/buffer.hpp>
#include <cmath>
#include <cv_bridge/cv_bridge.h>
#include <geometry_msgs/msg/point.hpp>
#include <geometry_msgs/msg/point_stamped.hpp>
#include <geometry_msgs/msg/polygon.hpp>
#include <limits>
#include <memory>
#include <opencv2/opencv.hpp>
#include <unordered_map>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include <ratio>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/image.hpp>

#include <vector>

namespace bg = boost::geometry;

class OptimizerNode : public rclcpp::Node {
public:
  OptimizerNode() : Node("optimizer_node") {
    // Declare parameters
    this->declare_parameter("opt.beta", 2.0);
    this->declare_parameter("opt.f_min", 0.0);
    this->declare_parameter("terrain_map.resolution",
                            std::vector<double>{0.1, 0.1});
    this->declare_parameter("debug.publish_debug_image", false);
    this->declare_parameter("opt.subgoal_erosion", 0.2);

    // Read parameters
    beta_ = this->get_parameter("opt.beta").as_double();
    f_min_ = this->get_parameter("opt.f_min").as_double();
    auto resolution_vector =
        this->get_parameter("terrain_map.resolution").as_double_array();
    if (resolution_vector.size() != 2) {
      RCLCPP_ERROR(this->get_logger(),
                   "terrain_map.resolution must be an array of 2 values "
                   "[x_resolution, y_resolution]");
      throw std::runtime_error("Invalid terrain_map.resolution parameter");
    }
    terrain_x_resolution_ = resolution_vector[0];
    terrain_y_resolution_ = resolution_vector[1];
    publish_debug_image_ =
        this->get_parameter("debug.publish_debug_image").as_bool();
    subgoal_erosion_ = this->get_parameter("opt.subgoal_erosion").as_double();

    // Create service clients
    spatial_data_client_ =
        this->create_client<trusses_custom_interfaces::srv::SpatialData>(
            "get_spatial_data");
    terrain_map_client_ = this->create_client<
        trusses_custom_interfaces::srv::GetTerrainMapWithUncertainty>(
        "get_terrain_map_with_uncertainty");

    // Create subscriber
    goal_point_sub_ = this->create_subscription<geometry_msgs::msg::PointStamped>(
        "goal_point", 10,
        std::bind(&OptimizerNode::goal_point_callback, this,
                  std::placeholders::_1));

    // Create publishers
    current_subgoal_pub_ = this->create_publisher<geometry_msgs::msg::Point>(
        "current_subgoal", 10);

    subgoal_marker_pub_ =
        this->create_publisher<visualization_msgs::msg::Marker>(
            "subgoal_marker", 10);
    polygons_pub_ =
        this->create_publisher<safe_bayesian_optimization::msg::PolygonArray>(
            "polygon_array", 10);

    envelope_pub_ = this->create_publisher<geometry_msgs::msg::Polygon>(
        "envelope_polygon", 10);

    // Create concave polygon markers publisher
    concave_markers_pub_ =
        this->create_publisher<visualization_msgs::msg::MarkerArray>(
            "/concave_markers", 10);

    // Create projected goal marker publisher
    projected_goal_marker_pub_ =
        this->create_publisher<visualization_msgs::msg::Marker>(
            "/projected_goal_marker", 10);

    // Only create debug image publisher if enabled
    if (publish_debug_image_) {
      debug_image_pub_ = this->create_publisher<sensor_msgs::msg::Image>(
          "debug_polygons_image", 10);
    }

    // Create timer to check spatial data periodically
    spatial_data_timer_ = this->create_wall_timer(
        std::chrono::seconds(5),
        std::bind(&OptimizerNode::check_spatial_data, this));

    RCLCPP_INFO(this->get_logger(), "Optimizer Node Initialized");
  }

private:
  Eigen::VectorXd mu_;  // Mean vector
  Eigen::VectorXd std_; // Standard deviation vector

  Eigen::Matrix<double, Eigen::Dynamic, 2> D_; // Parameter Set;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> S_;   // Safe set
  Eigen::Matrix<double, Eigen::Dynamic, 2> Q_; // Confidence intervals

  double beta_;
  double f_min_;
  double subgoal_erosion_;

  // Spatial data monitoring
  rclcpp::Client<trusses_custom_interfaces::srv::SpatialData>::SharedPtr
      spatial_data_client_;
  rclcpp::Client<trusses_custom_interfaces::srv::GetTerrainMapWithUncertainty>::
      SharedPtr terrain_map_client_;
  rclcpp::TimerBase::SharedPtr spatial_data_timer_;
  rclcpp::Subscription<geometry_msgs::msg::PointStamped>::SharedPtr goal_point_sub_;
  rclcpp::Publisher<geometry_msgs::msg::Point>::SharedPtr current_subgoal_pub_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr
      subgoal_marker_pub_;
  rclcpp::Publisher<safe_bayesian_optimization::msg::PolygonArray>::SharedPtr
      polygons_pub_;
  rclcpp::Publisher<geometry_msgs::msg::Polygon>::SharedPtr envelope_pub_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr
      concave_markers_pub_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr
      projected_goal_marker_pub_;
  rclcpp::Publisher<sensor_msgs::msg::Image>::SharedPtr debug_image_pub_;

  // Parameters from config
  double terrain_x_resolution_;
  double terrain_y_resolution_;
  bool publish_debug_image_;

  // Terrain dimensions from service response
  int terrain_width_cells_;
  int terrain_height_cells_;

  // Current goal point
  geometry_msgs::msg::PointStamped current_goal_;

  // Store eroded concave polygon for subgoal projection
  bg::model::polygon<bg::model::d2::point_xy<double>> eroded_concave_polygon_;

  void ComputeSets() {
    // Compute the confidence intervals
    ComputeConfidenceIntervals();

    // Debug print mu std, q and fmin

    // Update the safe set based on the confidence intervals
    UpdateSafeSet();
  }

  void UpdateSafeSet() { S_ = Q_.col(0).array() > f_min_; }

  void ComputeConfidenceIntervals() {
    const Eigen::VectorXd confidence = beta_ * std_;

    Q_.col(0) = mu_ - confidence;
    Q_.col(1) = mu_ + confidence;
  }

  std::vector<int> FindSafetyContourIndices() {
    if (D_.rows() == 0 || S_.rows() == 0) {
      RCLCPP_WARN(this->get_logger(), "D_ or S_ is empty");
      return {};
    }

    // count S_ that is true

    int safe_count = S_.count();
    RCLCPP_INFO(this->get_logger(), "Number of safe points: %d", safe_count);
    RCLCPP_INFO(this->get_logger(), "Total points in S_: %ld", S_.size());

    // Find grid bounds
    int min_x = static_cast<int>(D_.col(0).minCoeff());
    int max_x = static_cast<int>(D_.col(0).maxCoeff());
    int min_y = static_cast<int>(D_.col(1).minCoeff());
    int max_y = static_cast<int>(D_.col(1).maxCoeff());

    // Use stored grid dimensions from service response
    int width = terrain_width_cells_;
    int height = terrain_height_cells_;

    RCLCPP_INFO(this->get_logger(), "Using grid dimensions: %dx%d", width,
                height);

    // Create binary image from safety data
    cv::Mat safety_image = cv::Mat::zeros(height, width, CV_8UC1);

    for (int i = 0; i < D_.rows(); ++i) {
      // Get indices in the grid from D_, where the coordinates are not
      // necessarily integers. so we need to scale them

      int x = static_cast<int>((D_(i, 0) - min_x) / (max_x - min_x) * width);
      int y = static_cast<int>((D_(i, 1) - min_y) / (max_y - min_y) * height);

      if (x >= 0 && x < width && y >= 0 && y < height) {
        safety_image.at<uchar>(y, x) = S_(i) ? 255 : 0;
      }
    }

    // Find contours
    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarchy;
    cv::findContours(safety_image, contours, hierarchy, cv::RETR_EXTERNAL,
                     cv::CHAIN_APPROX_NONE);

    RCLCPP_INFO(this->get_logger(), "Found %lu contours", contours.size());

    // Convert contour points back to original indices using a hash map for O(1)
    // lookup
    std::unordered_map<int, std::unordered_map<int, int>> coord_to_index;
    for (int i = 0; i < D_.rows(); ++i) {
      int x = static_cast<int>((D_(i, 0) - min_x) / (max_x - min_x) * width);
      int y = static_cast<int>((D_(i, 1) - min_y) / (max_y - min_y) * height);
      if (x >= 0 && x < width && y >= 0 && y < height) {
        coord_to_index[x][y] = i;
      }
    }

    std::vector<int> contour_indices;
    contour_indices.reserve(contours.size() *
                            50); // Reserve reasonable capacity

    for (const auto &contour : contours) {
      for (const auto &point : contour) {
        // O(1) lookup using scaled coordinates
        auto x_it = coord_to_index.find(point.x);
        if (x_it != coord_to_index.end()) {
          auto y_it = x_it->second.find(point.y);
          if (y_it != x_it->second.end()) {
            contour_indices.push_back(y_it->second);
          }
        }
      }
    }

    RCLCPP_INFO(this->get_logger(), "Found %lu contour points",
                contour_indices.size());
    return contour_indices;
  }

  int GetNextSubgoal() {
    // Get frontier indices
    const std::vector<int> frontier_indices = FindSafetyContourIndices();

    if (frontier_indices.empty()) {
      RCLCPP_WARN(this->get_logger(), "No frontier points found");
      return -1;
    }

    // Extract frontier points
    const size_t num_frontiers = frontier_indices.size();
    Eigen::MatrixXd frontier_points(num_frontiers, 2);
    Eigen::VectorXd frontier_confidence_width(num_frontiers);

    for (size_t i = 0; i < num_frontiers; ++i) {
      const int idx = frontier_indices[i];
      frontier_points.row(i) = D_.row(idx);
      frontier_confidence_width(i) = Q_(idx, 1) - Q_(idx, 0);
    }

    const Eigen::VectorXd goal_eigen =
        (Eigen::VectorXd(2) << current_goal_.point.x, current_goal_.point.y).finished();
    const Eigen::VectorXd distances =
        (frontier_points.rowwise() - goal_eigen.transpose()).rowwise().norm();

    const Eigen::VectorXd scores =
        distances.array() / (frontier_confidence_width.array() + 1e-6);
    int best_index;
    scores.minCoeff(&best_index);

    return frontier_indices[best_index];
  }

  void check_spatial_data() {
    using namespace std::chrono_literals;

    if (!spatial_data_client_->wait_for_service(0s)) {
      RCLCPP_WARN(this->get_logger(),
                  "Spatial data service 'get_spatial_data' not available");
      return;
    }

    auto request = std::make_shared<
        trusses_custom_interfaces::srv::SpatialData::Request>();

    // Use callback-based async call
    spatial_data_client_->async_send_request(
        request,
        [this](rclcpp::Client<trusses_custom_interfaces::srv::SpatialData>::
                   SharedFuture future) {
          auto response = future.get();
          if (response->success) {
            RCLCPP_INFO(this->get_logger(),
                        "Received spatial data with %zu points",
                        response->values.size());

            // Process the spatial data here
            // response->x_coords, response->y_coords, response->values

            // Request terrain map with resolution (service will calculate
            // dimensions from data bounds)
            request_terrain_map();

          } else {
            RCLCPP_WARN(this->get_logger(), "Spatial data request failed: %s",
                        response->message.c_str());
          }
        });
  }

  void request_terrain_map() {
    if (!terrain_map_client_->wait_for_service(std::chrono::seconds(0))) {
      RCLCPP_WARN(this->get_logger(), "Terrain map service not available");
      return;
    }

    RCLCPP_INFO(this->get_logger(), "Requesting terrain map...");

    auto request =
        std::make_shared<trusses_custom_interfaces::srv::
                             GetTerrainMapWithUncertainty::Request>();
    // Use the resolution array field from the service interface
    request->resolution = {static_cast<float>(terrain_x_resolution_),
                           static_cast<float>(terrain_y_resolution_)};

    RCLCPP_INFO(this->get_logger(),
                "Requesting terrain map with resolution: [%.2f, %.2f] "
                "(service will calculate dimensions from data bounds)",
                terrain_x_resolution_, terrain_y_resolution_);

    // Use callback-based async call
    terrain_map_client_->async_send_request(
        request,
        [this](rclcpp::Client<
               trusses_custom_interfaces::srv::GetTerrainMapWithUncertainty>::
                   SharedFuture future) {
          auto response = future.get();
          if (response->success) {
            RCLCPP_INFO(this->get_logger(), "Received terrain map: %s",
                        response->message.c_str());

            // Store grid dimensions for later use
            terrain_width_cells_ = response->n_width_cells;
            terrain_height_cells_ = response->n_height_cells;

            // Process the terrain map data here
            process_terrain_map(response);

          } else {
            RCLCPP_WARN(this->get_logger(), "Terrain map request failed: %s",
                        response->message.c_str());
          }
        });
  }

  void process_terrain_map(
      const std::shared_ptr<const trusses_custom_interfaces::srv::
                                GetTerrainMapWithUncertainty::Response>
          response) {
    // Update the parameter set D_, mean mu_, and std_ with the new terrain map
    // data
    const size_t num_points = response->x_coords.size();

    D_.resize(num_points, 2);
    mu_.resize(num_points);
    std_.resize(num_points);
    Q_.resize(num_points, 2);
    S_.resize(num_points);

    for (size_t i = 0; i < num_points; ++i) {
      D_(i, 0) = response->x_coords[i];
      D_(i, 1) = response->y_coords[i];
      mu_(i) = response->values[i];
      std_(i) = response->uncertainties[i];
    }

    // Recompute sets with new data
    ComputeSets();

    publish_obstacle_polygons();

    // Check if the goal itself is safe by checking if it's within the eroded
    // polygon
    bg::model::d2::point_xy<double> goal_point(current_goal_.point.x,
                                               current_goal_.point.y);

    geometry_msgs::msg::Point subgoal;
    if (bg::within(goal_point, eroded_concave_polygon_)) {
      // Goal is safe, use it directly as the subgoal
      subgoal.x = current_goal_.point.x;
      subgoal.y = current_goal_.point.y;
      subgoal.z = 0.0;

      RCLCPP_INFO(
          this->get_logger(),
          "Goal is within eroded safe region, using goal as subgoal: (%f, %f)",
          subgoal.x, subgoal.y);
    } else {
      // Goal is not safe, find frontier subgoal and project it
      const int next_subgoal_index = GetNextSubgoal();
      if (next_subgoal_index >= 0) {
        // Get the raw subgoal from the terrain map
        double raw_x = response->x_coords[next_subgoal_index];
        double raw_y = response->y_coords[next_subgoal_index];

        // Convert boost polygon to polygeom_lib polygon format
        polygon eroded_poly_for_projection;
        for (const auto &boost_point : eroded_concave_polygon_.outer()) {
          point polygeom_point(boost_point.get<0>(), boost_point.get<1>());
          eroded_poly_for_projection.outer().push_back(polygeom_point);
        }

        bg::correct(eroded_poly_for_projection);

        RCLCPP_INFO(this->get_logger(),
                    "Next subgoal index: %d, coordinates: (%f, %f)",
                    next_subgoal_index, raw_x, raw_y);

        // Create point from raw subgoal coordinates
        point raw_subgoal_point(raw_x, raw_y);

        // Project the subgoal onto the eroded safe region using polydist
        ProjectionResultStruct projection_result =
            polydist(eroded_poly_for_projection, raw_subgoal_point);

        RCLCPP_INFO(this->get_logger(),
                    "Projected subgoal point: (%f, %f), distance: %f",
                    projection_result.projected_point.get<0>(),
                    projection_result.projected_point.get<1>(),
                    projection_result.dist);

        // Create final subgoal message from projected point
        subgoal.x = projection_result.projected_point.get<0>();
        subgoal.y = projection_result.projected_point.get<1>();
        subgoal.z = 0.0;

        // Publish visualization of the projection
        publish_projected_goal_marker(current_goal_.point, subgoal, projection_result.dist);

        RCLCPP_INFO(this->get_logger(),
                    "Published projected subgoal: (%f, %f) (distance: %f)",
                    subgoal.x, subgoal.y, projection_result.dist);
      } else {
        RCLCPP_WARN(this->get_logger(), "No valid subgoal found");
        return;
      }
    }

    current_subgoal_pub_->publish(subgoal);

    // Publish subgoal marker for visualization
    publish_subgoal_marker(subgoal);
  }

  void publish_obstacle_polygons() {

    // Use CGAL Alpha Shape instead of concaveman
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
    typedef CGAL::Alpha_shape_face_base_2<K> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds> Triangulation_2;
    typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
    typedef K::Point_2 Point_2;

    // Convert points to CGAL format - reserve capacity based on safe point
    // count
    int safe_count = S_.count();
    if (safe_count == 0) {
      RCLCPP_WARN(this->get_logger(),
                  "No safe points available for polygon generation");
      return;
    }

    std::vector<Point_2> cgal_points;
    cgal_points.reserve(safe_count);
    for (int i = 0; i < D_.rows(); ++i) {
      if (S_(i)) { // Only consider safe points
        cgal_points.emplace_back(D_(i, 0), D_(i, 1));
      }
    }

    // Create alpha shape
    Alpha_shape_2 alpha_shape(cgal_points.begin(), cgal_points.end());

    // Set alpha shape to regularized mode for better boundary extraction
    alpha_shape.set_mode(Alpha_shape_2::REGULARIZED);

    // Find optimal alpha value (can be adjusted based on your needs)
    // Using a small alpha to get a detailed boundary
    auto alpha_iterator = alpha_shape.find_optimal_alpha(1);
    if (alpha_iterator != alpha_shape.alpha_end()) {
      alpha_shape.set_alpha(*alpha_iterator);
    } else {
      // Fallback to a reasonable alpha value
      alpha_shape.set_alpha(0.1);
    }

    // Extract boundary edges and create ordered boundary
    std::vector<std::array<double, 2>> alpha_shape_boundary;
    alpha_shape_boundary.reserve(safe_count / 4); // Estimate boundary size

    // Collect boundary edges - reserve estimated capacity
    std::vector<typename Alpha_shape_2::Edge> boundary_edges;
    boundary_edges.reserve(safe_count / 2); // Conservative estimate
    for (auto edge_it = alpha_shape.alpha_shape_edges_begin();
         edge_it != alpha_shape.alpha_shape_edges_end(); ++edge_it) {
      boundary_edges.push_back(*edge_it);
    }

    // Create boost geometry polygon from alpha shape boundary
    bg::model::polygon<bg::model::d2::point_xy<double>> concave_polygon;

    if (!boundary_edges.empty()) {
      // Create a properly ordered boundary by traversing edges
      std::vector<Point_2> ordered_boundary;
      ordered_boundary.reserve(boundary_edges.size()); // Reserve capacity
      std::map<Point_2, std::vector<Point_2>> adjacency_map;

      // Build adjacency map from boundary edges
      for (const auto &edge : boundary_edges) {
        auto face = edge.first;
        int idx = edge.second;

        Point_2 p1 = face->vertex((idx + 1) % 3)->point();
        Point_2 p2 = face->vertex((idx + 2) % 3)->point();

        adjacency_map[p1].push_back(p2);
        adjacency_map[p2].push_back(p1);
      }

      // Traverse the boundary to create ordered polygon
      if (!adjacency_map.empty()) {
        std::set<Point_2> visited;
        Point_2 start_point = adjacency_map.begin()->first;
        Point_2 current = start_point;

        do {
          ordered_boundary.push_back(current);
          visited.insert(current);

          // Find next unvisited neighbor
          Point_2 next = current; // fallback
          for (const auto &neighbor : adjacency_map[current]) {
            if (visited.find(neighbor) == visited.end() ||
                (neighbor == start_point && ordered_boundary.size() > 2)) {
              next = neighbor;
              break;
            }
          }

          current = next;
        } while (current != start_point &&
                 visited.find(current) == visited.end() &&
                 ordered_boundary.size() < adjacency_map.size());
      }

      // Convert ordered boundary to boost geometry polygon
      for (const auto &point : ordered_boundary) {
        bg::append(concave_polygon.outer(),
                   bg::model::d2::point_xy<double>(CGAL::to_double(point.x()),
                                                   CGAL::to_double(point.y())));
        alpha_shape_boundary.push_back(
            {CGAL::to_double(point.x()), CGAL::to_double(point.y())});
      }
    } else {
      // Fallback: if no alpha shape boundary, use convex hull
      RCLCPP_WARN(this->get_logger(),
                  "No alpha shape boundary found, using all safe points");
      for (const auto &point : cgal_points) {
        bg::append(concave_polygon.outer(),
                   bg::model::d2::point_xy<double>(CGAL::to_double(point.x()),
                                                   CGAL::to_double(point.y())));
        alpha_shape_boundary.push_back(
            {CGAL::to_double(point.x()), CGAL::to_double(point.y())});
      }
    }

    bg::correct(concave_polygon);

    // Erode the concave polygon for safe subgoal placement
    bg::strategy::buffer::distance_symmetric<double> erosion_distance(
        -subgoal_erosion_);
    bg::strategy::buffer::join_round join_strategy;
    bg::strategy::buffer::end_round end_strategy;
    bg::strategy::buffer::point_circle point_strategy;
    bg::strategy::buffer::side_straight side_strategy;

    bg::model::multi_polygon<
        bg::model::polygon<bg::model::d2::point_xy<double>>>
        eroded_multi;
    bg::buffer(concave_polygon, eroded_multi, erosion_distance, side_strategy,
               join_strategy, end_strategy, point_strategy);

    // Store the eroded polygon (use the first polygon if buffer succeeded)
    if (!eroded_multi.empty()) {
      eroded_concave_polygon_ = eroded_multi[0];
      RCLCPP_INFO(this->get_logger(),
                  "Eroded concave polygon for subgoal projection");
    } else {
      // Fallback: use original polygon if erosion failed
      eroded_concave_polygon_ = concave_polygon;
      RCLCPP_WARN(this->get_logger(),
                  "Erosion failed, using original concave polygon");
    }

    // Publish concave polygon markers for visualization
    publish_concave_markers(concave_polygon);

    bg::model::box<bg::model::d2::point_xy<double>> envelope_box;
    bg::envelope(concave_polygon, envelope_box);

    // Publish envelope polygon
    geometry_msgs::msg::Polygon envelope_msg;
    geometry_msgs::msg::Point32 min_corner, max_corner;
    min_corner.x = envelope_box.min_corner().get<0>();
    min_corner.y = envelope_box.min_corner().get<1>();
    max_corner.x = envelope_box.max_corner().get<0>();
    max_corner.y = envelope_box.max_corner().get<1>();

    // Create rectangle from envelope box
    envelope_msg.points.push_back(min_corner);
    geometry_msgs::msg::Point32 bottom_right;
    bottom_right.x = max_corner.x;
    bottom_right.y = min_corner.y;
    envelope_msg.points.push_back(bottom_right);
    envelope_msg.points.push_back(max_corner);
    geometry_msgs::msg::Point32 top_left;
    top_left.x = min_corner.x;
    top_left.y = max_corner.y;
    envelope_msg.points.push_back(top_left);

    envelope_pub_->publish(envelope_msg);

    bg::model::multi_polygon<
        bg::model::polygon<bg::model::d2::point_xy<double>>>
        obstacle_polygons;

    bg::difference(envelope_box, concave_polygon, obstacle_polygons);

    safe_bayesian_optimization::msg::PolygonArray polygon_array_msg;
    for (const auto &polygon : obstacle_polygons) {
      geometry_msgs::msg::Polygon msg_polygon;
      for (const auto &point : polygon.outer()) {
        geometry_msgs::msg::Point32 p;
        p.x = point.get<0>();
        p.y = point.get<1>();
        msg_polygon.points.push_back(p);
      }
      polygon_array_msg.polygons.push_back(msg_polygon);
    }
    polygons_pub_->publish(polygon_array_msg);

    // Publish debug visualization image if enabled
    if (publish_debug_image_ && debug_image_pub_) {
      publish_debug_image(obstacle_polygons, alpha_shape_boundary);
    }
  }

  void
  publish_debug_image(const bg::model::multi_polygon<
                          bg::model::polygon<bg::model::d2::point_xy<double>>>
                          &obstacle_polygons,
                      const std::vector<std::array<double, 2>> &safe_hull) {
    // Calculate image bounds from data points
    if (D_.rows() == 0) {
      RCLCPP_WARN(this->get_logger(),
                  "No data points available for debug image");
      return;
    }

    double x_min = D_.col(0).minCoeff();
    double x_max = D_.col(0).maxCoeff();
    double y_min = D_.col(1).minCoeff();
    double y_max = D_.col(1).maxCoeff();

    // Image size - use higher resolution for better visualization
    int img_width = 800;
    int img_height = 600;

    // Create debug image
    cv::Mat debug_img = cv::Mat::zeros(img_height, img_width, CV_8UC3);

    // Helper function to convert world coordinates to image coordinates
    auto world_to_image = [&](double world_x, double world_y) -> cv::Point {
      int img_x =
          static_cast<int>((world_x - x_min) / (x_max - x_min) * img_width);
      int img_y = static_cast<int>((1.0 - (world_y - y_min) / (y_max - y_min)) *
                                   img_height);
      return cv::Point(img_x, img_y);
    };

    // Draw concave hull of safe set using pre-computed hull
    if (safe_hull.size() > 2) {
      std::vector<cv::Point> hull_img_points;
      hull_img_points.reserve(safe_hull.size()); // Reserve capacity

      for (const auto &point : safe_hull) {
        cv::Point img_point = world_to_image(point[0], point[1]);
        // Only add points within image bounds
        if (img_point.x >= 0 && img_point.x < img_width && img_point.y >= 0 &&
            img_point.y < img_height) {
          hull_img_points.push_back(img_point);
        }
      }

      if (hull_img_points.size() > 2) {
        const cv::Point *pts = hull_img_points.data();
        int npts = hull_img_points.size();
        cv::fillPoly(debug_img, &pts, &npts, 1, cv::Scalar(0, 150, 0));
        cv::polylines(debug_img, hull_img_points, true, cv::Scalar(0, 255, 0),
                      2);
      }
    }

    // Draw safe region as circles

    for (int i = 0; i < D_.rows(); ++i) {
      if (S_(i)) {
        cv::Point img_point = world_to_image(D_(i, 0), D_(i, 1));
        if (img_point.x >= 0 && img_point.x < img_width && img_point.y >= 0 &&
            img_point.y < img_height) {
          cv::circle(debug_img, img_point, 3, cv::Scalar(255, 0, 0), -1);
        }
      }
    }

    // Draw obstacle polygons in blue
    for (const auto &polygon : obstacle_polygons) {
      const auto &outer_ring = polygon.outer();
      if (outer_ring.size() < 3)
        continue; // Skip invalid polygons

      std::vector<cv::Point> cv_points;
      cv_points.reserve(outer_ring.size()); // Reserve capacity

      for (const auto &point : outer_ring) {
        cv::Point img_point = world_to_image(point.get<0>(), point.get<1>());
        cv_points.push_back(img_point);
      }

      if (cv_points.size() > 2) {
        const cv::Point *pts = cv_points.data();
        int npts = static_cast<int>(cv_points.size());
        cv::fillPoly(debug_img, &pts, &npts, 1, cv::Scalar(255, 100, 0));
        cv::polylines(debug_img, cv_points, true, cv::Scalar(255, 255, 0), 2);
      }
    }

    // // Draw current goal in magenta if available
    // if (current_goal_.point.x != 0.0 || current_goal_.point.y != 0.0) {
    //   cv::Point goal_img = world_to_image(current_goal_.point.x, current_goal_.point.y);
    //   if (goal_img.x >= 0 && goal_img.x < img_width &&
    //       goal_img.y >= 0 && goal_img.y < img_height) {
    //     cv::circle(debug_img, goal_img, 8, cv::Scalar(255, 0, 255), -1);
    //     cv::circle(debug_img, goal_img, 12, cv::Scalar(255, 255, 255), 2);
    //   }
    // }

    // Add legend
    cv::putText(debug_img, "Safe concave hull (green)", cv::Point(10, 30),
                cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 255, 0), 2);
    cv::putText(debug_img, "Obstacles (blue)", cv::Point(10, 55),
                cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(255, 100, 0), 2);
    cv::putText(debug_img, "Goal (magenta)", cv::Point(10, 80),
                cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(255, 0, 255), 2);

    // Convert OpenCV image to ROS image message
    std_msgs::msg::Header header;
    header.stamp = this->now();
    header.frame_id = "map";

    sensor_msgs::msg::Image::SharedPtr img_msg =
        cv_bridge::CvImage(header, "bgr8", debug_img).toImageMsg();

    debug_image_pub_->publish(*img_msg);

    RCLCPP_DEBUG(this->get_logger(), "Published debug visualization image");
  }

  void goal_point_callback(const geometry_msgs::msg::PointStamped::SharedPtr msg) {
    current_goal_.point = msg->point;
  }

  void publish_subgoal_marker(const geometry_msgs::msg::Point &subgoal) {
    visualization_msgs::msg::Marker marker;
    marker.header.frame_id = "map";
    marker.header.stamp = this->get_clock()->now();
    marker.ns = "subgoal";
    marker.id = 0;
    marker.type = visualization_msgs::msg::Marker::SPHERE;
    marker.action = visualization_msgs::msg::Marker::ADD;

    marker.pose.position.x = subgoal.x;
    marker.pose.position.y = subgoal.y;
    marker.pose.position.z = subgoal.z;
    marker.pose.orientation.w = 1.0;

    marker.scale.x = 0.3;
    marker.scale.y = 0.3;
    marker.scale.z = 0.3;

    marker.color.r = 0.0;
    marker.color.g = 1.0;
    marker.color.b = 0.0;
    marker.color.a = 1.0;

    subgoal_marker_pub_->publish(marker);
  }

  void publish_concave_markers(
      const bg::model::polygon<bg::model::d2::point_xy<double>>
          &concave_polygon) {
    auto marker_array = visualization_msgs::msg::MarkerArray();

    // Create a marker for the concave polygon
    auto marker = visualization_msgs::msg::Marker();
    marker.header.frame_id = "map";
    marker.header.stamp = this->now();
    marker.ns = "concave";
    marker.id = 0;
    marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
    marker.action = visualization_msgs::msg::Marker::ADD;

    // Set marker properties
    marker.scale.x = 0.06; // Line width
    marker.color.r = 0.0;
    marker.color.g = 0.0;
    marker.color.b = 1.0;
    marker.color.a = 0.8;

    // Convert concave polygon points to marker points
    for (const auto &boost_point : concave_polygon.outer()) {
      geometry_msgs::msg::Point point;
      point.x = boost_point.get<0>();
      point.y = boost_point.get<1>();
      point.z = 0.0;
      marker.points.push_back(point);
    }

    // Close the polygon by adding the first point at the end
    if (!marker.points.empty()) {
      marker.points.push_back(marker.points[0]);
    }

    marker_array.markers.push_back(marker);
    concave_markers_pub_->publish(marker_array);
  }

  void publish_projected_goal_marker(const geometry_msgs::msg::Point &original_goal,
                                   const geometry_msgs::msg::Point &projected_goal,
                                   double projection_distance) {
    visualization_msgs::msg::Marker line_marker;
    line_marker.header.frame_id = "map";
    line_marker.header.stamp = this->get_clock()->now();
    line_marker.ns = "projected_goal";
    line_marker.id = 0;
    line_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
    line_marker.action = visualization_msgs::msg::Marker::ADD;

    // Line from original goal to projected goal
    line_marker.points.push_back(original_goal);
    line_marker.points.push_back(projected_goal);

    line_marker.scale.x = 0.08; // Line width
    line_marker.color.r = 1.0;  // Red line
    line_marker.color.g = 0.0;
    line_marker.color.b = 0.0;
    line_marker.color.a = 0.9;

    projected_goal_marker_pub_->publish(line_marker);

    // Also publish arrow marker showing projection direction
    visualization_msgs::msg::Marker arrow_marker;
    arrow_marker.header.frame_id = "map";
    arrow_marker.header.stamp = this->get_clock()->now();
    arrow_marker.ns = "projected_goal";
    arrow_marker.id = 1;
    arrow_marker.type = visualization_msgs::msg::Marker::ARROW;
    arrow_marker.action = visualization_msgs::msg::Marker::ADD;

    arrow_marker.points.push_back(original_goal);
    arrow_marker.points.push_back(projected_goal);

    arrow_marker.scale.x = 0.05; // Shaft diameter
    arrow_marker.scale.y = 0.1;  // Head diameter
    arrow_marker.scale.z = 0.1;  // Head length
    arrow_marker.color.r = 0.8;
    arrow_marker.color.g = 0.2;
    arrow_marker.color.b = 0.0;
    arrow_marker.color.a = 0.8;

    // Small delay to ensure line is published first
    rclcpp::sleep_for(std::chrono::milliseconds(10));
    projected_goal_marker_pub_->publish(arrow_marker);
  }
};

int main(int argc, char *argv[]) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<OptimizerNode>();
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  executor.spin();
  rclcpp::shutdown();
  return 0;
}
