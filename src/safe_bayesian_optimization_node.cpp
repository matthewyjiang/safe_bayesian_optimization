#include "trusses_custom_interfaces/srv/get_terrain_map_with_uncertainty.hpp"
#include "trusses_custom_interfaces/srv/spatial_data.hpp"
#include <Eigen/Dense>
#include <geometry_msgs/msg/point.hpp>
#include <memory>
#include <opencv2/opencv.hpp>
#include <ratio>
#include <rclcpp/rclcpp.hpp>

#include <vector>

class OptimizerNode : public rclcpp::Node {
public:
  OptimizerNode() : Node("optimizer_node"), last_data_length_(0) {
    // Declare parameters
    this->declare_parameter("opt.beta", 2.0);
    this->declare_parameter("opt.f_min", 0.0);
    this->declare_parameter("terrain_map.width", 2.0);
    this->declare_parameter("terrain_map.height", 2.0);
    this->declare_parameter("terrain_map.width_cells", 20);
    this->declare_parameter("terrain_map.height_cells", 20);

    // Read parameters
    beta_ = this->get_parameter("opt.beta").as_double();
    f_min_ = this->get_parameter("opt.f_min").as_double();
    terrain_width_ = this->get_parameter("terrain_map.width").as_double();
    terrain_height_ = this->get_parameter("terrain_map.height").as_double();
    terrain_width_cells_ =
        this->get_parameter("terrain_map.width_cells").as_int();
    terrain_height_cells_ =
        this->get_parameter("terrain_map.height_cells").as_int();

    // Create service clients
    spatial_data_client_ =
        this->create_client<trusses_custom_interfaces::srv::SpatialData>(
            "get_spatial_data");
    terrain_map_client_ = this->create_client<
        trusses_custom_interfaces::srv::GetTerrainMapWithUncertainty>(
        "get_terrain_map_with_uncertainty");

    // Create publisher
    current_subgoal_pub_ = this->create_publisher<geometry_msgs::msg::Point>(
        "current_subgoal", 10);

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

  // Spatial data monitoring
  size_t last_data_length_;
  rclcpp::Client<trusses_custom_interfaces::srv::SpatialData>::SharedPtr
      spatial_data_client_;
  rclcpp::Client<trusses_custom_interfaces::srv::GetTerrainMapWithUncertainty>::
      SharedPtr terrain_map_client_;
  rclcpp::TimerBase::SharedPtr spatial_data_timer_;
  rclcpp::Publisher<geometry_msgs::msg::Point>::SharedPtr current_subgoal_pub_;

  // Parameters from config
  double terrain_width_;
  double terrain_height_;
  int terrain_width_cells_;
  int terrain_height_cells_;

  void ComputeSets() {
    // Compute the confidence intervals
    ComputeConfidenceIntervals();

    // Update the safe set based on the confidence intervals
    UpdateSafeSet();
  }

  void UpdateSafeSet() { S_ = Q_.col(0).array() > f_min_; }

  void ComputeConfidenceIntervals() {
    Eigen::VectorXd confidence = beta_ * std_;

    Q_.col(0) = mu_ - confidence;
    Q_.col(1) = mu_ + confidence;
  }

  std::vector<int> FindSafetyContourIndices() {
    if (D_.rows() == 0 || S_.rows() == 0) {
      RCLCPP_WARN(this->get_logger(), "D_ or S_ is empty");
      return {};
    }

    // Find grid bounds
    int min_x = static_cast<int>(D_.col(0).minCoeff());
    int max_x = static_cast<int>(D_.col(0).maxCoeff());
    int min_y = static_cast<int>(D_.col(1).minCoeff());
    int max_y = static_cast<int>(D_.col(1).maxCoeff());

    int width = max_x - min_x + 1;
    int height = max_y - min_y + 1;

    // Create binary image from safety data
    cv::Mat safety_image = cv::Mat::zeros(height, width, CV_8UC1);

    for (int i = 0; i < D_.rows(); ++i) {
      int x = static_cast<int>(D_(i, 0)) - min_x;
      int y = static_cast<int>(D_(i, 1)) - min_y;

      if (x >= 0 && x < width && y >= 0 && y < height) {
        safety_image.at<uchar>(y, x) = S_(i) ? 255 : 0;
      }
    }

    // Find contours
    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarchy;
    cv::findContours(safety_image, contours, hierarchy, cv::RETR_EXTERNAL,
                     cv::CHAIN_APPROX_NONE);

    // Convert contour points back to original indices
    std::vector<int> contour_indices;

    for (const auto &contour : contours) {
      for (const auto &point : contour) {
        int orig_x = point.x + min_x;
        int orig_y = point.y + min_y;

        // Find the index in D_ that matches this point
        for (int i = 0; i < D_.rows(); ++i) {
          if (static_cast<int>(D_(i, 0)) == orig_x &&
              static_cast<int>(D_(i, 1)) == orig_y) {
            contour_indices.push_back(i);
            break;
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
    std::vector<int> frontier_indices = FindSafetyContourIndices();

    if (frontier_indices.empty()) {
      RCLCPP_WARN(this->get_logger(), "No frontier points found");
      return -1;
    }

    // TODO: Implement selection algorithm on frontier_indices
    // For now, return the first frontier index as placeholder
    RCLCPP_INFO(this->get_logger(), "Found %lu frontier points for selection",
                frontier_indices.size());

    return frontier_indices[0];
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

            // Request terrain map after receiving spatial data
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
    // Use parameters from config
    request->width = terrain_width_;
    request->height = terrain_height_;
    request->n_width_cells = terrain_width_cells_;
    request->n_height_cells = terrain_height_cells_;

    RCLCPP_INFO(this->get_logger(),
                "Requesting terrain map with width: %f, height: %f, "
                "width_cells: %d, height_cells: %d",
                request->width, request->height, request->n_width_cells,
                request->n_height_cells);

    // Use callback-based async call
    terrain_map_client_->async_send_request(
        request,
        [this](rclcpp::Client<
               trusses_custom_interfaces::srv::GetTerrainMapWithUncertainty>::
                   SharedFuture future) {
          auto response = future.get();
          if (response->success) {
            RCLCPP_INFO(this->get_logger(),
                        "Received terrain map with %zu points",
                        response->values.size());

            // Process the terrain map data here
            // response->x_coords, response->y_coords, response->values,
            // response->uncertainties
            process_terrain_map(response);

          } else {
            RCLCPP_WARN(this->get_logger(), "Terrain map request failed: %s",
                        response->message.c_str());
          }
        });
  }

  void process_terrain_map(
      const std::shared_ptr<trusses_custom_interfaces::srv::
                                GetTerrainMapWithUncertainty::Response>
          response) {
    // Update the parameter set D_, mean mu_, and std_ with the new terrain map
    // data
    size_t num_points = response->x_coords.size();

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

    RCLCPP_INFO(this->get_logger(), "Updated terrain map data with %zu points",
                num_points);
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
