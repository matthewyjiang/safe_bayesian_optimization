#include <geometry_msgs/msg/point_stamped.hpp>
#include <rclcpp/rclcpp.hpp>
#include <visualization_msgs/msg/marker.hpp>

class GoalPointPublisher : public rclcpp::Node {
public:
  GoalPointPublisher() : Node("goal_point_publisher") {
    publisher_ =
        this->create_publisher<geometry_msgs::msg::PointStamped>("goal_point", 10);

    marker_pub_ = this->create_publisher<visualization_msgs::msg::Marker>(
        "goal_marker", 10);

    timer_ = this->create_wall_timer(
        std::chrono::seconds(1),
        std::bind(&GoalPointPublisher::publish_goal_point, this));

    RCLCPP_INFO(this->get_logger(), "Goal Point Publisher initialized");
  }

private:
  void publish_goal_point() {
    geometry_msgs::msg::PointStamped msg;
    msg.header.frame_id = "map";
    msg.header.stamp = this->get_clock()->now();
    msg.point.x = 8.0;
    msg.point.y = 4.0;
    msg.point.z = 0.0;

    publisher_->publish(msg);

    // Publish goal marker for visualization
    publish_goal_marker(msg.point);

    RCLCPP_INFO(this->get_logger(), "Published goal point: (%.2f, %.2f)", msg.point.x,
                msg.point.y);
  }

  void publish_goal_marker(const geometry_msgs::msg::Point &goal) {
    visualization_msgs::msg::Marker marker;
    marker.header.frame_id = "map";
    marker.header.stamp = this->get_clock()->now();
    marker.ns = "goal";
    marker.id = 0;
    marker.type = visualization_msgs::msg::Marker::SPHERE;
    marker.action = visualization_msgs::msg::Marker::ADD;

    marker.pose.position.x = goal.x;
    marker.pose.position.y = goal.y;
    marker.pose.position.z = goal.z;
    marker.pose.orientation.w = 1.0;

    marker.scale.x = 0.4;
    marker.scale.y = 0.4;
    marker.scale.z = 0.4;

    marker.color.r = 1.0;
    marker.color.g = 0.0;
    marker.color.b = 0.0;
    marker.color.a = 1.0;

    marker_pub_->publish(marker);
  }

  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr publisher_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr marker_pub_;
  rclcpp::TimerBase::SharedPtr timer_;
};

int main(int argc, char *argv[]) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<GoalPointPublisher>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
