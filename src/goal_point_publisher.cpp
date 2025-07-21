#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/point.hpp>

class GoalPointPublisher : public rclcpp::Node {
public:
  GoalPointPublisher() : Node("goal_point_publisher") {
    publisher_ = this->create_publisher<geometry_msgs::msg::Point>("goal_point", 10);
    
    timer_ = this->create_wall_timer(
        std::chrono::seconds(1),
        std::bind(&GoalPointPublisher::publish_goal_point, this));
    
    RCLCPP_INFO(this->get_logger(), "Goal Point Publisher initialized");
  }

private:
  void publish_goal_point() {
    geometry_msgs::msg::Point msg;
    msg.x = 0.0;
    msg.y = 0.0;
    msg.z = 0.0;
    
    publisher_->publish(msg);
    RCLCPP_INFO(this->get_logger(), "Published goal point: (%.2f, %.2f)", msg.x, msg.y);
  }
  
  rclcpp::Publisher<geometry_msgs::msg::Point>::SharedPtr publisher_;
  rclcpp::TimerBase::SharedPtr timer_;
};

int main(int argc, char *argv[]) {
  rclcpp::init(argc, argv);
  auto node = std::make_shared<GoalPointPublisher>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}