#include <reactive_planner_lib.h>
#include <cmath>

DiffeoTransformResult computeDiffeoTransform(
    std::vector<double> &robot_position, double robot_orientation,
    std::vector<std::vector<PolygonClass>> &diffeo_tree_array,
    DiffeoParamsClass &diffeo_params, rclcpp::Logger logger) {

  std::cout << "\n=== DEBUG: DIFFEOMORPHISM TRANSFORMATION STARTED ===" << std::endl;
  std::cout << "Initial robot position: [" << robot_position[0] << ", " << robot_position[1] << "]" << std::endl;
  std::cout << "Initial robot orientation: " << robot_orientation << " rad (" << (robot_orientation * 180.0 / M_PI) << " deg)" << std::endl;
  std::cout << "Number of obstacle trees to process: " << diffeo_tree_array.size() << std::endl;

  std::vector<double> RobotPositionTransformed = {robot_position[0],
                                                  robot_position[1]};
  std::vector<std::vector<double>> RobotPositionTransformedD = {{1.0, 0.0},
                                                                {0.0, 1.0}};
  std::vector<double> RobotPositionTransformedDD = {0.0, 0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0, 0.0};

  std::cout << "\nInitial transformation state:" << std::endl;
  std::cout << "  Position: [" << RobotPositionTransformed[0] << ", " << RobotPositionTransformed[1] << "]" << std::endl;
  std::cout << "  Jacobian: [[" << RobotPositionTransformedD[0][0] << ", " << RobotPositionTransformedD[0][1] << "], [" 
            << RobotPositionTransformedD[1][0] << ", " << RobotPositionTransformedD[1][1] << "]]" << std::endl;

  for (size_t i = 0; i < diffeo_tree_array.size(); i++) {
    std::cout << "\n--- DEBUG: PROCESSING OBSTACLE TREE " << i << " ---" << std::endl;
    std::cout << "  Tree " << i << " contains " << diffeo_tree_array[i].size() << " convex pieces" << std::endl;
    
    RCLCPP_INFO(logger, "Transform iteration %zu: position before [%.3f, %.3f]",
                i, RobotPositionTransformed[0], RobotPositionTransformed[1]);

    std::cout << "  Input to polygonDiffeoConvex:" << std::endl;
    std::cout << "    Position: [" << RobotPositionTransformed[0] << ", " << RobotPositionTransformed[1] << "]" << std::endl;

    OutputStructVector TempTransformation = polygonDiffeoConvex(
        RobotPositionTransformed, diffeo_tree_array[i], diffeo_params);

    std::vector<double> TempPositionTransformed = TempTransformation.Value;
    std::vector<std::vector<double>> TempPositionTransformedD =
        TempTransformation.Jacobian;
    std::vector<double> TempPositionTransformedDD =
        TempTransformation.JacobianD;

    std::cout << "  Output from polygonDiffeoConvex:" << std::endl;
    std::cout << "    Transformed position: [" << TempPositionTransformed[0] << ", " << TempPositionTransformed[1] << "]" << std::endl;
    std::cout << "    Jacobian: [[" << TempPositionTransformedD[0][0] << ", " << TempPositionTransformedD[0][1] << "], [" 
              << TempPositionTransformedD[1][0] << ", " << TempPositionTransformedD[1][1] << "]]" << std::endl;
    std::cout << "    Hessian: [" << TempPositionTransformedDD[0] << ", " << TempPositionTransformedDD[1] << ", " 
              << TempPositionTransformedDD[2] << ", " << TempPositionTransformedDD[3] << ", " 
              << TempPositionTransformedDD[4] << ", " << TempPositionTransformedDD[5] << ", " 
              << TempPositionTransformedDD[6] << ", " << TempPositionTransformedDD[7] << "]" << std::endl;

    double res1 =
        TempPositionTransformedD[0][0] * RobotPositionTransformedDD[0] +
        TempPositionTransformedD[0][1] * RobotPositionTransformedDD[4] +
        RobotPositionTransformedD[0][0] *
            (TempPositionTransformedDD[0] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[1] * RobotPositionTransformedD[1][0]) +
        RobotPositionTransformedD[1][0] *
            (TempPositionTransformedDD[2] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[3] * RobotPositionTransformedD[1][0]);
    double res2 =
        TempPositionTransformedD[0][0] * RobotPositionTransformedDD[1] +
        TempPositionTransformedD[0][1] * RobotPositionTransformedDD[5] +
        RobotPositionTransformedD[0][0] *
            (TempPositionTransformedDD[0] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[1] * RobotPositionTransformedD[1][1]) +
        RobotPositionTransformedD[1][0] *
            (TempPositionTransformedDD[2] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[3] * RobotPositionTransformedD[1][1]);
    double res3 =
        TempPositionTransformedD[0][0] * RobotPositionTransformedDD[2] +
        TempPositionTransformedD[0][1] * RobotPositionTransformedDD[6] +
        RobotPositionTransformedD[0][1] *
            (TempPositionTransformedDD[0] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[1] * RobotPositionTransformedD[1][0]) +
        RobotPositionTransformedD[1][1] *
            (TempPositionTransformedDD[2] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[3] * RobotPositionTransformedD[1][0]);
    double res4 =
        TempPositionTransformedD[0][0] * RobotPositionTransformedDD[3] +
        TempPositionTransformedD[0][1] * RobotPositionTransformedDD[7] +
        RobotPositionTransformedD[0][1] *
            (TempPositionTransformedDD[0] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[1] * RobotPositionTransformedD[1][1]) +
        RobotPositionTransformedD[1][1] *
            (TempPositionTransformedDD[2] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[3] * RobotPositionTransformedD[1][1]);
    double res5 =
        TempPositionTransformedD[1][0] * RobotPositionTransformedDD[0] +
        TempPositionTransformedD[1][1] * RobotPositionTransformedDD[4] +
        RobotPositionTransformedD[0][0] *
            (TempPositionTransformedDD[4] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[5] * RobotPositionTransformedD[1][0]) +
        RobotPositionTransformedD[1][0] *
            (TempPositionTransformedDD[6] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[7] * RobotPositionTransformedD[1][0]);
    double res6 =
        TempPositionTransformedD[1][0] * RobotPositionTransformedDD[1] +
        TempPositionTransformedD[1][1] * RobotPositionTransformedDD[5] +
        RobotPositionTransformedD[0][0] *
            (TempPositionTransformedDD[4] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[5] * RobotPositionTransformedD[1][1]) +
        RobotPositionTransformedD[1][0] *
            (TempPositionTransformedDD[6] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[7] * RobotPositionTransformedD[1][1]);
    double res7 =
        TempPositionTransformedD[1][0] * RobotPositionTransformedDD[2] +
        TempPositionTransformedD[1][1] * RobotPositionTransformedDD[6] +
        RobotPositionTransformedD[0][1] *
            (TempPositionTransformedDD[4] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[5] * RobotPositionTransformedD[1][0]) +
        RobotPositionTransformedD[1][1] *
            (TempPositionTransformedDD[6] * RobotPositionTransformedD[0][0] +
             TempPositionTransformedDD[7] * RobotPositionTransformedD[1][0]);
    double res8 =
        TempPositionTransformedD[1][0] * RobotPositionTransformedDD[3] +
        TempPositionTransformedD[1][1] * RobotPositionTransformedDD[7] +
        RobotPositionTransformedD[0][1] *
            (TempPositionTransformedDD[4] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[5] * RobotPositionTransformedD[1][1]) +
        RobotPositionTransformedD[1][1] *
            (TempPositionTransformedDD[6] * RobotPositionTransformedD[0][1] +
             TempPositionTransformedDD[7] * RobotPositionTransformedD[1][1]);

    std::cout << "  Chain rule hessian computation results:" << std::endl;
    std::cout << "    res1-8: [" << res1 << ", " << res2 << ", " << res3 << ", " << res4 << ", " 
              << res5 << ", " << res6 << ", " << res7 << ", " << res8 << "]" << std::endl;

    RobotPositionTransformedDD[0] = res1;
    RobotPositionTransformedDD[1] = res2;
    RobotPositionTransformedDD[2] = res3;
    RobotPositionTransformedDD[3] = res4;
    RobotPositionTransformedDD[4] = res5;
    RobotPositionTransformedDD[5] = res6;
    RobotPositionTransformedDD[6] = res7;
    RobotPositionTransformedDD[7] = res8;

    std::cout << "  Before chain rule composition:" << std::endl;
    std::cout << "    Old Jacobian: [[" << RobotPositionTransformedD[0][0] << ", " << RobotPositionTransformedD[0][1] << "], [" 
              << RobotPositionTransformedD[1][0] << ", " << RobotPositionTransformedD[1][1] << "]]" << std::endl;

    RobotPositionTransformedD = MatrixMatrixMultiplication(
        TempPositionTransformedD, RobotPositionTransformedD);
    RobotPositionTransformed = TempPositionTransformed;

    std::cout << "  After chain rule composition:" << std::endl;
    std::cout << "    New Jacobian: [[" << RobotPositionTransformedD[0][0] << ", " << RobotPositionTransformedD[0][1] << "], [" 
              << RobotPositionTransformedD[1][0] << ", " << RobotPositionTransformedD[1][1] << "]]" << std::endl;
    std::cout << "    Final position: [" << RobotPositionTransformed[0] << ", " << RobotPositionTransformed[1] << "]" << std::endl;

    RCLCPP_INFO(logger, "Transform iteration %zu: position after [%.3f, %.3f]",
                i, RobotPositionTransformed[0], RobotPositionTransformed[1]);
  }

  std::cout << "\n=== DEBUG: FINAL ORIENTATION CALCULATIONS ===" << std::endl;
  std::cout << "Original robot orientation: " << robot_orientation << " rad (" << (robot_orientation * 180.0 / M_PI) << " deg)" << std::endl;
  std::cout << "cos(orientation): " << cos(robot_orientation) << ", sin(orientation): " << sin(robot_orientation) << std::endl;

  // Find alpha1, alpha2, beta1, beta2
  double alpha1 = -(RobotPositionTransformedD[1][0] * cos(robot_orientation) +
                    RobotPositionTransformedD[1][1] * sin(robot_orientation));
  double beta1 =
      RobotPositionTransformedDD[0] * pow(cos(robot_orientation), 2) +
      (RobotPositionTransformedDD[1] + RobotPositionTransformedDD[2]) *
          sin(robot_orientation) * cos(robot_orientation) +
      RobotPositionTransformedDD[3] * pow(sin(robot_orientation), 2);
  double alpha2 = RobotPositionTransformedD[0][0] * cos(robot_orientation) +
                  RobotPositionTransformedD[0][1] * sin(robot_orientation);
  double beta2 =
      RobotPositionTransformedDD[4] * pow(cos(robot_orientation), 2) +
      (RobotPositionTransformedDD[5] + RobotPositionTransformedDD[6]) *
          sin(robot_orientation) * cos(robot_orientation) +
      RobotPositionTransformedDD[7] * pow(sin(robot_orientation), 2);

  std::cout << "Computed parameters:" << std::endl;
  std::cout << "  alpha1: " << alpha1 << std::endl;
  std::cout << "  alpha2: " << alpha2 << std::endl;
  std::cout << "  beta1: " << beta1 << std::endl;
  std::cout << "  beta2: " << beta2 << std::endl;

  // Find transformed orientation
  double y_component = RobotPositionTransformedD[1][0] * cos(robot_orientation) +
                       RobotPositionTransformedD[1][1] * sin(robot_orientation);
  double x_component = RobotPositionTransformedD[0][0] * cos(robot_orientation) +
                       RobotPositionTransformedD[0][1] * sin(robot_orientation);
  
  std::cout << "Orientation transformation:" << std::endl;
  std::cout << "  y_component: " << y_component << std::endl;
  std::cout << "  x_component: " << x_component << std::endl;
  
  double RobotOrientationTransformed = atan2(y_component, x_component);
  
  std::cout << "  Transformed orientation: " << RobotOrientationTransformed << " rad (" << (RobotOrientationTransformed * 180.0 / M_PI) << " deg)" << std::endl;

  std::cout << "\n=== DEBUG: DIFFEOMORPHISM TRANSFORMATION COMPLETED ===" << std::endl;
  std::cout << "TRANSFORMATION SUMMARY:" << std::endl;
  std::cout << "  Original position: [" << robot_position[0] << ", " << robot_position[1] << "]" << std::endl;
  std::cout << "  Final position:    [" << RobotPositionTransformed[0] << ", " << RobotPositionTransformed[1] << "]" << std::endl;
  std::cout << "  Position delta:    [" << (RobotPositionTransformed[0] - robot_position[0]) << ", " << (RobotPositionTransformed[1] - robot_position[1]) << "]" << std::endl;
  std::cout << "  Original orientation: " << (robot_orientation * 180.0 / M_PI) << " deg" << std::endl;
  std::cout << "  Final orientation:    " << (RobotOrientationTransformed * 180.0 / M_PI) << " deg" << std::endl;
  std::cout << "  Orientation delta:    " << ((RobotOrientationTransformed - robot_orientation) * 180.0 / M_PI) << " deg" << std::endl;

  DiffeoTransformResult result;
  result.transformed_position = RobotPositionTransformed;
  result.transformed_jacobian = RobotPositionTransformedD;
  result.transformed_hessian = RobotPositionTransformedDD;
  result.alpha1 = alpha1;
  result.alpha2 = alpha2;
  result.beta1 = beta1;
  result.beta2 = beta2;
  result.transformed_orientation = RobotOrientationTransformed;

  std::cout << "Result structure populated successfully!" << std::endl;
  return result;
}
