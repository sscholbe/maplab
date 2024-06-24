#include "map-optimization/vi-optimization-builder.h"

#include <gflags/gflags.h>
#include <vi-map-helpers/mission-clustering-coobservation.h>

#include "map-optimization/augment-loopclosure.h"
#include "map-optimization/optimization-state-fixing.h"

DEFINE_bool(
    ba_include_visual, true, "Whether or not to include visual error-terms.");
DEFINE_bool(
    ba_use_visual_outlier_rejection_solver, true,
    "Reject outlier landmarks during the solve?");
DEFINE_string(
    ba_feature_type, "Binary",
    "Type of features to use in bundle adjustment, will default to using all "
    "of them at once.");
DEFINE_bool(
    ba_include_inertial, true, "Whether or not to include IMU error-terms.");
DEFINE_bool(
    ba_include_wheel_odometry, false,
    "Whether or not to include wheel (relative pose) error-terms.");

DEFINE_bool(
    ba_fix_ncamera_intrinsics, true,
    "Whether or not to fix the intrinsics of the ncamera(s).");
DEFINE_bool(
    ba_fix_ncamera_extrinsics_rotation, true,
    "Whether or not to fix the rotation extrinsics of the ncamera(s).");
DEFINE_bool(
    ba_fix_ncamera_extrinsics_translation, true,
    "Whether or not to fix the translation extrinsics of the ncamera(s).");
DEFINE_bool(
    ba_fix_wheel_odometry_extrinsics, true,
    "Whether or not to fix the extrinsics of the wheel odometry sensor.");
DEFINE_bool(
    ba_fix_vertices, false,
    "Whether or not to fix vertex poses in optimization.");
DEFINE_bool(
    ba_fix_landmark_positions, false,
    "Whether or not to fix the positions of the landmarks.");
DEFINE_bool(
    ba_fix_accel_bias, false,
    "Whether or not to fix the bias of the IMU accelerometer.");
DEFINE_bool(
    ba_fix_gyro_bias, false,
    "Whether or not to fix the bias of the IMU gyroscope.");
DEFINE_bool(
    ba_fix_velocity, false,
    "Whether or not to fix the velocity of the vertices.");
DEFINE_double(
    ba_latitude, common::locations::kLatitudeZurichDegrees,
    "Latitude to estimate the gravity magnitude.");
DEFINE_double(
    ba_altitude_meters, common::locations::kAltitudeZurichMeters,
    "Altitude in meters to estimate the gravity magnitude.");

DEFINE_int32(
    ba_min_landmarks_per_frame, 0,
    "Minimum number of landmarks a frame must observe to be included in the "
    "problem.");

DEFINE_bool(
    ba_include_6dof_odometry, false,
    "Whether or not to include 6DoF odometry constraints (T_St_Stp1) for an "
    "odometry sensor frame S.");

DEFINE_bool(
    ba_fix_6dof_odometry_extrinsics, true,
    "Whether or not to fix the extrinsics of the 6DoF odometry sensor.");

DEFINE_bool(
    ba_include_absolute_pose_constraints, false,
    "Whether or not to include absolute 6DoF pose constraints (T_G_B), e.g. "
    "GPS or AprilTag detections.");

DEFINE_bool(
    ba_fix_absolute_pose_sensor_extrinsics, true,
    "Whether or not to fix the extrinsics of the absolute 6DoF pose "
    "constraints sensor (T_B_S), e.g. GPS to IMU calibration. This flag will "
    "only take effect if absolute pose constraints are enabled!");

DEFINE_bool(
    ba_absolute_pose_sensor_fix_mission_baseframes, true,
    "Whether or not to fix the mission baseframes (T_G_M) during "
    "optimization. This flag will only take effect if absolute pose "
    "constraints are enabled!");

DEFINE_bool(
    ba_include_loop_closure_edges, false,
    "Whether or not to add the loop closure edges present in the map to the "
    "optimization problem. The visual loop closure does not use these edges by "
    "default, but merges the landmarks (i.e. loop closures are part of the "
    "visual error terms). Pose graph relaxation on the other hand adds visual "
    "loop closures as edges by default.");

DEFINE_double(ba_appd_weight, 0.0, "The weight of the APPD drift");
DEFINE_double(ba_huber_delta, 0.25, "The weight of the VCReg");
DEFINE_string(
    ba_constant_intrinsics_drift, "",
    "Constant intrinsics (for drift cameras)");
DEFINE_string(
    ba_constant_intrinsics_base, "", "Constant intrinsics (for base camera)");

namespace map_optimization {
ViProblemOptions ViProblemOptions::initFromGFlags() {
  ViProblemOptions options;

  // Vertex constraints
  options.fix_vertices = FLAGS_ba_fix_vertices;

  // Inertial constraints
  options.add_inertial_constraints = FLAGS_ba_include_inertial;
  options.fix_gyro_bias = FLAGS_ba_fix_gyro_bias;
  options.fix_accel_bias = FLAGS_ba_fix_accel_bias;
  options.fix_velocity = FLAGS_ba_fix_velocity;

  common::GravityProvider gravity_provider(
      FLAGS_ba_altitude_meters, FLAGS_ba_latitude);
  options.gravity_magnitude = gravity_provider.getGravityMagnitude();

  // Visual constraints.
  options.add_visual_constraints = FLAGS_ba_include_visual;
  options.fix_intrinsics = FLAGS_ba_fix_ncamera_intrinsics;
  options.fix_extrinsics_rotation = FLAGS_ba_fix_ncamera_extrinsics_rotation;
  options.fix_extrinsics_translation =
      FLAGS_ba_fix_ncamera_extrinsics_translation;
  options.fix_landmark_positions = FLAGS_ba_fix_landmark_positions;
  options.min_landmarks_per_frame = FLAGS_ba_min_landmarks_per_frame;
  options.feature_type = vi_map::StringToFeatureType(FLAGS_ba_feature_type);

  // Wheel odometry constraints
  options.add_wheel_odometry_constraints = FLAGS_ba_include_wheel_odometry;
  options.fix_wheel_extrinsics = FLAGS_ba_fix_wheel_odometry_extrinsics;

  // 6DoF odometry constraints.
  options.add_6dof_odometry_constraints = FLAGS_ba_include_6dof_odometry;
  options.fix_6dof_odometry_extrinsics = FLAGS_ba_fix_6dof_odometry_extrinsics;

  // Absolute 6DoF pose constraints
  options.add_absolute_pose_constraints =
      FLAGS_ba_include_absolute_pose_constraints;
  options.fix_absolute_pose_sensor_extrinsics =
      FLAGS_ba_fix_absolute_pose_sensor_extrinsics;
  options.fix_baseframes = FLAGS_ba_absolute_pose_sensor_fix_mission_baseframes;

  // Loop closure constraints (can be from an external source)
  options.add_loop_closure_edges = FLAGS_ba_include_loop_closure_edges;

  options.solver_options = initSolverOptionsFromFlags();

  options.enable_visual_outlier_rejection =
      FLAGS_ba_use_visual_outlier_rejection_solver;
  options.visual_outlier_rejection_options =
      map_optimization::OutlierRejectionSolverOptions::initFromFlags();

  options.printToConsole();

  return options;
}

#define NUM_GRID 6

struct APPDCostFunctor3 {
  APPDCostFunctor3(double weight, unsigned width, unsigned height)
      : weight_{weight}, width_{width}, height_{height} {}

  template <typename T>
  inline void unproject(
      const T* const I0, const T* const I1, const T* const I2, const T& x,
      const T& y, T& x_out, T& y_out, T& z_out) const {
    // intrinsics are delta to C0
    T avg_fx = I0[0] + (I1[0] + I2[0]) / 2.0;
    T avg_fy = I0[1] + (I1[1] + I2[1]) / 2.0;
    T avg_cx = I0[2] + (I1[2] + I2[2]) / 2.0;
    T avg_cy = I0[3] + (I1[3] + I2[3]) / 2.0;
    x_out = (x - avg_cx) / avg_fx;
    y_out = (y - avg_cy) / avg_fy;
    z_out = T(1.0);
    T div = T(
        1.0);  // / ceres::sqrt(x_out * x_out + y_out * y_out + z_out * z_out);
    x_out *= div;
    y_out *= div;
    z_out *= div;
  }

  template <typename T>
  inline void project(
      const T* const I0, const T* const D0, const T* const I, const T* const D,
      const T& x, const T& y, const T& z, T& x_out, T& y_out,
      T& scale_out) const {
    // intrinsics and distortion are delta to C0
    T fx = I0[0] + I[0];
    T fy = I0[1] + I[1];
    T cx = I0[2] + I[2];
    T cy = I0[3] + I[3];
    T k1 = D0[0] + D[0];
    T k2 = D0[1] + D[1];
    T k3 = D0[2] + D[2];
    T k4 = D0[3] + D[3];

    T x_z = x;  // / z;
    T y_z = y;  // / z;
    T x2 = x_z * x_z;
    T y2 = y_z * y_z;
    T r = ceres::sqrt(x2 + y2);

    T theta = ceres::atan(r);
    T theta2 = theta * theta;
    T theta4 = theta2 * theta2;
    T theta6 = theta2 * theta4;
    T theta8 = theta4 * theta4;
    T thetad =
        theta * (1.0 + k1 * theta2 + k2 * theta4 + k3 * theta6 + k4 * theta8);
    T scaling = thetad / r;

    scale_out = scaling;
    x_out = x_z * fx + cx;
    y_out = y_z * fy + cy;
  }

  template <typename T>
  bool operator()(
      const T* const I0, const T* const D0, const T* const I1,
      const T* const D1, const T* const I2, const T* const D2,
      T* residual) const {
    int rc = 0;
    for (int ix = 0; ix <= NUM_GRID - 1; ix++) {
      for (int iy = 0; iy <= NUM_GRID - 1; iy++) {
        T ux, uy, uz;
        unproject(
            I0, I1, I2, T(width_ * ix / (NUM_GRID - 1)),
            T(height_ * iy / (NUM_GRID - 1)), ux, uy, uz);
        T dx1, dy1, sc1;
        project(I0, D0, I1, D1, ux, uy, uz, dx1, dy1, sc1);
        T dx2, dy2, sc2;
        project(I0, D0, I2, D2, ux, uy, uz, dx2, dy2, sc2);
        residual[rc] = weight_ * (dx1 - dx2) / T(NUM_GRID * NUM_GRID);
        residual[rc + 1] = weight_ * (dy1 - dy2) / T(NUM_GRID * NUM_GRID);
        residual[rc + 2] =
            weight_ * T(450.0) * (sc1 - sc2) / T(NUM_GRID * NUM_GRID);
        residual[rc + 3] =
            weight_ * T(450.0) * (sc1 - sc2) / T(NUM_GRID * NUM_GRID);
        rc += 4;
      }
    }

    return true;
  }

 private:
  const double weight_;
  const double width_;
  const double height_;
};

std::vector<std::string> splitStringByComma(const std::string& str) {
  std::vector<std::string> result;
  std::stringstream ss(str);
  std::string item;

  while (std::getline(ss, item, ',')) {
    result.push_back(item);
  }

  return result;
}

OptimizationProblem* constructOptimizationProblem(
    const vi_map::MissionIdSet& mission_ids, const ViProblemOptions& options,
    vi_map::VIMap* map) {
  CHECK(map);
  CHECK(options.isValid());

  OptimizationProblem* problem = new OptimizationProblem(map, mission_ids);
  size_t num_visual_constraints_added = 0u;
  if (options.add_visual_constraints) {
    num_visual_constraints_added = addLandmarkTerms(
        options.feature_type, options.fix_landmark_positions,
        options.fix_intrinsics, options.fix_extrinsics_rotation,
        options.fix_extrinsics_translation, options.min_landmarks_per_frame,
        problem);
    if (num_visual_constraints_added == 0u) {
      LOG(WARNING)
          << "WARNING: Visual constraints enabled, but none "
          << "were found, adapting DoF settings of optimization problem...";
    }

    ceres_error_terms::ProblemInformation* problem_information =
        CHECK_NOTNULL(problem->getProblemInformationMutable());

    //
    // Add the drift terms
    //

    std::shared_ptr<ceres::LossFunction> loss_function(nullptr); /*
        new ceres::HuberLoss(FLAGS_ba_huber_delta)*/

    vi_map::SensorManager& sensor_manager = map->getSensorManager();

    static double* delta_0_intrinsics = new double[24];
    static double* delta_0_distortion = new double[24];
    problem_information->setParameterBlockConstant(delta_0_intrinsics);
    problem_information->setParameterBlockConstant(delta_0_distortion);

    for (const vi_map::MissionId& mission_id : mission_ids) {
      vi_map::VIMission& mission = map->getMission(mission_id);
      if (mission.drift_ncamera_ids.empty()) {
        continue;
      }

      LOG(INFO) << "Using camera drift for mission " << mission_id.shortHex()
                << " with weight " << FLAGS_ba_appd_weight;

      aslam::NCamera::Ptr ncam_r = map->getMissionNCameraPtr(mission_id);
      aslam::NCamera::Ptr ncam_a = ncam_r;

      // Might bug if different resolutions
      std::shared_ptr<ceres::CostFunction> cost_function2(
          new ceres::AutoDiffCostFunction<
              APPDCostFunctor3, NUM_GRID * NUM_GRID * 4, 4, 4, 4, 4, 4, 4>(
              new APPDCostFunctor3(
                  FLAGS_ba_appd_weight, ncam_r->getCamera(0).imageWidth(),
                  ncam_r->getCamera(0).imageHeight())));

      std::vector<std::string> fixs_base =
          splitStringByComma(FLAGS_ba_constant_intrinsics_base);
      std::vector<std::string> fixs_drift =
          splitStringByComma(FLAGS_ba_constant_intrinsics_drift);

      std::vector<int> fix_base_intrinsics;
      std::vector<int> fix_base_distortion;
      std::vector<int> fix_drift_intrinsics;
      std::vector<int> fix_drift_distortion;
      for (std::string& f : fixs_base) {
        if (f == "fx") {
          fix_base_intrinsics.push_back(0);
        } else if (f == "fy") {
          fix_base_intrinsics.push_back(1);
        } else if (f == "cx") {
          fix_base_intrinsics.push_back(2);
        } else if (f == "cy") {
          fix_base_intrinsics.push_back(3);
        } else if (f == "k1") {
          fix_base_distortion.push_back(0);
        } else if (f == "k2") {
          fix_base_distortion.push_back(1);
        } else if (f == "k3") {
          fix_base_distortion.push_back(2);
        } else if (f == "k4") {
          fix_base_distortion.push_back(3);
        }
      }
      for (std::string& f : fixs_drift) {
        if (f == "fx") {
          fix_drift_intrinsics.push_back(0);
        } else if (f == "fy") {
          fix_drift_intrinsics.push_back(1);
        } else if (f == "cx") {
          fix_drift_intrinsics.push_back(2);
        } else if (f == "cy") {
          fix_drift_intrinsics.push_back(3);
        } else if (f == "k1") {
          fix_drift_distortion.push_back(0);
        } else if (f == "k2") {
          fix_drift_distortion.push_back(1);
        } else if (f == "k3") {
          fix_drift_distortion.push_back(2);
        } else if (f == "k4") {
          fix_drift_distortion.push_back(3);
        }
      }

      std::shared_ptr<ceres::LocalParameterization> par_base_intrinsics(
          new ceres::SubsetParameterization(4, fix_base_intrinsics));
      std::shared_ptr<ceres::LocalParameterization> par_base_distortion(
          new ceres::SubsetParameterization(4, fix_base_distortion));
      std::shared_ptr<ceres::LocalParameterization> par_drift_intrinsics(
          new ceres::SubsetParameterization(4, fix_drift_intrinsics));
      std::shared_ptr<ceres::LocalParameterization> par_drift_distortion(
          new ceres::SubsetParameterization(4, fix_drift_distortion));

      /*std::shared_ptr<ceres::LocalParameterization> fix_principal(
          new SingleVariableConstantParameterization());*/

      /*pose_graph::VertexIdList vertex_ids;
      map->getAllVertexIdsInMissionAlongGraph(mission_id, &vertex_ids);

      size_t total_keypoints = 0;
      size_t num_cameras = 0;

      std::vector<size_t> keypoints;

      for (size_t ncam_idx = 0; ncam_idx < mission.drift_ncamera_ids.size();
           ncam_idx++) {
        auto& vertex = map->getVertex(vertex_ids.at(ncam_idx + 1));
        aslam::NCamera::Ptr ncam = sensor_manager.getSensorPtr<aslam::NCamera>(
            mission.drift_ncamera_ids.at(ncam_idx));

        for (size_t cam_idx = 0; cam_idx < ncam->getNumCameras(); cam_idx++) {
          const size_t num_keypoints =
              vertex.getVisualFrame(cam_idx).getNumKeypointMeasurements();
          total_keypoints += num_keypoints;
          keypoints.push_back(num_keypoints);
          num_cameras++;
        }
      }

      double avg_keypoints = total_keypoints / double(num_cameras);
      std::sort(keypoints.begin(), keypoints.end());

      size_t median_keypoints = keypoints.at(keypoints.size() / 2);
      LOG(INFO) << "Median keypoints: " << median_keypoints
                << ", average keypoints: " << avg_keypoints;*/

      for (size_t ncam_idx = 0; ncam_idx < mission.drift_ncamera_ids.size();
           ncam_idx++) {
        aslam::NCamera::Ptr ncam_b =
            sensor_manager.getSensorPtr<aslam::NCamera>(
                mission.drift_ncamera_ids.at(ncam_idx));
        // auto& vertex = map->getVertex(vertex_ids.at(ncam_idx + 1));

        // For each individual camera, add the drift terms
        for (size_t cam_idx = 0; cam_idx < ncam_b->getNumCameras(); cam_idx++) {
          aslam::Camera& cam_r = ncam_r->getCameraMutable(cam_idx);
          aslam::Camera& cam_a = ncam_a->getCameraMutable(cam_idx);
          aslam::Camera& cam_b = ncam_b->getCameraMutable(cam_idx);

          /*const size_t num_keypoints =
              vertex.getVisualFrame(cam_idx).getNumKeypointMeasurements();*/

          // Make distortion constant
          /*problem_information->setParameterBlockConstantIfPartOfTheProblem(
              cam_a.getDistortionMutable()->getParametersMutable());
          problem_information->setParameterBlockConstantIfPartOfTheProblem(
              cam_b.getDistortionMutable()->getParametersMutable());*/

          // Make principal point constant
          /*problem_information->setParameterization(
              cam_b.getParametersMutable(), fix_principal);*/

          /* Make k2-k4 constant */
          /*problem_information->setParameterization(
              cam_a.getDistortionMutable()->getParametersMutable(), fix_k2k3k4);
          problem_information->setParameterization(
              cam_b.getDistortionMutable()->getParametersMutable(),
          fix_k2k3k4);*/

          // Make intrinsics constant
          /*problem_information->setParameterBlockConstantIfPartOfTheProblem(
              cam_a.getParametersMutable());
          problem_information->setParameterBlockConstantIfPartOfTheProblem(
              cam_b.getParametersMutable());*/

          bool a_is_ref = ncam_a->getId() == ncam_r->getId();

          if (a_is_ref) {
            problem_information->setParameterization(
                cam_a.getDistortionMutable()->getParametersMutable(),
                par_base_distortion);
            problem_information->setParameterization(
                cam_a.getParametersMutable(), par_base_intrinsics);
          }

          problem_information->setParameterization(
              cam_b.getDistortionMutable()->getParametersMutable(),
              par_drift_distortion);
          problem_information->setParameterization(
              cam_b.getParametersMutable(), par_drift_intrinsics);

          /*std::shared_ptr<ceres::CostFunction> cost_function(
              new ceres::AutoDiffCostFunction<
                  PinholeAPPDCostFunctor, 8, 4, 4, 4>(
                  new PinholeAPPDCostFunctor(
                      FLAGS_ba_appd_weight, cam_b.imageWidth(),
                      cam_b.imageHeight())));
          problem_information->addResidualBlock(
              ceres_error_terms::ResidualType::kGenericPrior, cost_function,
              loss_function,
              {cam_r.getParametersMutable(),
               a_is_ref ? delta_0_intrinsics : cam_a.getParametersMutable(),
               cam_b.getParametersMutable()});*/

          problem_information->addResidualBlock(
              ceres_error_terms::ResidualType::kGenericPrior, cost_function2,
              loss_function,
              {cam_r.getParametersMutable(),
               cam_r.getDistortionMutable()->getParametersMutable(),
               a_is_ref ? delta_0_intrinsics : cam_a.getParametersMutable(),
               a_is_ref ? delta_0_distortion
                        : cam_a.getDistortionMutable()->getParametersMutable(),
               cam_b.getParametersMutable(),
               cam_b.getDistortionMutable()->getParametersMutable()});

          // VCREG
          /*std::shared_ptr<ceres::CostFunction> cost_function3(
              new ceres::AutoDiffCostFunction<
                  VarianceCovarianceRegularization, 1, 4, 4>(
                  new VarianceCovarianceRegularization(FLAGS_ba_vcreg_weight)));
          problem_information->addResidualBlock(
              ceres_error_terms::ResidualType::kGenericPrior, cost_function3,
              loss_function,
              {a_is_ref ? delta_0_intrinsics : cam_a.getParametersMutable(),
               cam_b.getParametersMutable()});*/
        }

        ncam_a = ncam_b;
      }
    }
  }
  size_t num_inertial_constraints_added = 0u;
  if (options.add_inertial_constraints) {
    num_inertial_constraints_added = addInertialTerms(
        options.fix_gyro_bias, options.fix_accel_bias, options.fix_velocity,
        options.gravity_magnitude, problem);
    if (num_inertial_constraints_added == 0u) {
      LOG(WARNING)
          << "WARNING: Inertial constraints enabled, but none "
          << "were found, adapting DoF settings of optimization problem...";
    }
  }

  size_t num_wheel_odometry_constraints_added = 0u;
  if (options.add_wheel_odometry_constraints) {
    num_wheel_odometry_constraints_added =
        addWheelOdometryTerms(options.fix_wheel_extrinsics, problem);
    if (num_wheel_odometry_constraints_added == 0u) {
      LOG(WARNING)
          << "WARNING: Wheel odometry constraints enabled, but none "
          << "were found, adapting DoF settings of optimization problem...";
    }
  }

  size_t num_6dof_odometry_constraints_added = 0u;
  if (options.add_6dof_odometry_constraints) {
    num_6dof_odometry_constraints_added =
        add6DoFOdometryTerms(options.fix_6dof_odometry_extrinsics, problem);
    if (num_6dof_odometry_constraints_added == 0u) {
      LOG(WARNING)
          << "WARNING: 6DoF odometry constraints enabled, but none "
          << "were found, adapting DoF settings of optimization problem...";
    }
  }

  size_t num_absolute_6dof_constraints_added = 0u;
  if (options.add_absolute_pose_constraints) {
    num_absolute_6dof_constraints_added = addAbsolutePoseConstraintsTerms(
        options.fix_absolute_pose_sensor_extrinsics, problem);
    if (num_absolute_6dof_constraints_added == 0u) {
      LOG(WARNING)
          << "WARNING: Absolute 6DoF constraints enabled, but none "
          << "were found, adapting DoF settings of optimization problem...";
    }
  }

  size_t num_lc_edges = 0u;
  if (options.add_loop_closure_edges) {
    num_lc_edges = numLoopclosureEdges(*map);
    if (num_lc_edges == 0u) {
      LOG(WARNING) << "WARNING: Loop closure edges are enabled, but none "
                   << "were found.";
    } else {
      size_t actually_added_lc_error_terms =
          augmentOptimizationProblemWithLoopclosureEdges(problem);
      CHECK(actually_added_lc_error_terms == num_lc_edges)
          << "The pose graph has " << num_lc_edges << "loop closure edges, but "
          << actually_added_lc_error_terms
          << "were added to the optimization problem!";
    }
  }

  if (options.fix_vertices) {
    LOG(INFO) << "Fixing vertex positions.";
    fixAllVerticesInProblem(problem);
  }

  // We analyze the available constraints in each mission clusters, i.e.
  // missions that are connected through either visual constraints or loop
  // closure edges, and determine the gauge fixes.
  const std::vector<vi_map::MissionIdSet>& mission_clusters =
      problem->getMissionCoobservationClusters();
  const size_t num_clusters = mission_clusters.size();
  std::vector<MissionClusterGaugeFixes> mission_cluster_gauge_fixes(
      num_clusters);
  CHECK_EQ(mission_clusters.size(), mission_cluster_gauge_fixes.size());

  for (size_t cluster_idx = 0u; cluster_idx < num_clusters; ++cluster_idx) {
    MissionClusterGaugeFixes& mission_cluster_gauge_fix =
        mission_cluster_gauge_fixes[cluster_idx];
    const vi_map::MissionIdSet& mission_cluster = mission_clusters[cluster_idx];

    const size_t cluster_num_absolute_6dof_present =
        vi_map_helpers::getNumAbsolute6DoFConstraintsForMissionCluster(
            *map, mission_cluster);
    const size_t cluster_num_absolute_6dof_used = std::min(
        cluster_num_absolute_6dof_present, num_absolute_6dof_constraints_added);
    const bool cluster_has_inertial =
        vi_map_helpers::hasInertialConstraintsInAllMissionsInCluster(
            *map, mission_cluster) &&
        (num_inertial_constraints_added > 0u);
    const bool cluster_has_visual =
        vi_map_helpers::hasVisualConstraintsInAllMissionsInCluster(
            *map, mission_cluster) &&
        (num_visual_constraints_added > 0u);
    const bool cluster_has_wheel_odometry =
        vi_map_helpers::hasWheelOdometryConstraintsInAllMissionsInCluster(
            *map, mission_cluster) &&
        (num_wheel_odometry_constraints_added > 0u);
    const bool cluster_has_6dof_odometry =
        vi_map_helpers::has6DoFOdometryConstraintsInAllMissionsInCluster(
            *map, mission_cluster) &&
        (num_6dof_odometry_constraints_added > 0u);

    // Note that if there are lc edges they always are within the cluster,
    // because otherwise the other mission would have been part of the cluster.
    const bool cluster_has_lc_edges =
        vi_map_helpers::hasLcEdgesInMissionCluster(*map, mission_cluster) &&
        (num_lc_edges > 0u);

    CHECK(
        cluster_has_inertial || cluster_has_visual ||
        cluster_has_wheel_odometry || cluster_has_6dof_odometry)
        << "Either inertial, visual or wheel odometry constraints need to be "
           "available to form a stable graph.";

    // Determine observability of scale, global position and global orientation.
    const bool scale_is_observable =
        cluster_has_inertial || cluster_has_wheel_odometry ||
        cluster_has_6dof_odometry ||
        (cluster_has_visual && cluster_num_absolute_6dof_used > 1u);

    const bool global_position_is_observable =
        cluster_num_absolute_6dof_used > 0u;

    const bool global_yaw_is_observable = cluster_num_absolute_6dof_used > 0u;

    const bool global_roll_pitch_is_observable =
        cluster_num_absolute_6dof_used > 0u || cluster_has_inertial;

    std::stringstream ss;
    ss << "\nMission Cluster: ";
    for (const vi_map::MissionId& mission_id : mission_cluster) {
      ss << "\n\t" << mission_id;
    }

    ss << "\n\nConstraints:";
    ss << "\n\tInertial constraints:\t\t"
       << ((cluster_has_inertial) ? "on" : "off");
    ss << "\n\tVisual constraints:\t\t"
       << ((cluster_has_visual) ? "on" : "off");
    ss << "\n\tWheel odometry constraints:\t"
       << ((cluster_has_wheel_odometry) ? "on" : "off");
    ss << "\n\t6DoF odometry constraints:\t"
       << ((cluster_has_6dof_odometry) ? "on" : "off");
    ss << "\n\tAbsolute 6DoF constraints:\t"
       << ((cluster_num_absolute_6dof_used > 0) ? "on" : "off");
    ss << "\n\tLoop closure edge constraints:\t"
       << ((cluster_has_lc_edges > 0) ? "on" : "off");

    ss << "\n\nIs observable:";
    ss << "\n\tScale: " << ((scale_is_observable) ? "yes" : "no");
    ss << "\n\tGlobal position: "
       << ((global_position_is_observable) ? "yes" : "no");
    ss << "\n\tGlobal yaw: " << ((global_yaw_is_observable) ? "yes" : "no");
    ss << "\n\tGlobal roll/pitch: "
       << ((global_roll_pitch_is_observable) ? "yes" : "no");
    ss << "\n";
    VLOG(1) << ss.str();

    // Apply gauge fixes for this cluster.
    mission_cluster_gauge_fix.position_dof_fixed =
        !global_position_is_observable;
    mission_cluster_gauge_fix.scale_fixed = !scale_is_observable;

    // Currently there is no scenario where yaw is observable but roll and
    // pitch, this could change if magnetometer are added. This check should
    // make sure that the logic below does not do weird stuff it add such a
    // scenario, but forget to change code here.
    CHECK(!(global_yaw_is_observable && !global_roll_pitch_is_observable));

    mission_cluster_gauge_fix.rotation_dof_fixed =
        global_roll_pitch_is_observable
            ? (global_yaw_is_observable ? FixedRotationDoF::kNone
                                        : FixedRotationDoF::kYaw)
            : FixedRotationDoF::kAll;
  }

  // Merge with already applied fixes (if necessary).
  const std::vector<MissionClusterGaugeFixes>* already_applied_cluster_fixes =
      problem->getAppliedGaugeFixesForInitialVertices();
  if (already_applied_cluster_fixes) {
    std::vector<MissionClusterGaugeFixes> merged_fixes;
    mergeGaugeFixes(
        mission_cluster_gauge_fixes, *already_applied_cluster_fixes,
        &merged_fixes);
    problem->applyGaugeFixesForInitialVertices(merged_fixes);
  } else {
    problem->applyGaugeFixesForInitialVertices(mission_cluster_gauge_fixes);
  }

  // Only case were we do NOT fix the baseframe is if there are absolute 6dof
  // constraints and the option to fix them has been disabled.
  if (num_absolute_6dof_constraints_added == 0u || options.fix_baseframes) {
    fixAllBaseframesInProblem(problem);
    VLOG(1) << "Baseframes fixed: yes";
  } else {
    VLOG(1) << "Baseframes fixed: no";
  }

  return problem;
}

}  // namespace map_optimization
