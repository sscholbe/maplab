#include "map-optimization-plugin/optimizer-plugin.h"

#include <ceres/ceres.h>
#include <console-common/basic-console-plugin.h>
#include <console-common/console.h>
#include <map-manager/map-manager.h>
#include <map-optimization/outlier-rejection-solver.h>
#include <map-optimization/solver-options.h>
#include <map-optimization/vi-map-optimizer.h>
#include <map-optimization/vi-map-relaxation.h>
#include <map-optimization/vi-optimization-builder.h>
#include <vi-map-helpers/vi-map-queries.h>
#include <vi-map/vi-map.h>
#include <visualization/viwls-graph-plotter.h>

DECLARE_string(map_mission);
DECLARE_string(map_mission_list);
DEFINE_string(
    relax_external_loop_closures_file, "",
    "yaml file containing the loop "
    "closures for the pose graph relaxation.");
DEFINE_string(
    dump_camera_file, "cams.txt", "File where to store the camera parameters.");

namespace map_optimization_plugin {
OptimizerPlugin::OptimizerPlugin(
    common::Console* console, visualization::ViwlsGraphRvizPlotter* plotter)
    : common::ConsolePluginBaseWithPlotter(console, plotter) {
  addCommand(
      {"optimize_visual", "optv"},
      [this]() -> int {
        map_optimization::ViProblemOptions options =
            map_optimization::ViProblemOptions::initFromGFlags();

        // The only difference to the default behaviour of the optimization is
        // to disable inertial constraints.
        options.add_inertial_constraints = false;

        return optimize(options);
      },
      "Visual optimization over the selected missions "
      "(per default all).",
      common::Processing::Sync);
  addCommand(
      {"optimize_visual_inertial", "optvi"},
      [this]() -> int {
        map_optimization::ViProblemOptions options =
            map_optimization::ViProblemOptions::initFromGFlags();

        // No special flags need to be set, since this is the default behaviour
        // of the optimization.

        return optimize(options);
      },
      "Visual-inertial optimization over the selected missions "
      "(per default all).",
      common::Processing::Sync);
  addCommand(
      {"relax"}, [this]() -> int { return relaxMap(); }, "nRelax posegraph.",
      common::Processing::Sync);
  addCommand(
      {"relax_missions_independently", "srelax"},
      [this]() -> int { return relaxMapMissionsSeparately(); },
      "Relax missions separately.", common::Processing::Sync);
  addCommand(
      {"relax_external"},

      [this]() -> int {
        std::string selected_map_key;

        if (!getSelectedMapKeyIfSet(&selected_map_key)) {
          return common::kStupidUserError;
        }

        if (FLAGS_relax_external_loop_closures_file == "") {
          LOG(ERROR) << "You need to provide the yaml file containing the loop "
                        "closures. Usage: --external_loop_closures_file "
                        "/path/to/file.yaml";
          return common::kStupidUserError;
        }

        vi_map::VIMapManager map_manager;
        vi_map::VIMapManager::MapWriteAccess map =
            map_manager.getMapWriteAccess(selected_map_key);

        const data_import_export::LoopClosureEdges edges =
            data_import_export::loopClosureEdgesFromYamlFile(
                FLAGS_relax_external_loop_closures_file);
        addExternalLoopClosureEdgesToMap(edges, map.get());
        relaxMap(map.get());
        return common::kSuccess;
      },

      "This command runs pose graph relaxation with loop closure constraints "
      "provided in a yaml file. After relaxation, the loop closure constraints "
      "are removed from the map again. Set the path to the yaml file "
      "containing the loop closure constraints with the flag "
      "--external_loop_closures_file",
      common::Processing::Sync);
  addCommand(
      {"create_camera_drift"}, [this]() -> int { return createCameraDrift(); },
      "Creates the cameras for a drift.", common::Processing::Sync);
  addCommand(
      {"dump_cameras"}, [this]() -> int { return dumpCameras(); },
      "Dumps the camera parameters into --dump_camera_file.",
      common::Processing::Sync);
}

void OptimizerPlugin::addExternalLoopClosureEdgesToMap(
    const data_import_export::LoopClosureEdges& edges, vi_map::VIMap* map) {
  CHECK_NOTNULL(map);

  vi_map_helpers::VIMapQueries vi_map_queries(*map);

  constexpr uint64_t kTimestampDifferenceToleranceNs =
      static_cast<uint64_t>(20e6);  // 20 ms

  for (const data_import_export::LoopClosureEdge& edge : edges) {
    pose_graph::VertexId id_from;
    uint64_t timestamp_difference_from;
    const bool found_from = vi_map_queries.getClosestVertexIdByTimestamp(
        edge.from.timestamp_ns, kTimestampDifferenceToleranceNs, &id_from,
        &timestamp_difference_from);

    pose_graph::VertexId id_to;
    uint64_t timestamp_difference_to;
    const bool found_to = vi_map_queries.getClosestVertexIdByTimestamp(
        edge.to.timestamp_ns, kTimestampDifferenceToleranceNs, &id_to,
        &timestamp_difference_to);

    if (!found_from || !found_to) {
      std::string msg = "Could not find a close enough vertex for:\n";
      std::string camera_msg = "Camera pose with timestamp ";

      if (!found_from) {
        msg += camera_msg + std::to_string(edge.from.timestamp_ns) + "\n";
      }
      if (!found_to) {
        msg += camera_msg + std::to_string(edge.to.timestamp_ns) + "\n";
      }
      LOG(ERROR) << msg << "Skipping this loop closure.";
      continue;
    }

    LOG(INFO) << "Found vertices corresponding to camera poses. \n"
                 "Timestamp differences in Nanoseconds: camera from - "
              << timestamp_difference_from << " camera_to - "
              << timestamp_difference_to << std::endl;
    LOG(INFO) << "Adding loop closure from vertex " << id_from << " to vertex "
              << id_to << std::endl;

    const vi_map::Vertex& vertex_source = map->getVertex(id_from);
    const vi_map::Vertex& vertex_target = map->getVertex(id_to);

    // Loop closure constraints are not necessarily defined on poses existing in
    // the pose graph. This function calculates the transformation of the loop
    // closure constraint between two vertex poses given two poses and the
    // corresponding loop closure transformation.
    const pose::Transformation T_loop_closure = adaptTransformation(
        edge.from.pose, edge.to.pose, vertex_source.get_T_M_I(),
        vertex_target.get_T_M_I(), edge.T_from_to);

    // Add loop closure edge to map.
    pose_graph::EdgeId edge_id;
    aslam::generateId(&edge_id);
    CHECK(edge_id.isValid());

    vi_map::Edge::UniquePtr loop_closure_edge(new vi_map::LoopClosureEdge(
        edge_id, id_from, id_to, edge.switch_variable,
        edge.switch_variable_variance, T_loop_closure, edge.covariance));
    map->addEdge(std::move(loop_closure_edge));
  }
}

void OptimizerPlugin::relaxMap(vi_map::VIMap* map) {
  CHECK_NOTNULL(map);
  map_optimization::VIMapRelaxation relaxation(
      getPlotterUnsafe(), kSignalHandlerEnabled);
  vi_map::MissionIdList mission_id_list;
  map->getAllMissionIds(&mission_id_list);

  vi_map::MissionIdSet mission_ids(
      mission_id_list.begin(), mission_id_list.end());
  ceres::Solver::Options solver_options =
      map_optimization::initSolverOptionsFromFlags();
  relaxation.solveRelaxation(solver_options, mission_ids, map);
}

pose::Transformation OptimizerPlugin::adaptTransformation(
    const pose::Transformation& T_M_S, const pose::Transformation& T_M_T,
    const pose::Transformation& T_M_S2, const pose::Transformation& T_M_T2,
    const pose::Transformation& T_S_T) {
  const pose::Transformation T_S2_S = T_M_S2.inverse() * T_M_S;
  const pose::Transformation T_T_T2 = T_M_T.inverse() * T_M_T2;
  const pose::Transformation T_S2_T2 = T_S2_S * T_S_T * T_T_T2;
  return T_S2_T2;
}

int OptimizerPlugin::optimize(
    const map_optimization::ViProblemOptions& options) {
  // Select map and missions to optimize.
  std::string selected_map_key;
  if (!getSelectedMapKeyIfSet(&selected_map_key)) {
    return common::kStupidUserError;
  }
  vi_map::VIMapManager map_manager;
  vi_map::VIMapManager::MapWriteAccess map =
      map_manager.getMapWriteAccess(selected_map_key);

  vi_map::MissionIdList missions_to_optimize_list;
  if (!FLAGS_map_mission.empty()) {
    if (!FLAGS_map_mission_list.empty()) {
      LOG(ERROR) << "Please provide only one of --map_mission and "
                 << "--map_mission_list.";
      return common::kStupidUserError;
    }
    vi_map::MissionId mission_id;
    if (!map->hexStringToMissionIdIfValid(FLAGS_map_mission, &mission_id)) {
      LOG(ERROR) << "The given mission id \"" << FLAGS_map_mission
                 << "\" is not valid.";
      return common::kStupidUserError;
    }
    missions_to_optimize_list.emplace_back(mission_id);
  } else if (!FLAGS_map_mission_list.empty()) {
    if (!vi_map::csvIdStringToIdList(
            FLAGS_map_mission_list, &missions_to_optimize_list)) {
      LOG(ERROR) << "The provided CSV mission id list is not valid!";
      return common::kStupidUserError;
    }
  } else {
    map->getAllMissionIds(&missions_to_optimize_list);
  }
  vi_map::MissionIdSet missions_to_optimize(
      missions_to_optimize_list.begin(), missions_to_optimize_list.end());

  map_optimization::VIMapOptimizer optimizer(
      getPlotterUnsafe(), kSignalHandlerEnabled);
  bool success = optimizer.optimize(options, missions_to_optimize, map.get());

  if (!success) {
    return common::kUnknownError;
  }
  return common::kSuccess;
}

int OptimizerPlugin::relaxMap() {
  std::string selected_map_key;
  if (!getSelectedMapKeyIfSet(&selected_map_key)) {
    return common::kStupidUserError;
  }

  map_optimization::VIMapRelaxation relaxation(
      getPlotterUnsafe(), kSignalHandlerEnabled);
  vi_map::VIMapManager map_manager;
  vi_map::VIMapManager::MapWriteAccess map =
      map_manager.getMapWriteAccess(selected_map_key);

  vi_map::MissionIdList mission_id_list;
  map.get()->getAllMissionIds(&mission_id_list);

  relaxation.findLoopClosuresAndSolveRelaxation(mission_id_list, map.get());

  return common::kSuccess;
}

int OptimizerPlugin::relaxMapMissionsSeparately() {
  std::string selected_map_key;
  if (!getSelectedMapKeyIfSet(&selected_map_key)) {
    return common::kStupidUserError;
  }

  map_optimization::VIMapRelaxation relaxation(
      getPlotterUnsafe(), kSignalHandlerEnabled);
  vi_map::VIMapManager map_manager;
  vi_map::VIMapManager::MapWriteAccess map =
      map_manager.getMapWriteAccess(selected_map_key);

  vi_map::MissionIdList mission_id_list;
  map.get()->getAllMissionIds(&mission_id_list);

  for (const vi_map::MissionId& mission_id : mission_id_list) {
    relaxation.findLoopClosuresAndSolveRelaxation({mission_id}, map.get());
  }

  return common::kSuccess;
}

int OptimizerPlugin::createCameraDrift() {
  std::string selected_map_key;
  if (!getSelectedMapKeyIfSet(&selected_map_key)) {
    return common::kStupidUserError;
  }

  vi_map::VIMapManager map_manager;
  vi_map::VIMapManager::MapWriteAccess map =
      map_manager.getMapWriteAccess(selected_map_key);

  vi_map::MissionIdList mission_id_list;
  map.get()->getAllMissionIds(&mission_id_list);

  vi_map::SensorManager& sensor_manager = map->getSensorManager();

  for (const vi_map::MissionId& mission_id : mission_id_list) {
    if (map->numVerticesInMission(mission_id) < 2) {
      LOG(INFO) << "Mission " << mission_id.shortHex()
                << " contains less than two vertices. Skipping.";
      continue;
    }
    vi_map::VIMission& mission = map->getMission(mission_id);
    if (!mission.hasNCamera()) {
      LOG(INFO) << "Mission " << mission_id.shortHex()
                << " contains no camera. Skipping.";
      continue;
    }

    const aslam::NCamera::Ptr reference_camera =
        map->getMissionNCameraPtr(mission_id);

    pose_graph::VertexIdList vertex_ids;
    map->getAllVertexIdsInMissionAlongGraph(mission_id, &vertex_ids);

    bool first_vertex = true;
    for (const pose_graph::VertexId& vertex_id : vertex_ids) {
      if (first_vertex) {
        first_vertex = false;
        continue;
      }
      vi_map::Vertex& vertex = map->getVertex(vertex_id);

      aslam::NCamera::UniquePtr clone_ptr(reference_camera->cloneWithNewIds());
      aslam::SensorId clone_id = clone_ptr->getId();
      if (sensor_manager.isBaseSensor(reference_camera->getId())) {
        // TODO(sscholbe): Can a camera be a base sensor?
        sensor_manager.addSensorAsBase<aslam::NCamera>(std::move(clone_ptr));
      } else {
        sensor_manager.addSensor<aslam::NCamera>(
            std::move(clone_ptr),
            sensor_manager.getBaseSensorId(reference_camera->getId()),
            sensor_manager.getSensor_T_B_S(reference_camera->getId()));
      }
      mission.drift_ncamera_ids.push_back(clone_id);

      // Since the cameras only model a drift, their parameters are a delta
      // Set intrinsics and distortion parameters initially to 0
      aslam::NCamera::Ptr ncamera =
          sensor_manager.getSensorPtr<aslam::NCamera>(clone_id);
      for (size_t cam_idx = 0u; cam_idx < ncamera->numCameras(); ++cam_idx) {
        aslam::Camera& cam = ncamera->getCameraMutable(cam_idx);
        std::fill(
            cam.getParametersMutable(),
            cam.getParametersMutable() + cam.getParameterSize(), 0.0);
        if (cam.getDistortion().getType() !=
            aslam::Distortion::Type::kNoDistortion) {
          std::fill(
              cam.getDistortionMutable()->getParametersMutable(),
              cam.getDistortionMutable()->getParametersMutable() +
                  cam.getDistortion().getParameterSize(),
              0.0);
        }
      }
      vertex.setNCameras(sensor_manager.getSensorPtr<aslam::NCamera>(clone_id));
    }
  }

  return common::kSuccess;
}

int OptimizerPlugin::dumpCameras() {
  std::string selected_map_key;
  if (!getSelectedMapKeyIfSet(&selected_map_key)) {
    return common::kStupidUserError;
  }
  if (FLAGS_dump_camera_file.empty()) {
    LOG(ERROR) << "Path not specified.";
    return common::kStupidUserError;
  }

  std::ofstream of(FLAGS_dump_camera_file);
  of << "# mission_id ncamera_id camera_idx width height intrinsics* "
        "distortion*"
     << std::endl;

  vi_map::VIMapManager map_manager;
  vi_map::VIMapManager::MapWriteAccess map =
      map_manager.getMapWriteAccess(selected_map_key);

  vi_map::SensorManager& sensor_manager = map->getSensorManager();

  vi_map::MissionIdList mission_id_list;
  map->getAllMissionIds(&mission_id_list);

  for (const vi_map::MissionId& mission_id : mission_id_list) {
    vi_map::VIMission& mission = map->getMission(mission_id);
    if (!mission.hasNCamera()) {
      LOG(INFO) << "Mission " << mission_id.shortHex()
                << " contains no camera. Skipping.";
      continue;
    }

    std::vector<aslam::NCamera::Ptr> cameras;
    const aslam::NCamera::Ptr reference_camera =
        map->getMissionNCameraPtr(mission_id);
    cameras.push_back(reference_camera);
    for (const aslam::SensorId& ncamera_id : mission.drift_ncamera_ids) {
      cameras.push_back(
          sensor_manager.getSensorPtr<aslam::NCamera>(ncamera_id));
    }

    for (const aslam::NCamera::Ptr& ncam : cameras) {
      for (size_t idx = 0; idx < ncam->getNumCameras(); idx++) {
        const aslam::Camera& cam = ncam->getCamera(idx);
        of << std::fixed << std::setprecision(0);
        of << mission_id.shortHex() << " " << ncam->getId().shortHex() << " "
           << idx << std::setprecision(8);
        // intrinsics
        of << " " << cam.imageWidth() << " " << cam.imageHeight();
        const Eigen::VectorXd& intrinsics = cam.getParameters();
        for (size_t i = 0; i < intrinsics.size(); i++) {
          of << " " << intrinsics(i);
        }
        // distortion
        const Eigen::VectorXd& dist = cam.getDistortion().getParameters();
        for (size_t i = 0; i < dist.size(); i++) {
          of << " " << dist(i);
        }
        of << std::endl;
      }
    }
  }

  return common::kSuccess;
}

}  // namespace map_optimization_plugin

MAPLAB_CREATE_CONSOLE_PLUGIN_WITH_PLOTTER(
    map_optimization_plugin::OptimizerPlugin);
