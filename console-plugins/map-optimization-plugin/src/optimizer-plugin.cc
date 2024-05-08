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
DEFINE_double(
    camera_segment_duration, 60.0,
    "Duration of the camera segments in seconds.");

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
      {"create_camera_segments"},
      [this]() -> int { return createCameraSegments(); },
      "This command creates new cameras after --camera_segment_duration "
      "seconds",
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

int OptimizerPlugin::createCameraSegments() {
  std::string selected_map_key;
  if (!getSelectedMapKeyIfSet(&selected_map_key)) {
    return common::kStupidUserError;
  }
  if (FLAGS_camera_segment_duration < 1) {
    LOG(ERROR) << "The segment duration is too short.";
    return common::kStupidUserError;
  }

  vi_map::VIMapManager map_manager;
  vi_map::VIMapManager::MapWriteAccess map =
      map_manager.getMapWriteAccess(selected_map_key);

  vi_map::SensorManager& sensor_manager = map->getSensorManager();
  if (!sensor_manager.hasSensorOfType(vi_map::SensorType::kNCamera)) {
    LOG(INFO) << "Map does not contain any camera.";
    return common::kSuccess;
  }

  vi_map::MissionIdList mission_id_list;
  map->getAllMissionIds(&mission_id_list);

  for (const vi_map::MissionId& mission_id : mission_id_list) {
    if (map->numVerticesInMission(mission_id) == 0) {
      LOG(INFO) << "Mission " << mission_id.shortHex()
                << " contains no vertices. Skipping.";
      continue;
    }
    vi_map::VIMission& mission = map->getMission(mission_id);
    if (!mission.hasNCamera()) {
      LOG(INFO) << "Mission " << mission_id.shortHex()
                << " contains no camera. Skipping.";
      continue;
    }
    if (!mission.getSegmentNCameraIds().empty()) {
      // TODO(sscholbe): How should we handle calling the command multiple
      // times?
      LOG(INFO) << "Mission " << mission_id.shortHex()
                << " already contains camera segments. Skipping.";
      continue;
    }

    // Use the original mission camera as our reference
    const aslam::NCamera::Ptr reference_camera =
        map->getMissionNCameraPtr(mission_id);

    pose_graph::VertexIdList vertex_ids;
    map->getAllVertexIdsInMissionAlongGraph(mission_id, &vertex_ids);

    // Represents the timestamps for a segment in nanoseconds, excluding the end
    typedef std::pair<size_t, size_t> Segment;

    // Compute the segments based on the desired duration
    std::vector<Segment> segments;
    size_t segment_duration_ns =
        aslam::time::secondsToNanoSeconds(FLAGS_camera_segment_duration);
    size_t start_time_ns =
        map->getVertex(vertex_ids.front()).getMinTimestampNanoseconds();
    size_t end_time_ns =
        map->getVertex(vertex_ids.back()).getMinTimestampNanoseconds();
    for (size_t time_ns = start_time_ns; time_ns <= end_time_ns;
         time_ns += segment_duration_ns) {
      segments.push_back(
          std::make_pair(time_ns, time_ns + segment_duration_ns));
    }

    std::vector<aslam::NCamera::Ptr> cameras;

    // The first segment will use the reference camera
    cameras.push_back(reference_camera);

    // The other segments will have their copy of the reference camera
    for (size_t i = 1; i < segments.size(); i++) {
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
      cameras.push_back(sensor_manager.getSensorPtr<aslam::NCamera>(clone_id));
      mission.addSegmentNCameraId(clone_id);
    }

    // Traverse the vertices and assign them to the corresponding segment camera
    std::vector<Segment>::iterator current_segment = segments.begin();
    std::vector<aslam::NCamera::Ptr>::iterator current_camera = cameras.begin();
    for (const pose_graph::VertexId& vertex_id : vertex_ids) {
      vi_map::Vertex& vertex = map->getVertex(vertex_id);
      size_t time_ns = vertex.getMinTimestampNanoseconds();
      while (time_ns >= current_segment->second) {
        current_segment++;
        current_camera++;
      }
      vertex.setNCameras(*current_camera);
    }

    LOG(INFO) << "Created " << segments.size() << " segments for mission "
              << mission_id.shortHex();
  }

  return common::kSuccess;
}

}  // namespace map_optimization_plugin

MAPLAB_CREATE_CONSOLE_PLUGIN_WITH_PLOTTER(
    map_optimization_plugin::OptimizerPlugin);
