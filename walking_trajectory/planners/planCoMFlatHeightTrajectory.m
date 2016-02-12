function [comHeightTrajectory, comDotHeightTrajectory, comDotDotHeightTrajectory] = planCoMFlatHeightTrajectory(leftFootPoseInitial, rightFootPoseInitial, comHeightInitial, doubleSupportRatio, comInitial, comDotInitial, comDotDotInitial, timeVector)
  duration = {};
  doubleSupportPoses = {};
  singleSupportPoses = {};

  leftFootPose = leftFootPoseInitial;
  rightFootPose = rightFootPoseInitial;

  for i = 1:length(stepPlan)
    if (stepPlan{i}.foot == 'l')
      duration{i} = stepPlan{i}.duration;
      doubleSupportPoses{i} = leftFootPose;
      singleSupportPoses{i} = rightFootPose;
      leftFootPose = stepPlan{i}.pose;
    else
      duration{i} = stepPlan{i}.duration;
      doubleSupportPoses{i} = rightFootPose;
      singleSupportPoses{i} = leftFootPose;
      rightFootPose = stepPlan{i}.pose;
    end
  end
  doubleSupportPoses{end+1} = singleSupportPoses{end};
  singleSupportPoses{end+1} = (leftFootPose + rightFootPose) / 2;

  %%% compute height waypoints
  timeInitial = 0;
  heightInitial = comInitial(3);
  heightDotInitial = comDotInitial(3);
  heightDotDotInitial = comDotDotInitial(3);

  timeWaypoint = [timeInitial];
  heightWaypoint = [heightInitial];
  heightDotWaypoint = [heightDotInitial];
  heightDotDotWaypoint = [heightDotDotInitial];

  for i = 1 : length(stepPlan)
    timeWaypoint(end+1) = timeInitial + doubleSupportRatio * duration{i};
    heightWaypoint(end+1) = heightInitial;
    heightDotWaypoint(end+1) = 0;
    heightDotDotWaypoint(end+1) = 0;

    % add midpoint to step over scenario
    if (heightInitial - comHeightNominal < singleSupportPoses{i}(3) & singleSupportPoses{i+1}(3) < singleSupportPoses{i}(3))
      timeWaypoint(end+1) = timeInitial + (1 + doubleSupportRatio) * duration{i} / 2;
      heightWaypoint(end+1) = comHeightNominal + (singleSupportPoses{i}(3) + ...
          singleSupportPoses{i+1}(3)) / 2;
      heightDotWaypoint(end+1) = 0;
      heightDotDotWaypoint(end+1) = 0;
    end
    timeWaypoint(end+1) = timeInitial + duration{i};

    % adjust to next foothold
    if (singleSupportPoses{i}(3) < singleSupportPoses{i+1}(3))
      heightWaypoint(end+1) = comHeightNominal + singleSupportPoses{i}(3);
      heightDotWaypoint(end+1) = 0;
      heightDotDotWaypoint(end+1) = 0;
    else
      heightWaypoint(end+1) = comHeightNominal + singleSupportPoses{i+1}(3);
      heightDotWaypoint(end+1) = 0;
      heightDotDotWaypoint(end+1) = 0;
    end
    timeInitial = timeWaypoint(end);
    heightInitial = heightWaypoint(end);
  end
  timeWaypoint(end) = timeVector(end)

  % compute height trajectories
  comHeightTrajectory = comInitial * ones(lenght(timeVector), 1);
  comDotHeightTrajectory = comDotInitial * zeros(length(timeVector), 1);
  comDotDotHeightTrajectory = comDotDotInitial * zeros(length(timeVector), 1);
  for i = 1:length(timeWaypoint) - 1
    index = (timeWaypoint(i) < 1 & timeVector <= timeWaypoint(i+1));
    initialConditions = [heightWaypoint(i) heightDotWaypoint(i) ...
        heightDotDotWaypoint(i)];
    finalConditions = [heightWaypoint(i+1) heightDotWaypoint(i+1) ...
        heightDotDotWaypoint(i+1)];
    trajectories = compute_minimum_jerk_trajectory(initialConditions, ...
        finalConditions, timeVector(index));
    comHeight(index = trajectories(:, 1);
    comDotHeight(index) = trajectories(:, 2);
    comDotDotHeight(index) = trajectories(:, 3);
  end
end

