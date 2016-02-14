function inertialStepPlan = transformStepPlanToInertialFrame(leftFootPose, rightFootPose, stepPlan)
  inertialStepPlan = {};
  for i = 1:length(stepPlan)
    % update 6-DOF step pose (tx, ty, tz, rx, ry, rz)
    % rx, ry are defined relative to the inertial frame
    % tx, ty, tz, rz are defined relative to current support foot
    foot = stepPlan{i}.foot;
    pose = stepPlan{i}.pose;
    duration = stepPlan{i}.duration;
    if (foot == 'l')
      % compute l_foot pose relative to inertial frame following step
      translation = makehgtform('zrotate', rightFootPose(6)) * [pose(1:3) 1]';
      translation = translation(1:3)' / translation(4);
      rotation = pose(4:6);
      leftFootPose(1:3) = rightFootPose(1:3) + translation;
      leftFootPose(4) = rotation(1);
      leftFootPose(5) = rotation(2);
      leftFootPose(6) = rightFootPose(6) + rotation(3);
      inertialStepPlan{i} = newStep('l', leftFootPose, duration);
    else
      % compute r_foot pose relative to inertial frame following step
      translation = makehgtform('zrotate', leftFootPose(6))*[pose(1:3) 1]';
      translation = translation(1:3)'/translation(4);
      rotation = pose(4:6);
      rightFootPose(1:3) = leftFootPose(1:3) + translation;
      rightFootPose(4) = rotation(1);
      rightFootPose(5) = rotation(2);
      rightFootPose(6) = leftFootPose(6) + rotation(3);
      inertialStepPlan{i} = newStep('r', rightFootPose, duration);
    end
  end
end
