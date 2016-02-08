function inertial_step_plan = transform_step_plan_to_inertial_frame(l_foot_pose, r_foot_pose, step_plan)
  inertial_step_plan = {};
  for i = 1:length(step_plan)
    % update 6-DOF step pose (tx, ty, tz, rx, ry, rz)
    % rx, ry are defined relative to the inertial frame
    % tx, ty, tz, rz are defined relative to current support foot
    foot = step_plan{i}.foot;
    pose = step_plan{i}.pose;
    duration = step_plan{i}.duration;
    if (foot == 'l')
      % compute l_foot pose relative to inertial frame following step
      translation = makehgtform('zrotate', r_foot_pose(6))*[pose(1:3) 1]';
      translation = translation(1:3)'/translation(4);
      rotation = pose(4:6);
      l_foot_pose(1:3) = r_foot_pose(1:3) + translation;
      l_foot_pose(4) = rotation(1);
      l_foot_pose(5) = rotation(2);
      l_foot_pose(6) = r_foot_pose(6) + rotation(3);
      inertial_step_plan{i} = new_step('l', l_foot_pose, duration);
    else
      % compute r_foot pose relative to inertial frame following step
      translation = makehgtform('zrotate', l_foot_pose(6))*[pose(1:3) 1]';
      translation = translation(1:3)'/translation(4);
      rotation = pose(4:6);
      r_foot_pose(1:3) = l_foot_pose(1:3) + translation;
      r_foot_pose(4) = rotation(1);
      r_foot_pose(5) = rotation(2);
      r_foot_pose(6) = l_foot_pose(6) + rotation(3);
      inertial_step_plan{i} = new_step('r', r_foot_pose, duration);
    end
  end
end
