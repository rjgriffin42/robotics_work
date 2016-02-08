function zmp = plan_zmp_toe_off(l_foot_pose0, r_foot_pose0, step_plan, ds_ratio, to_ratio, hs_ratio, zmp0, t)
  
  % toe off
  toe_offset = [0.125, 0, 0, 1];
  heel_offset = [-0.125, 0, 0, 1];
  toe_off_min_angle = 0.60; % > 0 radians
  toe_off_max_angle = .122; % < pi/2 radians
  heel_strike_min_angle = 0.60; % > 0 radians
  heel_strike_max_angle = 1.22; % < pi/2 radians

  zmp = zeros(length(t), 3);
  duration = {};
  ds_poses = {};
  ss_poses = {};
  ds_position = {};
  hs_position = {};
  to_position = {};

  l_foot_pose = l_foot_pose0;
  r_foot_pose = r_foot_pose0;

  for i = 1:length(step_plan)
    duration{i} = step_plan{i}.duration;
    % compute sole poses
    if (step_plan{i}.foot == 'l')
      ds_poses{i} = l_foot_pose;
      ss_poses{i} = r_foot_pose;
      l_foot_pose = step_plan{i}.pose;
    else
      ds_poses{i} = r_foot_pose;
      ss_poses{i} = l_foot_pose;
      r_foot_pose = step_plan{i}.pose;
    end
  end
  ds_poses{end+1} = ss_poses{end};
  ss_poses{end+1} = (l_foot_pose + r_foot_pose)/2;
