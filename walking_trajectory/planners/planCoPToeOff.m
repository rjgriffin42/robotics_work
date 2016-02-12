function zmpTrajectory = planCoPToeOff(leftFootPoseInitial, rightFootPoseInitial, stepPlan, doubleSupportRatio, toeOffRatio, heelStrikeRatio, copInitial, t)
  
  % toe off
  toeOffset = [0.125, 0, 0, 1];
  heelOffset = [-0.125, 0, 0, 1];
  toeOffMinimumAngle = 0.60; % > 0 radians
  toeOffMaximumAngle = 0.122; % < pi/2 radians
  heelStrikeMinimumAngle = 0.60; % > 0 radians
  heelStrikeMaximumAngle = 1.22; % < pi/2 radians

  zmpTrajectory = zeros(length(t), 3);
  duration = {};
  doubleSupportPoses = {};
  singleSupportPoses = {};
  doubleSupportPostion = {};
  heelStrikePosition = {};
  toeOffPosition = {};

  leftFootPose = leftFootPoseInitial;
  rightFootPose = rightFootPoseInitial;

  for i = 1:length(stepPlan)
    duration{i} = stepPlan{i}.duration;
    % compute sole poses
    if (stepPlan{i}.foot == 'l')
      doubleSupportPoses{i} = leftFootPose;
      singleSupportPoses{i} = rightFootPose;
      leftFootPose = stepPlan{i}.pose;
    else
      doubleSupportPoses{i} = rightFootPose;
      singleSupportPoses{i} = leftFootPose;
      rightFootPose = stepPlan{i}.pose;
    end
  end
  doubleSupportPoses{end+1} = singleSupportPoses{end};
  singleSupportPoses{end+1} = (leftFootPose + rightFootPose)/2;

  for i = 1:length(stepPlan)
    % compute toe positions
    if i == 1
      doubleSupportPostions{i} = copInitial;
    else
      doubleSupportPostions{i} = toeOffPosition{i-1};
    end

    singleSupportStepOffset = (singleSupportPoses{i + 1}(1:3) -
                              singleSupportPoses{i}(1:3))';
    singleSupportHeelOffset = makehgtform('xrotate', singleSupportPoses{i}(4)) * ...
                              makehgtform('yrotate', singleSupportPoses{i}(5)) * ...
                              makehgtform('zrotate', singleSupportPoses{i}(6)) * ...
                              heelOffset';
    singleSupportHeelOffset = singleSupportHeelOffset(1:3) / ...
                              singleSupportHeelOffset(4);
    cosThetaHeel = singleSupportStepOffset' * singleSupportHeelOffset / ...
                   (sqrt(singleSupportStepOffset' * singleSupportStepOffset) * ...
                   (sqrt(singleSupportHeelOffset' * singleSupportHeelOffset);
    cosThetaHeel = cosThetaHeel - cos(heelStrikeMaximumAngle);
    cosThetaHeel = cosThetaHeel / (cos(heelStrikeMinimumAngle) - ...
                                   cos(heelStrikeMaximumAngle);
    cosThetaHeel = max(min(cosThetaHeel, 1), 0);
    toeOffPosition{i} = singleSupportPoses{i}(1:3) + ...
                        toeOffRatio * cosThetaToe * singleSupportToeOffset';
  end

  timeInitial = 0;
  for i = 1:length(stepPlan)
    doubleSupportIndex = (t >= timeInitial & ...
                          t < timeInitial + doubleSupportRatio * duration{i});
    singleSupportIndex = (t >= timeInitial + doubleSupportRatio * duration{i} & ...
                          t < timeInitial + duration{i});

    for j = 1:3
      xInitial = [doubleSupportPosition{i}(j) 0];
      xFinal = [heelStrikePosition{i}(j) 0];
      X = computeCubicSplineTrajectory(xInitial, xFinal, t(doubleSupportIndex));
      zmpTrajectory(doubleSupportIndex, j) = X(:, 1);
    end
    for j = 1:3
      xInitial = [heelStrikePosition{i}(j) 0];
      xFinal = [toeOffPosition{i}(j) 0];
      X = computeCubicSplineTrajectory(xInitial, xFinal, t(doubleSupportIndex));
      zmpTrajectory(singleSupportIndex, j) = X(:,1);
    end
    tInitial = tInitial + duration{i};
  end
  zmpTrajectory(end, :) = toeOffPosition{i}(1:3);
end
