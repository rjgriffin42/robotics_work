function copTrajectory = planCoPToeOff(leftFootPoseInitial, rightFootPoseInitial, stepPlan, doubleSupportRatio, toeOffRatio, heelStrikeRatio, copInitial, t)
  
  % toe off
  toeOffset = [0.125, 0, 0, 1];
  heelOffset = [-0.125, 0, 0, 1];
  toeOffMinimumAngle = 0.60; % > 0 radians
  toeOffMaximumAngle = 0.122; % < pi/2 radians
  heelStrikeMinimumAngle = 0.60; % > 0 radians
  heelStrikeMaximumAngle = 1.22; % < pi/2 radians

  copTrajectory = zeros(length(t), 3);
  duration = {};
  doubleSupportPoses = {};
  singleSupportPoses = {};
  doubleSupportPositions = {};
  heelStrikePositions = {};
  toeOffPositions = {};

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
      doubleSupportPositions{i} = copInitial;
    else
      doubleSupportPositions{i} = toeOffPositions{i-1};
    end

    singleSupportStepOffset = (singleSupportPoses{i + 1}(1:3) - ...
                              singleSupportPoses{i}(1:3))';
                          
    singleSupportHeelOffset = makehgtform('xrotate', singleSupportPoses{i}(4)) * ...
                              makehgtform('yrotate', singleSupportPoses{i}(5)) * ...
                              makehgtform('zrotate', singleSupportPoses{i}(6)) * ...
                              heelOffset';
    singleSupportHeelOffset = singleSupportHeelOffset(1:3) / ...
                              singleSupportHeelOffset(4);
                          
    cosThetaHeel = singleSupportStepOffset' * singleSupportHeelOffset / ...
                   (sqrt(singleSupportStepOffset' * singleSupportStepOffset) * ...
                    sqrt(singleSupportHeelOffset' * singleSupportHeelOffset));
    cosThetaHeel = cosThetaHeel - cos(heelStrikeMaximumAngle);
    cosThetaHeel = cosThetaHeel / (cos(heelStrikeMinimumAngle) - ...
                                   cos(heelStrikeMaximumAngle));
    cosThetaHeel = max(min(cosThetaHeel, 1), 0);
    
    heelStrikePositions{i} = singleSupportPoses{i}(1:3) +...
        heelStrikeRatio * cosThetaHeel * singleSupportHeelOffset';
    
    singleSupportToeOffset = makehgtform('xrotate', singleSupportPoses{i}(4)) * ...
                             makehgtform('yrotate', singleSupportPoses{i}(5)) * ...
                             makehgtform('zrotate', singleSupportPoses{i}(6)) * ...
                             toeOffset';
    singleSupportToeOffset = singleSupportToeOffset(1:3) / singleSupportToeOffset(4);
    cosThetaToe = singleSupportStepOffset' * singleSupportToeOffset / ...
        (sqrt(singleSupportStepOffset' * singleSupportStepOffset) * ...
        sqrt(singleSupportToeOffset' * singleSupportToeOffset));
    cosThetaToe = cosThetaToe - cos(toeOffMaximumAngle);
    cosThetaToe = cosThetaToe / (cos(toeOffMinimumAngle) - cos(toeOffMaximumAngle));
    cosThetaToe = max(min(cosThetaToe, 1), 0);
    
    toeOffPositions{i} = singleSupportPoses{i}(1:3) + ...
                        toeOffRatio * cosThetaToe * singleSupportToeOffset';
  end

  timeInitial = 0;
  for i = 1:length(stepPlan)
    doubleSupportIndex = (t >= timeInitial & ...
                          t < timeInitial + doubleSupportRatio * duration{i});
    singleSupportIndex = (t >= timeInitial + doubleSupportRatio * duration{i} & ...
                          t < timeInitial + duration{i});

    for j = 1:3
      xInitial = [doubleSupportPositions{i}(j) 0];
      xFinal = [heelStrikePositions{i}(j) 0];
      copTrajectories = computeCubicSplineTrajectory(xInitial, xFinal, t(doubleSupportIndex));
      copTrajectory(doubleSupportIndex, j) = copTrajectories(:, 1);
    end
    for j = 1:3
      xInitial = [heelStrikePositions{i}(j) 0];
      xFinal = [toeOffPositions{i}(j) 0];
      copTrajectories = computeCubicSplineTrajectory(xInitial, xFinal, t(singleSupportIndex));
      copTrajectory(singleSupportIndex, j) = copTrajectories(:, 1);
    end
    timeInitial = timeInitial + duration{i};
  end
  copTrajectory(end, :) = toeOffPositions{i}(1:3);
end
