function [dcmTrajectory, dcmDotTrajectory, vrpTrajectory, copTrajectory] = ...
    planDynamicDCMSpline(cmpTrajectory, leftFootPose, rightFootPose, footstepPlan, ...
    omegaTrajectory, omegaDotTrajectory, dcmInitial, dcmDotInitial,...
    numberOfKnots, plannerDT, Q, R, V, D)

  gravity = 9.81;
  timeVector = footstepPlan.timeVector;

  % compute vrp reference
  zInitial = 9.81 * ones(length(timeVector), 1) ./ (omegaTrajectory.^2 - ...
      omegaDotTrajectory);
  vrpTrajectory = cmpTrajectory + [zeros(length(timeVector), 1) ...
      zeros(length(timeVector), 1) zInitial];

  stepPlan = footstepPlan.stepPlan;
  % compute final DCM position
  for i = 1:length(stepPlan)
    duration{i} = stepPlan{i}.duration;
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
  singleSupportPoses{end+1} = (leftFootPose + rightFootPose) / 2;
  
  % assemble spline coefficient matrix
  totalTime = 0;
  for step = 1:length(stepPlan)
      totalTime = footstepPlan.doubleSupportDuration(step) + totalTime;
      totalTime = footstepPlan.singleSupportDuration(step) + totalTime;
  end
  neededKnots = floor((totalTime / plannerDT) / (numberOfKnots - 1));
  numberOfKnots = (totalTime / plannerDT) / neededKnots + 1;
  segmentTime = totalTime / (numberOfKnots - 1);
  for segment = 1:numberOfKnots-1
      segmentDuration(segment) = segmentTime;
  end
  [betaMatrices, betaDotMatrices, betaDotDotMatrices] = ...
      computeCubicProjection(segmentDuration);

  % create discretized projection
  rowIndex = 0;
  for segment = 1:numberOfKnots-1
      if segment == 1
          t = 0:plannerDT:segmentDuration(segment);
      else
          t = plannerDT:plannerDT:segmentDuration(segment);
      end
      T = [t'.^0, t'.^1, t'.^2, t'.^3];
      Phi(rowIndex+1:rowIndex+size(T,1),:) = T * betaMatrices{segment};
      PhiDot(rowIndex+1:rowIndex+size(T,1),:) = T(:,1:3) * betaDotMatrices{segment};
      PhiDotDot(rowIndex+1:rowIndex+size(T,1),:) = T(:,1:2) * betaDotDotMatrices{segment};
      rowIndex = rowIndex + size(T,1);
  end

     % add constraints on initial and final value
     CE = zeros(2, 2*numberOfKnots);
     CE(1,1) = 1;
     CE(2,numberOfKnots) = 1;
     ce(1,:) = dcmInitial(1:2);
     ce(2,:) = [singleSupportPoses{end}(1) singleSupportPoses{end}(2)];
     
     G00 = R * PhiDotDot' * PhiDotDot + V * PhiDot' * PhiDot + ...
         D * (PhiDot - omegaTrajectory(1) * Phi)' * (PhiDot - omegaTrajectory(1)* Phi);
     G01 = D * (PhiDot - omegaTrajectory(1) * Phi)' * (omegaTrajectory(1) * Phi);
     G11 = Q * Phi' * Phi + D * omegaTrajectory(1)^2 * Phi' * Phi;
     
     g1 =  -Q * cmpTrajectory(:,1:2)' * Phi;
     g0 = zeros(size(g1));

     
     G = [G00, G01;...
          G01', G11];
     g = [g0 g1];
     
     H = [G CE'; CE zeros(size(CE,1))];
     h = [g -ce'];

     % solve problem
     x = -inv(H) * h';
     dcmKnots = x(1:numberOfKnots,:);
     copKnots = x(numberOfKnots+1:2*numberOfKnots,:);
     size(dcmKnots)
     size(copKnots)
  
  dcmTrajectory(:,1:2) = Phi * dcmKnots;
  dcmDotTrajectory(:,1:2) = PhiDot * dcmKnots;
  copTrajectory = Phi * copKnots;
  for i = 1:length(dcmTrajectory)
      vrpTrajectory(i,1:2) = dcmTrajectory(i,:) - dcmDotTrajectory(i,:) / omegaTrajectory(i);
  end
end
