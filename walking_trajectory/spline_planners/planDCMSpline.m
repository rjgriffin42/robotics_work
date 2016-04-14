function [dcmTrajectory, dcmDotTrajectory, vrpTrajectory] = ...
    planDCMSpline(cmpTrajectory, leftFootPose, rightFootPose, footstepPlan, ...
    omegaTrajectory, omegaDotTrajectory, dcmInitial, dcmDotInitial,...
    numberOfKnots, plannerDT, Q, R, V)

  planAgainstDCM = false;
  useDynamicConstraint = false;
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

  if (planAgainstDCM)
     % integrate DCM in reverse time
      dcmTrajectory = zeros(length(timeVector), 3);
      dcmDotTrajectory = zeros(length(timeVector), 3);
      dcmTrajectory(end, :) = [singleSupportPoses{end}(1) singleSupportPoses{end}(2) ...
          vrpTrajectory(end, 3)];
      dcmDotTrajectory(end, :) = [0 0 0];
      for i = (length(timeVector)-1):-1:1
        dt = timeVector(i) - timeVector(i+1);
       k1 = (omegaTrajectory(i+1) - omegaDotTrajectory(i+1) / omegaTrajectory(i+1)) ...
            * (dcmTrajectory(i+1, :) - vrpTrajectory(i+1, :));
       k2 = (omegaTrajectory(i) - omegaDotTrajectory(i) / omegaTrajectory(i)) * ...
           (dcmTrajectory(i+1, :) + dt * k1 - vrpTrajectory(i, :));
       dcmDotTrajectory(i, :) = k1 / 2 + k2 / 2;
       dcmTrajectory(i, :) = dcmTrajectory(i+1, :) + dt * dcmDotTrajectory(i, :);
     end
  end
  
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
  
  if (planAgainstDCM)
     % add constraints on initial and final value
     CE = zeros(2, numberOfKnots);
     CE(1,1) = 1;
     CE(2,end) = 1;
     ce(1,:) = dcmInitial(1:2);
     ce(2,:) = dcmTrajectory(end,1:2);

     G = Q * Phi' * Phi + R * PhiDotDot' * PhiDotDot + V * PhiDot' * PhiDot;
     g = -Q * dcmTrajectory(:,1:2)' * Phi;

     H = [G CE'; CE zeros(size(CE,1))];
     h = [g -ce'];

     % solve problem
     x = -inv(H) * h';
     knots = x(1:numberOfKnots,:);
  elseif(useDynamicConstraint)
     % add constraints on initial and final value
     CE = zeros(2 + size(Phi,1), numberOfKnots);
     CE(1,1) = 1;
     CE(2,end) = 1;
     ce(1,:) = dcmInitial(1:2);
     ce(2,:) = [singleSupportPoses{end}(1) singleSupportPoses{end}(2)];
     
     % add constraints on dynamics
     CE(3:(2+size(Phi,1)),:) = [(PhiDot - omegaTrajectory(1) * Phi)];
     ce(3:(2+size(Phi,1)),:) = omegaTrajectory(1) * cmpTrajectory(:,1:2);
     
     G = R * PhiDotDot' * PhiDotDot + V * PhiDot' * PhiDot;
     g = zeros(2, size(Phi, 2));
     
     H = [G CE'; CE zeros(size(CE,1))];
     size(g)
     size(ce)
     h = [g -ce'];

     % solve problem
     x = -inv(H) * h';
     knots = x(1:numberOfKnots,:)
  else
     % add constraints on initial and final value
     CE = zeros(2, numberOfKnots);
     CE(1,1) = 1;
     CE(2,end) = 1;
     ce(1,:) = dcmInitial(1:2);
     ce(2,:) = [singleSupportPoses{end}(1) singleSupportPoses{end}(2)];
     
     G = R * PhiDotDot' * PhiDotDot + V * PhiDot' * PhiDot + ...
         Q * (PhiDot - omegaTrajectory(1) * Phi)' * (PhiDot - omegaTrajectory(1)* Phi);
     g = omegaTrajectory(1) * Q * cmpTrajectory(:,1:2)' * (PhiDot - omegaTrajectory(1) * Phi);
     
     H = [G CE'; CE zeros(size(CE,1))];
     size(g)
     size(ce)
     h = [g -ce'];

     % solve problem
     x = -inv(H) * h';
     knots = x(1:numberOfKnots,:)
  end
  
  dcmTrajectory(:,1:2) = Phi * knots;
  dcmDotTrajectory(:,1:2) = PhiDot * knots;
  size(dcmTrajectory)
  size(dcmDotTrajectory)
  for i = 1:length(dcmTrajectory)
      vrpTrajectory(i,1:2) = dcmTrajectory(i,:) - dcmDotTrajectory(i,:) / omegaTrajectory(i);
  end
end
