function [copTrajectory] = planCoPSpline(cmpTrajectory, leftFootPose,...
    rightFootPose, footstepPlan, numberOfKnots, plannerDT, Q, R, V)

  gravity = 9.81;
  timeVector = footstepPlan.timeVector;

  stepPlan = footstepPlan.stepPlan;
  
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
      rowIndex = rowIndex + size(T,1);
  end

  G = Q * Phi' * Phi;
  g = -Q * cmpTrajectory(:,1:2)' * Phi;

  % solve problem
  x = -inv(G) * g';
  knots = x(1:numberOfKnots,:);
  
  copTrajectory(:,1:2) = Phi * knots;
end