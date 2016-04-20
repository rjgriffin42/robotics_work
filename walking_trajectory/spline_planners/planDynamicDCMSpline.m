function [dcmTrajectory, dcmDotTrajectory, vrpTrajectory, copTrajectory] = ...
    planDynamicDCMSpline(cmpTrajectory, leftFootPose, rightFootPose, footstepPlan, ...
    omegaTrajectory, omegaDotTrajectory, dcmInitial, dcmDotInitial,...
    numberOfKnots, plannerDT, Q, R, V, D, Q_goal)

  gravity = 9.81;
  timeVector = footstepPlan.timeVector;
  stepPlan = footstepPlan.stepPlan;
  numberOfSteps = length(stepPlan);
 
  
  knotsPerPhase = 21;
  alpha = 0.5;

  % compute vrp reference
  zInitial = 9.81 * ones(length(timeVector), 1) ./ (omegaTrajectory.^2 - ...
      omegaDotTrajectory);
  vrpTrajectory = cmpTrajectory + [zeros(length(timeVector), 1) ...
      zeros(length(timeVector), 1) zInitial];

 
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
  segmentsPerPhase = knotsPerPhase - 1;
  numberOfKnots = 2 * segmentsPerPhase * numberOfSteps + 1;
  segmentIndex = 0;
  for step = 1:numberOfSteps
      segmentTime = footstepPlan.doubleSupportDuration(step) / segmentsPerPhase;
      for segment = segmentIndex+1:segmentIndex+segmentsPerPhase
          segmentDuration(segment) = segmentTime;
      end
      segmentIndex = segmentIndex + segmentsPerPhase;
      
      segmentTime = footstepPlan.singleSupportDuration(step) / segmentsPerPhase;
      for segment = segmentIndex+1:segmentIndex+segmentsPerPhase
          segmentDuration(segment) = segmentTime;
      end
      segmentIndex = segmentIndex + segmentsPerPhase;
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
  
  planLength = length(stepPlan);
  for i = 1:planLength
    vrpKnots{i+1} = stepPlan{i}.pose(1:3);
  end
  
  vrpKnots{1} = [0, 0.1, 0];
  % add in first one for the initial foot position
  dcmFinal = 0.5 * (stepPlan{numberOfSteps}.pose(1:2) + stepPlan{numberOfSteps-1}.pose(1:2));
  dcmInitialGoal{numberOfSteps+1} = dcmFinal;
  vrpKnots{numberOfSteps+1} = dcmFinal;
  
  % compute knot points for DCM trajectory at heel strike
  for i = numberOfSteps:-1:1
    duration{i} = stepPlan{i}.duration;
    dcmEndGoal{i} = dcmInitialGoal{i+1};

    doubleSupportDuration{i} = footstepPlan.doubleSupportDuration(i);
    singleSupportDuration{i} = footstepPlan.singleSupportDuration(i);

    alpha = exp(omegaTrajectory(1) * duration{i});
    
    dcmInitialGoal{i} = (dcmEndGoal{i} - vrpKnots{i}(1:2)) / alpha + vrpKnots{i}(1:2);
  end

     
     
 
     
     % add constraints on initial and final value
     CE = zeros(2, 2 * numberOfKnots);
     CE(1,1) = 1;
     CE(2,numberOfKnots) = 1;
     ce(1,:) = dcmInitial(1:2);
     ce(2,:) = dcmFinal;
     
     
     % add constraints on DCM value at heel strike
     S = zeros(numberOfSteps, numberOfKnots);
     for i = 1:numberOfSteps
         index = 2 * i * segmentsPerPhase + 1;
         S(i, index) = 1;
         dcmGoal(i,:) = dcmInitialGoal{i+1};
         %CE(i+2,(i * 2 * (knotsPerPhase-1))+1) = 1;
         %ce(i+2,:) = dcmInitialGoal{i + 1};
     end
     
     G00 = R * PhiDotDot' * PhiDotDot + V * PhiDot' * PhiDot + ...
         D * (PhiDot - omegaTrajectory(1) * Phi)' * (PhiDot - omegaTrajectory(1)* Phi) ...
         + Q_goal * S' * S;
     G01 = D * (PhiDot - omegaTrajectory(1) * Phi)' * (omegaTrajectory(1) * Phi);
     G11 = Q * Phi' * Phi + D * omegaTrajectory(1)^2 * Phi' * Phi;
     
     g1 =  -Q * cmpTrajectory(:,1:2)' * Phi;
     g0 = -Q_goal * dcmGoal' * S;
     
     G = [G00, G01;...
          G01', G11];
     g = [g0 g1];
     
     H = [G CE'; CE zeros(size(CE,1))];
     h = [g -ce'];

     % solve problem
     x = -inv(H) * h';
     dcmKnots = x(1:numberOfKnots,:);
     copKnots = x(numberOfKnots+1:2*numberOfKnots,:);
  
  dcmTrajectory(:,1:2) = Phi * dcmKnots;
  dcmDotTrajectory(:,1:2) = PhiDot * dcmKnots;
  copTrajectory = Phi * copKnots;
  for i = 1:length(dcmTrajectory)
      vrpTrajectory(i,1:2) = dcmTrajectory(i,:) - dcmDotTrajectory(i,:) / omegaTrajectory(i);
  end
end
