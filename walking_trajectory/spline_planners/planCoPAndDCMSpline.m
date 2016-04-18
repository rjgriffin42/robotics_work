function [dcmTrajectory, dcmDotTrajectory, vrpTrajectory, cmpTrajectory] = ...
    planCoPAndDCMSpline(cmpTrajectory, leftFootPose, rightFootPose, ...
    footstepPlan, omegaTrajectory, omegaDotTrajectory,...
    dcmInitial, dcmDotInitial, numberOfKnots, plannerDT, Qp, R, D)


  footDimensions = [0.5, 0.1;...
                    0.5, -0.1;...
                   -0.5, -0.1;...
                   -0.5, 0.1];

  useDynamicsAsObjective = true;
  gravity = 9.81;
  timeVector = footstepPlan.timeVector;
  stepPlan = footstepPlan.stepPlan;
  numberOfSteps = length(stepPlan);
  knotsPerPhase = 11;

  cmpInitial = dcmInitial;
  dcmEnd = 0.5 * (stepPlan{end}.pose(1:2) + stepPlan{end-1}.pose(1:2));
  cmpEnd = cmpTrajectory(end,:);

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
  
  % compute vrp reference
  zInitial = 9.81 * ones(length(timeVector), 1) ./ (omegaTrajectory.^2 - ...
      omegaDotTrajectory);
  vrpTrajectory = cmpTrajectory + [zeros(length(timeVector), 1) ...
      zeros(length(timeVector), 1) zInitial];
  
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
  for i = 1:length(dcmTrajectory)
      vrpTrajectory(i,:) = dcmTrajectory(i,:) - dcmDotTrajectory(i,:) ./ omegaTrajectory(i);
  end

  % assemble spline coefficient matrix
  %{
  totalTime = 0;
  for step = 1:length(stepPlan)
      totalTime = footstepPlan.doubleSupportDuration(step) + totalTime;
      totalTime = footstepPlan.singleSupportDuration(step) + totalTime;
  end
  neededKnots = floor((totalTime / plannerDT) / (numberOfKnots - 1));
  numberOfKnots = (totalTime / plannerDT) / neededKnots + 1;
  segmentTime = totalTime / (numberOfKnots - 1);
  segmentDuration = 0;
  for segment = 1:numberOfKnots-1
      segmentDuration(segment) = segmentTime;
  end
  %}
  
  segmentIndex = 0;
  numberOfSegments = knotsPerPhase - 1;
  for step = 1:numberOfSteps
      phaseTime = footstepPlan.doubleSupportDuration(step)
      segmentTime = phaseTime / numberOfSegments;
      for segment = (segmentIndex+1):(segmentIndex+numberOfSegments);
          segmentDuration(segment) = segmentTime;
      end
      segmentIndex = segmentIndex + numberOfSegments;
      
      phaseTime = footstepPlan.singleSupportDuration(step)
      segmentTime = phaseTime / numberOfSegments;
      for segment = (segmentIndex+1):(segmentIndex+numberOfSegments)
          segmentDuration(segment) = segmentTime;
      end
      segmentIndex = segmentIndex + numberOfSegments;
  end
  numberOfKnots = segmentIndex + 1;
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
  size(Phi)
  size(cmpTrajectory)

  % compute vrp reference
  zInitial = 9.81 * ones(length(timeVector), 1) ./ (omegaTrajectory.^2 - ...
      omegaDotTrajectory);
  vrpTrajectory = cmpTrajectory + [zeros(length(timeVector), 1) ...
      zeros(length(timeVector), 1) zInitial];

  if (useDynamicsAsObjective)
     % add constraints on initial and final value
     CE = zeros(4, 2*numberOfKnots);
     CE(1,1) = 1;               % constrain dcm initial
     CE(2,numberOfKnots) = 1;   % constrain dcm final
     CE(3,numberOfKnots+1) = 1; % constrain cmp initial
     CE(4,end) = 1;             % constrain cmp final
     ce(1,:) = dcmInitial(1:2); % dcm initial
     ce(2,:) = dcmEnd;          % dcm final
     ce(3,:) = cmpInitial(1:2); % cmp initial
     ce(4,:) = cmpEnd(1:2);   % cmp final
     
     % assemble cost matrices
     G00 = R * PhiDotDot' * PhiDotDot;        % minimize dcm acceleration
     G00 = G00 + D * (PhiDot - omegaTrajectory(1) * Phi)' * (PhiDot - omegaTrajectory(1) * Phi); % dynamic constraint
     G11 = (Qp + omegaTrajectory(1)^2 * D) * Phi' * Phi;
     G01 = D * (omegaTrajectory(1) * Phi)' * (PhiDot - omegaTrajectory(1) * Phi);
     g0 = zeros(numberOfKnots, 2)';
     g1 = -Qp * cmpTrajectory(:,1:2)' * Phi;
     G = [G00 G01; G01' G11];
     g = [g0 g1];
  
     % restructure to include constraints
     H = [G CE'; CE zeros(size(CE,1))];
     h = [g -ce'];

     % solve problem
     x = -inv(H) * h';
     knots = x(1:2*numberOfKnots,:);
     dcmKnots = knots(1:numberOfKnots,:);
     copKnots = knots(numberOfKnots+1:2*numberOfKnots,:);
  else
     % add constraints on initial and final value
     CE = zeros(2+size(Phi,1), 2*numberOfKnots);
     CE(1,1) = 1;
     CE(2,end) = 1;
     CE(3:(2+size(Phi,1)),:) = [(PhiDot - omegaTrajectory(1) * Phi) (omegaTrajectory(1) * Phi)];
     ce(1,:) = dcmInitial(1:2);
     ce(2,:) = dcmEnd;
     ce(3:(2+size(Phi,1)),:) = 0;
     size(CE)
     size(ce)
     
     % assemble cost matrices
     G00 = R * PhiDotDot' * PhiDotDot;        % minimize dcm acceleration
     G11 = Qp * Phi' * Phi;
     G01 = zeros(numberOfKnots, numberOfKnots);
     g0 = zeros(numberOfKnots, 2)';
     g1 = -Qp * cmpTrajectory(:,1:2)' * Phi;
     G = [G00 G01; G01' G11];
     g = [g0 g1];
  
     % restructure to include constraints
     H = [G CE'; CE zeros(size(CE,1))];
     h = [g -ce'];

     % solve problem
     x = -inv(H) * h';
     knots = x(1:2*numberOfKnots,:);
     dcmKnots = knots(1:numberOfKnots,:);
     copKnots = knots(numberOfKnots+1:2*numberOfKnots,:);
  end
  
  dcmTrajectory(:,1:2) = Phi * dcmKnots;
  dcmDotTrajectory(:,1:2) = PhiDot * dcmKnots;
  cmpTrajectory(:,1:2) = Phi * copKnots;
 
  for i = 1:length(dcmTrajectory)
      vrpTrajectory(i,1:3) = dcmTrajectory(i,:) - dcmDotTrajectory(i,:) / omegaTrajectory(i);
  end
end