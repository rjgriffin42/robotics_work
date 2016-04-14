% Author: Robert Griffin (rjgriffin42@gmail.com)
% Date: 2/7/14

include;
clear all;

% define constants
plannerDT = 0.005;
mass = 75; % kg

% define gait parameters
doubleSupportRatio = 0.2; % double support ratio
initialDoubleSupportDuration = 0.2;
comHeightNominal = 0.85; % nominal com height;
toeOffRatio = 0.7;
heelStrikeRatio = 0.0;

% define LQR parameters
Q = 1e-2;
R = 1e-4;
F = 1e6;
numberOfKnots = 100;

% define initial conditions
leftFootPoseInitial = [0, 0.1, 0, 0, 0, 0];
rightFootPoseInitial = [0,-0.1, 0, 0, 0, 0];
comInitial = [0, 0, 0.85];
comDotInitial = [0, 0, 0];
comDotDotInitial = [0, 0, 0];
copInitial = [comDotInitial(1) comDotInitial(2) ...
              (leftFootPoseInitial(3) + rightFootPoseInitial(3))];
omegaInitial = 1/sqrt(9.81/comDotInitial(3));
dcmInitial = comDotInitial + omegaInitial * comDotInitial;
dcmDotInitial = zeros(1,3);

% define step plan
stepPlan = forwardStepPlan();

% compute discrete time vector
timeVector = computeStepPlanTimeVector(stepPlan, plannerDT);

% compute inertial step plan
stepPlan = transformStepPlanToInertialFrame(leftFootPoseInitial, ...
              rightFootPoseInitial, stepPlan);
          
% compute footstep plan
footstepPlan = computeFootstepPlan(stepPlan, doubleSupportRatio, plannerDT, initialDoubleSupportDuration);

% plan cop trajectory
copTrajectory = planDiscreteCOPToeOff(leftFootPoseInitial, rightFootPoseInitial, ...
    stepPlan, doubleSupportRatio, toeOffRatio, heelStrikeRatio, copInitial, ...
    timeVector);

% plan com height trajectory
[comHeightTrajectory, comDotHeightTrajectory, comDotDotHeightTrajectory] = ...
    planDiscreteCoMFlatHeightTrajectory(leftFootPoseInitial, ...
    rightFootPoseInitial, stepPlan, comHeightNominal, doubleSupportRatio, ...
    comInitial, comDotInitial, comDotDotInitial, timeVector);

% plan angular momentum rate of change
tauTrajectory = zeros(size(copTrajectory));
for i = 1:length(timeVector)
  tauTrajectory(i, 1) = 0;
  tauTrajectory(i, 2) = 0;
end

% plan cmp trajectory
cmpTrajectory = planDiscreteCMPTrajectory(copTrajectory, tauTrajectory, ...
    comDotDotHeightTrajectory, mass, timeVector);

% plan omega trajectory
[omegaTrajectory, omegaDotTrajectory] = ...
    planDiscreteOmegaTrajectory(cmpTrajectory, comHeightTrajectory, ...
    comDotHeightTrajectory, comDotDotHeightTrajectory, timeVector);

% plan vrp trajectory
vrpTrajectory = planDiscreteVRPTrajectory(cmpTrajectory, omegaTrajectory, ...
    omegaDotTrajectory, timeVector);

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

% plan dcm trajectory
[dcmTrajectory, dcmDotTrajectory, vrpTrajectory, cmpTrajectoryNew] = ...
    planCoPAndDCMSpline(cmpTrajectory, leftFootPoseInitial, ...
    rightFootPoseInitial, footstepPlan, omegaTrajectory, omegaDotTrajectory, ...
    dcmInitial, dcmDotInitial, numberOfKnots, plannerDT, 1e10, R, 1e3);

[dcmTrajectoryOld, dcmDotTrajectoryOld, vrpTrajectory] = ...
    planDiscreteDCMHybrid(cmpTrajectory, leftFootPoseInitial, ...
    rightFootPoseInitial, stepPlan, omegaTrajectory, omegaDotTrajectory, ...
    dcmInitial, dcmDotInitial, timeVector, Q, R, F);

figure;
subplot(3,1,1)
plot(timeVector, dcmTrajectory(:,1), timeVector, dcmTrajectoryOld(:,1), '--',...
    timeVector, cmpTrajectoryNew(:,1), '-.', timeVector, vrpTrajectory(:,1), ':');
subplot(3,1,2)
plot(timeVector, dcmTrajectory(:,2), timeVector, dcmTrajectoryOld(:,2), '--',...
    timeVector, cmpTrajectoryNew(:,2), '-.', timeVector, vrpTrajectory(:,2), ':');
       
subplot(3,1,3)
plot(dcmTrajectory(:,1), dcmTrajectory(:,2), ...
     dcmTrajectoryOld(:,1), dcmTrajectoryOld(:,2), '--')
 
figure;
subplot(2,1,1)
plot(timeVector, cmpTrajectoryNew(:,1), timeVector, cmpTrajectory(:,1), '--')
subplot(2,1,2)
plot(timeVector, cmpTrajectoryNew(:,2), timeVector, cmpTrajectory(:,2), '--')

% plan com trajectory
[comTrajectory, comDotTrajectory] = planDiscreteCoMGivenDCM(dcmTrajectoryOld, ...
    omegaTrajectory, comInitial, timeVector);
