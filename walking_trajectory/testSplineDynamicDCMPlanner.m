% Author: Robert Griffin (rjgriffin42@gmail.com)
% Date: 2/7/14

include;
clear all;

% define constants
plannerDT = 0.005;
mass = 75; % kg
plan_alpha = 0.5;

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

Qspline = 1e5;   % cmp tracking
Rspline = 1e-2;  % dcm acceleration
Vspline = 0;     % dcm velocity
Dspline = 1e3;   % dcm dynamics
numberOfKnots = 100;

% define initial conditions
leftFootPoseInitial = [0, 0.1, 0, 0, 0, 0];
rightFootPoseInitial = [0,-0.1, 0, 0, 0, 0];
comInitial = [0, 0, 0.85];
comDotInitial = [0, 0, 0];
comDotDotInitial = [0, 0, 0];
copInitial = [comDotInitial(1) comDotInitial(2) ...
              (leftFootPoseInitial(3) + rightFootPoseInitial(3))];
omegaInitial = sqrt(9.81/comInitial(3));
dcmInitial = comDotInitial + 1 / omegaInitial * comDotInitial;
dcmDotInitial = zeros(1,3);

% define step plan
stepPlan = forwardStepPlan();

% compute inertial step plan
stepPlan = transformStepPlanToInertialFrame(leftFootPoseInitial, ...
              rightFootPoseInitial, stepPlan);
footstepPlan = computeFootstepPlan(stepPlan, doubleSupportRatio, plannerDT, initialDoubleSupportDuration);
timeVector = footstepPlan.timeVector;

% plan com height trajectory
[comHeightTrajectory, comDotHeightTrajectory, comDotDotHeightTrajectory] = ...
    planDiscreteCoMFlatHeightTrajectory(leftFootPoseInitial, ...
    rightFootPoseInitial, stepPlan, comHeightNominal, doubleSupportRatio, ...
    comInitial, comDotInitial, comDotDotInitial, timeVector);

% plan cop trajectory
copTrajectory = planDiscreteCOPToeOff(leftFootPoseInitial, rightFootPoseInitial, ...
    stepPlan, doubleSupportRatio, toeOffRatio, heelStrikeRatio, copInitial, ...
    timeVector);

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

% plan dcm trajectory with splines
[dcmTrajectorySpline, dcmDotTrajectorySpline, vrpTrajectorySpline, copTrajectorySpline] = ...
    planDynamicDCMSpline(cmpTrajectory, leftFootPoseInitial, ...
    rightFootPoseInitial, footstepPlan, omegaTrajectory, omegaDotTrajectory, ...
    dcmInitial, dcmDotInitial, numberOfKnots, plannerDT, Qspline, Rspline, Vspline, Dspline);

% plan continuous double support dcm trajectory
[dcmTrajectoryCDS, dcmDotTrajectoryCDS, vrpTrajectoryCDS] = ...
    planCDSClosedFormDCM(footstepPlan, doubleSupportRatio, plan_alpha,...
    dcmInitial, omegaInitial, plannerDT, initialDoubleSupportDuration);

% plan using the discrete method
[dcmTrajectoryDiscrete, dcmDotTrajectoryDiscrete, vrpTrajectoryDiscrete] = ...
    planDiscreteDCMHybrid(cmpTrajectory, leftFootPoseInitial, ...
    rightFootPoseInitial, stepPlan, omegaTrajectory, omegaDotTrajectory, ...
    dcmInitial, dcmDotInitial, timeVector, Q, R, F);

figure;
subplot(3,1,1)
plot(timeVector, dcmTrajectorySpline(:,1), timeVector, dcmTrajectoryCDS(:,1), '--',...
    timeVector, dcmTrajectoryDiscrete(:,1), '-.');

subplot(3,1,2)
plot(timeVector, dcmTrajectorySpline(:,2), timeVector, dcmTrajectoryCDS(:,2), '--', ...
    timeVector, dcmTrajectoryDiscrete(:,2), '-.');

subplot(3,1,3)
plot(dcmTrajectorySpline(:,1), dcmTrajectorySpline(:,2), vrpTrajectory(:,1), vrpTrajectory(:,2))

figure;
subplot(3,1,1)
plot(timeVector, copTrajectorySpline(:,1), timeVector, cmpTrajectory(:,1), '--');

subplot(3,1,2)
plot(timeVector, copTrajectorySpline(:,2), timeVector, cmpTrajectory(:,2), '--');