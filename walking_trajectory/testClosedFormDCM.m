% Author: Robert Griffin (rjgriffin42@gmail.com)
% Date: 2/7/14

include;
clear all;

% define constants
plannerDT = 0.005;
mass = 75; % kg

% define gait parameters
doubleSupportRatio = 0.2; % double support ratio
comHeightNominal = 0.85; % nominal com height;
toeOffRatio = 0.7;
heelStrikeRatio = 0.0;

% define LQR parameters
Q = 1e-2;
R = 1e-4;
F = 1e6;

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

% compute discrete time vector
timeVector = computeStepPlanTimeVector(stepPlan, plannerDT);

% compute inertial step plan
stepPlan = transformStepPlanToInertialFrame(leftFootPoseInitial, ...
              rightFootPoseInitial, stepPlan);
footstepPlan = computeFootstepPlan(stepPlan, doubleSupportRatio, plannerDT);

% plan com height trajectory
[comHeightTrajectory, comDotHeightTrajectory, comDotDotHeightTrajectory] = ...
    planDiscreteCoMFlatHeightTrajectory(leftFootPoseInitial, ...
    rightFootPoseInitial, stepPlan, comHeightNominal, doubleSupportRatio, ...
    comInitial, comDotInitial, comDotDotInitial, timeVector);

% plan dcm trajectory
[dcmTrajectory, vrpTrajectory] = ...
    planClosedFormDCM(footstepPlan, doubleSupportRatio, omegaInitial, timeVector);

% plan com trajectory
%[comTrajectory, comDotTrajectory] = planDiscreteCoMGivenDCM(dcmTrajectory, ...
%    omegaTrajectory, comInitial, timeVector);

figure;
subplot(3,1,1)
plot(timeVector, dcmTrajectory(:,1), timeVector, vrpTrajectory(:,1))
subplot(3,1,2)
plot(timeVector, dcmTrajectory(:,2), timeVector, vrpTrajectory(:,2))
subplot(3,1,3)
plot(dcmTrajectory(:,1), dcmTrajectory(:,2), vrpTrajectory(:,1), vrpTrajectory(:,2))
