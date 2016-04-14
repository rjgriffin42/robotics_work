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
initialDoubleSupportDuration = 0.5;
comHeightNominal = 0.85; % nominal com height;
toeOffRatio = 0.7;
heelStrikeRatio = 0.0;

% define LQR parameters
Q = 1e2;
R = 1e-1;

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
zmpInitial = dcmInitial(1:2);
dcmDotInitial = zeros(1,3);

% define step plan
stepPlan = forwardStepPlan();

% compute inertial step plan
stepPlan = transformStepPlanToInertialFrame(leftFootPoseInitial, ...
              rightFootPoseInitial, stepPlan);
footstepPlan = computeFootstepPlan(stepPlan, doubleSupportRatio, plannerDT, initialDoubleSupportDuration);
timeVector = footstepPlan.timeVector;

[zmpTrajectory, zmpDefined] = ...
    planLQRZMPTrajectories(footstepPlan, comHeightNominal, zmpInitial, Q, R);

%{
subplot(2,1,1)
plot(zmpTrajectory(1,:), zmpTrajectory(2,:))
subplot(2,1,2)
plot(zmpDefined(1,:), zmpDefined(2,:))
%}

