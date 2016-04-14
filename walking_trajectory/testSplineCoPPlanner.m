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
numberOfKnots = 50;

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
copTrajectoryOld = planDiscreteCOPToeOff(leftFootPoseInitial, rightFootPoseInitial, ...
    stepPlan, doubleSupportRatio, toeOffRatio, heelStrikeRatio, copInitial, ...
    timeVector);

% plan new cop trajectory
[copTrajectory] = ...
    planCoPSpline(copTrajectoryOld, leftFootPoseInitial, ...
    rightFootPoseInitial, footstepPlan, numberOfKnots, plannerDT, 1e-1, ...
    R, 1e-5);

subplot(2,1,1)
plot(timeVector, copTrajectory(:,1), timeVector, copTrajectoryOld(:,1), '--');%, ...
    %timeVector, dcmDotTrajectory(:,1), timeVector, dcmDotTrajectoryOld(:,1), '--')
subplot(2,1,2)
plot(timeVector, copTrajectory(:,2), timeVector, copTrajectoryOld(:,2), '--');%, ...
        %timeVector, dcmDotTrajectory(:,2), timeVector, dcmDotTrajectoryOld(:,2), '--')

