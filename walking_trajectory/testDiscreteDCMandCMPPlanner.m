% Author: Robert Griffin (rjgriffin42@gmail.com)
% Date: 2/7/14

include;
clear all;

% define constants
plannerDT = 0.005;
mass = 75; % kg

% define gait parameters
doubleSupportRatio = 0.15; % double support ratio
comHeightNominal = 0.85; % nominal com height;
toeOffRatio = 0.7;
heelStrikeRatio = 0.0;

% define LQR parameters
Q = 1;
R = 1e-4;
F = 1e6;
P = 1e0;
C = 1e3;

Q_old = 1e-2;
R_old = 1e-4;
F_old = 1e6;

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
N = min(length(timeVector), length(timeVector));

% compute inertial step plan
stepPlan = transformStepPlanToInertialFrame(leftFootPoseInitial, ...
              rightFootPoseInitial, stepPlan);

footstepPlan = computeFootstepPlan(stepPlan, doubleSupportRatio, plannerDT);

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

% plan dcm trajectory
[dcmOldTrajectory, dcmDotOldTrajectory, vrpOldTrajectory] = ...
    planDiscreteDCMHybrid(copTrajectory, leftFootPoseInitial, ...
    rightFootPoseInitial, stepPlan, omegaTrajectory, omegaDotTrajectory, ...
    dcmInitial, dcmDotInitial, timeVector, Q_old, R_old, F_old);

% plan com trajectory
[comOldTrajectory, comDotOldTrajectory] = ...
    planDiscreteCoMGivenDCM(dcmOldTrajectory, omegaTrajectory, comInitial, ...
    timeVector);

% plan dynamic dcm trajectory
[dcmTrajectory, dcmDotTrajectory, vrpTrajectory] = ...
    planDiscreteDCMandCMP(copTrajectory, leftFootPoseInitial, ...
    rightFootPoseInitial, footstepPlan, omegaTrajectory, omegaDotTrajectory, ...
    dcmInitial, dcmDotInitial, timeVector, Q, R, F, P, C, N);

% plan com trajectory
[comTrajectory, comDotTrajectory] = planDiscreteCoMGivenDCM(dcmTrajectory, ...
    omegaTrajectory, comInitial, timeVector);

figure;
subplot(2,1,1)
plot(timeVector, dcmOldTrajectory(:,1), 'b', ...
     timeVector, vrpOldTrajectory(:,1), 'r', ...
     timeVector, comOldTrajectory(:,1), 'g', ...
     timeVector, dcmTrajectory(:,1), 'b--', ...
     timeVector, vrpTrajectory(:,1), 'r--', ...
     timeVector, comTrajectory(:,1), 'g--')
subplot(2,1,2)
plot(timeVector, dcmOldTrajectory(:,2), 'b', ...
     timeVector, vrpOldTrajectory(:,2), 'r', ...
     timeVector, comOldTrajectory(:,2), 'g', ...
     timeVector, dcmTrajectory(:,2), 'b--', ...
     timeVector, vrpTrajectory(:,2), 'r--', ...
     timeVector, comTrajectory(:,2), 'g--')
