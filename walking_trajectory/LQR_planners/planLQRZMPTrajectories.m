%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implements the paper
% Russ Kuindersma, Scott Kuindersma, Robin Deits, and Kanako Miura, "A closed-form
%    form solution for real-time ZMP gait generation and feedback stabilization"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [S1] = planLQRZMPTrajectories(footstepPlan, nominalCoMHeight, zmp0, Q, R)
  
  stepPlan = footstepPlan.stepPlan;
  numberOfSteps = length(stepPlan);

  zmpInitial{1} = zmp0;
  zmpFinal{1} = [0 0.1];

  for stepIndex = 2:numberOfSteps
    zmpInitial{stepIndex} = stepPlan{stepIndex-1}.pose(1:2);
    zmpFinal{stepIndex} = stepPlan{stepIndex}.pose(1:2);
  end

  % compute spline coefficients
  for stepIndex = 1:numberOfSteps
  	xInitial = [zmpInitial{stepIndex}' [0; 0]];
    xFinal = [zmpFinal{stepIndex}' [0; 0]];

  	c{2*stepIndex - 1} = computeContinuousCubicCoefficients(xInitial, xFinal, ...
  		footstepPlan.doubleSupportDuration(stepIndex));

  	c{2*stepIndex - 1} = computeContinuousCubicCoefficients(xFinal, xFinal, ...
  		footstepPlan.singleSupportDuration(stepIndex));
  end

  gravity = -9.81;

  % define ZMP state-space dynamics
  A = [zeros(2,2), eye(2,2); zeros(2,4)];
  B = [zeros(2,2); eye(2,2)];
  C = [eye(2,2), zeros(2,2)];
  D = -nominalCoMHeight / gravity * eye(2,2);

  %ybar = desiredZMPTrajectory - desiredZMPTrajectory(end);

  Q1 = C' * Q * C;
  %q2 = -2 * C' * Q * ybar;
  %q3 = Q * ybar' * ybar;
  R1 = R + Q * D' * D;
  %r2 = -2 * D * Q * ybar;
  N = C' * Q * D;

  % find the solution to the algebraic ricatti equation
  S1 = care(A, B, Q1, R1, N);
  S11 = S1(1:2, 1:2);
  S12 = S1(1:2, 3:4);
  S22 = S1(3:4, 3:4);

  % compute linear time varying LQR term
  NB = N' + B' * S1;
  B2 = 2 * (C' - NB' * R1^-1 * D) * Q;

  s2 = computeLQRTimeVaryingLinearTerm(S11, S12, S22, B2, Q, D, R1, c);

  zmpTrajectory = 0;
end