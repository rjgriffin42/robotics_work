%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implements the paper
% Russ Kuindersma, Scott Kuindersma, Robin Deits, and Kanako Miura, "A closed-form
%    form solution for real-time ZMP gait generation and feedback stabilization"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [zmpTrajectory] = planLQRZMPTrajectories(desiredZMPTrajectory, nominalCoMHeight, Q, R)
  gravity = -9.81;

  % define ZMP state-space dynamics
  A = [zeros(2,2), eye(2,2); zeros(2,4)];
  B = [zeros(2,2); eye(2,2)];
  C = [eye(2,2), zeros(2,2)];
  D = -nominalCoMHeight / gravity * eye(2,2);

  ybar = desiredZMPTrajectory - desiredZMPTrajectory(end);

  Q1 = C' * Q * C;
  q2 = -2 * C' * Q * ybar;
  q3 = Q * ybar' * ybar;
  R1 = R + Q * D' * D;
  r2 = -2 * D * Q * ybar;
  N = C' * Q * D;
end