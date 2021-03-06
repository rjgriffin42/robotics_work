%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implements the paper
% Russ Kuindersma, Scott Kuindersma, Robin Deits, and Kanako Miura, "A closed-form
%    form solution for real-time ZMP gait generation and feedback stabilization"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [zmpTrajectory, zmpDefined] = planLQRZMPTrajectories(footstepPlan, nominalCoMHeight, zmp0, Q, R)
  
  plannerDT = 0.005;
  stepPlan = footstepPlan.stepPlan;
  numberOfSteps = length(stepPlan);

  zmpInitial{1} = zmp0;
  zmpFinal{1} = [0 0.1];


  for stepIndex = 2:numberOfSteps
    zmpInitial{stepIndex} = stepPlan{stepIndex-1}.pose(1:2);
    zmpFinal{stepIndex} = stepPlan{stepIndex}.pose(1:2);
  end

  timeKnots{1} = 0;
  % compute spline coefficients
  for stepIndex = 1:numberOfSteps
  	xKnots{1} = [zmpInitial{stepIndex}' [0; 0]];
    xKnots{2} = [zmpFinal{stepIndex}' [0; 0]];
    timeKnots{2*stepIndex} = footstepPlan.doubleSupportDuration(stepIndex)...
        + timeKnots{2*stepIndex-1};
    timeKnots{2*stepIndex+1} = footstepPlan.singleSupportDuration(stepIndex)...
        + timeKnots{2*stepIndex};

  	coefficients_temp = computeContinuousCubicCoefficients(xKnots, ...
  		footstepPlan.doubleSupportDuration(stepIndex));
    coefficients{2*stepIndex - 1} = coefficients_temp{1};

  	coefficients_temp = computeContinuousCubicCoefficients(xKnots, ...
  		footstepPlan.singleSupportDuration(stepIndex));
    coefficients{2*stepIndex} = coefficients_temp{1};
  end
  
  % compute cubic spline ZMP trajectory
  lineIndex = 0;
  for segment = 1:length(coefficients)
      t = timeKnots{segment}:plannerDT:(timeKnots{segment+1} - plannerDT);
      for col = 1:length(t)
          timeVector(col+lineIndex) = t(col);
          zmp(:,col+lineIndex) = [0;0];
          for index = 1:4
              zmp(:,col+lineIndex) = zmp(:,col+lineIndex) + ...
                  coefficients{segment}(:,index) * (t(col) - timeKnots{segment})^(index-1);
          end
      end
      lineIndex = lineIndex + length(t);
  end
  zmp_bar(1,:) = zmp(1,:) - zmpFinal{end}(1);
  zmp_bar(2,:) = zmp(2,:) - zmpFinal{end}(2);
  zmpDefined = zmp;

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
  r2 = -2 * D * Q * zmp_bar;
  N = C' * Q * D;
  R1inv = inv(R1);

  % find the solution to the algebraic ricatti equation
  S1 = care(A, B, Q1, R1, N);
  S11 = S1(1:2, 1:2);
  S12 = S1(1:2, 3:4);
  S22 = S1(3:4, 3:4);
  
  [K, S] = lqr(A,B,Q1,R1,N);

  K = -K;
  
  % compute linear time varying LQR term
  NB = N' + B' * S;
  A2 = NB' * R1inv * B' - A';
  B2 = 2 * (C' - NB' * R1inv * D) * Q;

  A2inv = inv(A2);
  
  n = numberOfSteps * 2;
  k = 3;
  
  [s2, k2, alpha, beta] = computeLQRTimeVaryingLinearTerm(S11, S12, S22, B2, Q, B, D, R1,...
     coefficients, timeKnots, plannerDT);
  

  
  Ay = [A + B*K, -.5 * B * R1inv * B'; zeros(4), A2];
  Ayinv = inv(Ay);
  By = [B*R1inv * D * Q; B2];
  
  a = zeros(4, n);
  b = zeros(4,n,k);
  
  for j = 1:n
      b(:,j,k) = -Ayinv(1:4,:) * By * coefficients{j}(:,k);
      for i = k-1:-1:1
          b(:,j,i) = Ayinv(1:4,:) * (i*[b(:,j,i+1); beta(:,j,i+1)] - ...
              By * coefficients{j}(:,i));
      end
      a(:,j) = x - b(:,j,1);
      dt = timeKnots{j+1} - timeKnots{j};
      x = [eye(4), zeros(4)] * expm(Ay * dt) * [a(:,j); alpha(:,j)] + ...
          squeeze(b(:,j,:)) * (dt.^(0:k-1)');
  end
  
  rs = 1/2 * (r2 + B' * s2);
  
  u = -inv(R1) * (NB * x + rs);
  zmpTrajectory = C * x + D * u;
end
