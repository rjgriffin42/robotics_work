%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based off of the work presented in
%    Michael A Hopkins, Dennis Hong, and Alex Leonessa, "Humanoid Locomotion on
%      Uneven Terrain Using the Time-Varying Divergent Component of Motion."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dcmTrajectory, dcmDotTrajectory, vrpTrajectory] = ...
    planDCMHybrid(cmpTrajectory, leftFootPose, rightFootPose, stepPlan, ...
    omegaTrajectory, omegaDotTrajectory, dcmInitial, dcmDotInitial, timeVector,...
    Q, R, F)

  gravity = 9.81;

  % compute vrp reference
  zInitial = 9.81 * ones(length(timeVector), 1) ./ (omegaTrajectory.^2 - ...
      omegaDotTrajectory);
  vrpTrajectory = cmpTrajectory + [zeros(length(timeVector), 1) ...
      zeros(length(timeVector), 1) zInitial];

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

  % initialize state matrices
  N = min(200, length(timeVector));
  timeVector = timeVector';
  dt = timeVector(2) - timeVector(1);
  A = [1, dt; 0 1];                       % state = dcm position and velocity
  B = [dt^2/2; dt];                       % input = dcm acceleration
  CC = zeros(N, 2*N);                     % output = vrp
  for i = 1:N
    CC(i, 2*i-1:2*i) = [1, -omegaTrajectory(i) / (omegaTrajectory(i)^2 - ...
        omegaDotTrajectory(i))];
  end

  % compute transition matrices
  X0 = zeros(2*N, 2);             % output transition matrix for t in (1, tf)
  Xu = zeros(2*N, N);             % output transition matrix for t in (1, tf)
  Y0 = zeros(N, 2);
  Yu = zeros(N, N);
  Ak = eye(size(A));
  for i = 1:N
    Ak_minus_B = Ak * B;
    Ak = A * Ak;
    X0((2*i-1):2*i, :) = Ak;
    for j = i:N
      Xu(2*j-1:2*j, (j-i)+1) = Ak_minus_B;
    end
  end
  Y0 = CC * X0;
  Yu = CC * Xu;

  % compute optimal control trajectory using QP differential
  xInitial = [dcmInitial; dcmDotInitial];
  xFinal = [dcmTrajectory(N, :); dcmDotTrajectory(N, :)];
  Asolve = Q*Yu'*Yu + F*Xu(end-1:end, :)'*Xu(end-1:end, :) + R*eye(N,N);
  for i = 1:2
    bsolve = -(Q*Yu'*(Y0*xInitial(:, i) - vrpTrajectory(1:N, i)) + ...
        F*Xu(end-1:end, :)' * (X0(end-1:end, :) * xInitial(:, i) - xFinal(:, i)));
    U(:, i) = Asolve \ bsolve;
  end

  % compute dcm and dcmdot trajectories
  for i = 1:2
    X = Xu*U(:,i) + X0*xInitial(:, i);
    Y = Yu*U(:,i) + Y0*xInitial(:, i);
    for j = 1:N
      dcmTrajectory(j, i) = X(2*j-1);
      dcmDotTrajectory(j, i) = X(2*j);
      vrpTrajectory(j, i) = Y(j);
    end
  end
end
