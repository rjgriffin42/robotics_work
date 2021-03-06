function [dcmTrajectory, dcmDotTrajectory, vrpTrajectory] = ...
    planDCMandCMP(copTrajectory, leftFootPose, rightFootPose, footstepPlan, ...
    omegaTrajectory, omegaDotTrajectory, dcmInitial, dcmDotInitial, timeVector,...
    Q, R, F, P, C, N)

  stepPlan = footstepPlan.stepPlan;
  gravity = 9.81;

  % compute vrp reference. Assuming CMP and CoP equivalence for now.
  zInitial = 9.81 * ones(length(timeVector), 1) ./ (omegaTrajectory.^2 - ...
      omegaDotTrajectory);
  vrpTrajectory = copTrajectory + [zeros(length(timeVector), 1) ...
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
  timeVector = timeVector';
  dt = timeVector(2) - timeVector(1);
  A = [1, dt; 0 1];                       % state = dcm position and velocity
  B = [dt^2/2; dt];                       % input = dcm acceleration
  Ck = zeros(N, 2*N);                     % output = vrp
  for i = 1:N
    Ck(i, 2*i-1:2*i) = [1, -omegaTrajectory(i) / (omegaTrajectory(i)^2 - ...
        omegaDotTrajectory(i))];
  end

  % compute transition matrices
  PhiX0 = zeros(2*N, 2);             % dcm transition matrix for t in (1, tf)
  PhiXu = zeros(2*N, N);             % dcm transition matrix for t in (1, tf)
  PhiY0 = zeros(N, 2);               % output transition matrix for t in (1, tf)
  PhiYu = zeros(N, N);               % output transition matrix for t in (1, tf)
  Ak = eye(size(A));
  for i = 1:N
    Ak_minus_B = Ak * B;
    Ak = A * Ak;
    PhiX0((2*i-1):2*i, :) = Ak;
    for j = i:N
      PhiXu(2*j-1:2*j, (j-i)+1) = Ak_minus_B;
    end
  end
  PhiY0 = Ck * PhiX0;
  PhiYu = Ck * PhiXu;

  % compute desired step location weight
  index = 1;
  planLength = length(stepPlan);
  heelStrikeIndices = footstepPlan.heelStrikeIndices;
  for i = 1:(planLength-1)
    for j = heelStrikeIndices(i):footstepPlan.singleSupportIndices(i+1)
      footPose(index, :) = stepPlan{i}.pose(1:2);
      PhiX0_filtered(index, :) = PhiX0(2*j-1, :);
      PhiXu_filtered(index, :) = PhiXu(2*j-1, :);
      index = index + 1;
    end
  end
  footPose(index, :) = stepPlan{length(stepPlan)}.pose(1:2);
  PhiX0_filtered(index, :) = PhiX0(2*heelStrikeIndices(planLength)-1, :);
  PhiXu_filtered(index, :) = PhiXu(2*heelStrikeIndices(planLength)-1, :);



  % compute optimal control trajectory using QP differential
  xInitial = [dcmInitial; dcmDotInitial];
  xFinal = [dcmTrajectory(N, :); dcmDotTrajectory(N, :)];
  Asolve = Q*PhiYu'*PhiYu + F*PhiXu(end-1:end, :)'*PhiXu(end-1:end, :) + R*eye(N,N);
  for i = 1:2
    bsolve = -(Q*PhiYu'*(PhiY0*xInitial(:, i) - vrpTrajectory(1:N, i)) + ...
        F*PhiXu(end-1:end, :)' * (PhiX0(end-1:end, :) * xInitial(:, i) - xFinal(:, i)));
    U(:, i) = Asolve \ bsolve;
  end

  % compute different optimal control using QP differential
  xInitial = [dcmInitial; dcmDotInitial];
  xFinal = [dcmTrajectory(N, :); dcmDotTrajectory(N, :)];
  Gsolve = zeros(2*N, 2*N);
  gsolve = zeros(2*N, 1);
  Gsolve(1:N, 1:N) = Q*PhiYu'*PhiYu + F*PhiXu(end-1:end, :)'*PhiXu(end-1:end, :) ...
      + R*eye(N, N) + C*PhiXu_filtered'*PhiXu_filtered;
  Gsolve(1:N, N+1:2*N) = -Q*PhiYu';
  Gsolve(N+1:2*N, 1:N) = -Q*PhiYu;
  Gsolve(N+1:2*N, N+1:2*N) = (Q + P) * eye(N, N);
  for i = 1:2
    gsolve(1:N, 1) = Q*PhiYu'*(PhiY0*xInitial(:, i) - copTrajectory(1:N, i)) + ...
        F*PhiXu(end-1:end, :)'*(PhiX0(end-1:end, :)*xInitial(:, i) - xFinal(:, i))...
        + C*PhiXu_filtered'*(PhiX0_filtered*xInitial(:,i) - footPose(:, i));
    gsolve(N+1:2*N, 1) = Q*(copTrajectory(1:N, i) - PhiY0*xInitial(:, i));
    U_alt(:, i) = Gsolve \ (-gsolve);
  end

  knot_x = PhiX0_filtered*xInitial(:,1:2) + PhiXu_filtered*U_alt(1:N,1:2);

  % compute dcm and dcmdot trajectories using new method
  for i = 1:2
    X = PhiXu*U_alt(1:N,i) + PhiX0*xInitial(:, i);
    Y = PhiYu*U_alt(1:N,i) + PhiY0*xInitial(:, i);
    for j = 1:N
      dcmTrajectory(j, i) = X(2*j-1);
      dcmDotTrajectory(j, i) = X(2*j);
      vrpTrajectory(j, i) = Y(j);
    end
  end
end
