function [omegaTrajectory, omegaDotTrajectory] = ...
    planOmegaTrajectory(cmpTrajectory, comHeightTrajectory, comDotHeightTrajectory,...
    comDotDotHeightTrajectory, timeVector)

  gravity = 9.81;

  % integrate in reverse time using second order runge-kutta (Huen's method)
  omegaTrajectory = sqrt(gravity/(comHeightTrajectory(end)-cmpTrajectory(end,3))) *...
      ones(length(timeVector), 1);
  omegaDotTrajectory = omegaTrajectory(end)^2 - gravity/(comHeightTrajectory(end)...
      - cmpTrajectory(end, 3)) * ones(length(timeVector), 1);
  for k = (length(timeVector)-1):-1:1
    dt = timeVector(k) - timeVector(k+1);
    k1 = omegaTrajectory(k+1)^2 - (comDotDotHeightTrajectory(k+1) + gravity) / ...
        (comHeightTrajectory(k+1) - cmpTrajectory(k+1, 3));
    k2 = (omegaTrajectory(k+1) + dt * k1)^2 - (comDotDotHeightTrajectory(k) + ...
        gravity) / (comHeightTrajectory(k) - cmpTrajectory(k, 3));
    omegaDotTrajectory(k) = k1 / 2 + k2 / 2;
    omegaTrajectory(k) = omegaTrajectory(k+1) + dt * omegaDotTrajectory(k); 
  end
end
