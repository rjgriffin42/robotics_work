function vrpTrajectory = planVRPTrajectory(cmpTrajectory, omegaTrajectory, ...
    omegaDotTrajectory, timeVector)

    gravity = 9.81;

    heightInitial = gravity * ones(length(timeVector), 1) ./ ...
        (omegaTrajectory.^2 - omegaDotTrajectory);
    vrpTrajectory = cmpTrajectory + [zeros(length(timeVector), 1), ...
        zeros(length(timeVector), 1) heightInitial];
end
