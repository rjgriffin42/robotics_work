function cmpTrajectory = planCMPTrajectory(copTrajectory, tauTrajectory, comDotDotHeightTrajectory, mass, timeVector)

  gravity = 9.81

  cmpTrajectory = zeros(size(copTrajectory));
  for i = 1:length(timeVector)
    verticalForce = mass * (comDotDotHeightTrajectory(i) + gravity;
    cmpTrajectory(i, 1) = copTrajectory(i, 1) + tauTrajectory(i, 2) / verticalForce;
    cmpTrajectory(i, 2) = copTrajectory(i, 2) - tauTrajectory(i, 2) / verticalForce;
    cmpTrajectory(i, 3) = copTrajectory(i, 3);
  end
end

