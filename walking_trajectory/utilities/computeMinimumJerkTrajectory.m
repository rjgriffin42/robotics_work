function trajectories = computeMinimumJerkTrajectory(xInitial, xFinal, timeVector)
  % Time Matrix
  duration = timeVector(end) - timeVector(1);
  timeVector = (timeVector - timeVector(1)) / duration;
  if (size(timeVector, 2) > size(timeVector, 1))
    timeVector = timeVector';
  end
  timeArray = [timeVector.^0 timeVector.^1 timeVector.^2 timeVector.^3 timeVector.^4 timeVector.^5];
  
  % minimum jerk coefficient matrix
  xInitial = xInitial .* [1 duration duration^2];
  xFinal = xFinal .* [1 duration duration^2];

  coefficients = computeMinimumJerkCoefficients(xInitial, xFinal);
  B = [coefficients(1)    coefficients(2) / duration   2*coefficients(3) / duration^2;
       coefficients(2)  2*coefficients(3) / duration   6*coefficients(4) / duration^2;
       coefficients(3)  3*coefficients(4) / duration  12*coefficients(5) / duration^2;
       coefficients(4)  4*coefficients(5) / duration  20*coefficients(6) / duration^2;
       coefficients(5)  5*coefficients(6) / duration        0 / duration^2;
       coefficients(6)       0 / duration        0 / duration^2];

  % compute position, velocity, and acceleration trajectories
  trajectories = timeArray * B;
end
