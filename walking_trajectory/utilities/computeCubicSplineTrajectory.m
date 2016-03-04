function trajectories = computeCubicSplineTrajectory(xInitial, xFinal, timeVector)
  duration = timeVector(end) - timeVector(1);
  timeVector = (timeVector - timeVector(1)) / duration;
  if (size(timeVector, 2) > size(timeVector, 1))
    timeVector = timeVector';
  end
  timeArray = [timeVector.^0 timeVector.^1 timeVector.^2 timeVector.^3];

  t1 = duration;
  xInitial = xInitial .* [1 t1];
  xFinal = xFinal .* [1 t1];
  cubicSplineCoefficients = computeCubicSplineCoefficients(xInitial, xFinal);
  cubicSplineMatrix = [cubicSplineCoefficients(1)    cubicSplineCoefficients(2) / t1;
                       cubicSplineCoefficients(2)  2*cubicSplineCoefficients(3) / t1;
                       cubicSplineCoefficients(3)  3*cubicSplineCoefficients(4) / t1;
                       cubicSplineCoefficients(4)                                 0];

  trajectories = timeArray * cubicSplineMatrix; % = [position velocity]
end
