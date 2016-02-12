function trajectories = computeCubicSplineTrajectory(xInitial, xFinal, t)
  duration = t(end) - t(1);
  t = (t - t(1)) / duration;
  if (size(t, 2) > size(t, 1))
    t = t';
  end
  time = [t.^0 t.^1 t.^2 t.^3];

  t1 = duration;
  xInitial = xInitial. * [1 t1];
  xFinal = xFinal * [1 t1];
  cubicSplineCoefficients = computeCubicSplineCoefficients(xInitial, xFinal);
  cubicSplineMatrix = [cubicSplineCoefficients(1)    cubicSplineCoefficients(2) / t1;
                       cubicSplineCoefficients(2)  2*cubicSplineCoefficients(3) / t1;
                       cubicSplineCoefficients(3)  3*cubicSplineCoefficients(4) / t1;
                       cubicSplineCoefficients(4)                                 0

  trajectories = time * cubicSplineMatrix; % = [position velocity]
