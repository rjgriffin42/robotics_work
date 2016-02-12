function trajectories = compute_minimum_jerk_trajectory(x_initial, x_final, t)
  % Time Matrix
  duration = t(end) - t(1);
  t = (t - t(1)) / duration;
  if (size(t, 2) > size(t, 1))
    t = t';
  end
  times = [t.^0 t.^1 t.^2 t.^3 t.^4 t.^5];
  
  % minimum jerk coefficient matrix
  x_initial = x_initial. * [1 duration duration^2];
  x_final = x_final. * [1 duration duration^2];

  b = compute_minimum_jerk_coefficients(x_initial, x_final);
  B = [b(1)    b(2) / duration   2*b(3) / duration^2;
       b(2)  2*b(3) / duration   6*b(4) / duration^2;
       b(3)  3*b(4) / duration  12*b(5) / duration^2;
       b(4)  4*b(5) / duration  20*b(6) / duration^2;
       b(5)  5*b(6) / duration        0 / duration^2;
       b(6)       0 / duration        0 / duration^2];

  % compute position, velocity, and acceleration trajectories
  trajectories = times * B;
end
