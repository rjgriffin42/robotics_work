function b = compute_minimum_jerk_coefficients(x_initial, x_final)
  % returns 'b_k' coefficients for minimum jerk trajectory computation
  % see "Jerk-Bounded Manipulator Trajectory Planning: Design for
  % Real-Time Applications," S. Macfarlane and E. Croft
  % x_initial = [x(0), xdot(0), xdotdot(0)]
  % x_final   = [x(0), xdot(0), xdotdot(0)]

  A = [1.0  0.0  0.0  0.0  0.0  0.0;
       0.0  1.0  0.0  0.0  0.0  0.0;
       0.0  0.0  0.5  0.0  0.0  0.0;
      -10  -6.0 -1.5  10  -4.0  0.5;
       15   8.0  1.5 -15   7.0 -1.0;
      -6.0 -3.0 -0.5  6   -3.0  0.5];
  x = [x_initial(1) x_initial(2) x_initial(3) x_final(1) x_final(2) x_final(3)]';
  b = A*x;
end
