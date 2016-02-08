function t = compute_step_plan_time_vector(step_plan, dt)
  T = 0;
  for i = 1:length(step_plan)
    T = T + step_plan{i}.duration;
  end
  t = [0:dt:T];
end
