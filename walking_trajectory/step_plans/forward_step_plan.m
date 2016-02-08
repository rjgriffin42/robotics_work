function step_plan = forward_step_plan()
  step_plan = {};
  step_plan{1} = new_step('r', [0.25,-0.2, 0.0, 0.0, 0.0, 0], 1.0);
  step_plan{2} = new_step('l', [0.25, 0.2, 0.0, 0.0, 0.0, 0], 1.0);
  step_plan{3} = new_step('r', [0.25,-0.2, 0.0, 0.0, 0.0, 0], 1.0);
  step_plan{4} = new_step('l', [0.25, 0.2, 0.0, 0.0, 0.0, 0], 1.0);
end
