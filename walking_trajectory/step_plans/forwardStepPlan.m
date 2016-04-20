function stepPlan = forwardStepPlan()
  stepPlan = {};
  stepPlan{1} = newStep('r', [0.9,-0.2, 0.0, 0.0, 0.0, 0], 0.5);
  stepPlan{2} = newStep('l', [0.9, 0.2, 0.0, 0.0, 0.0, 0], 0.5);
  stepPlan{3} = newStep('r', [0.9,-0.2, 0.0, 0.0, 0.0, 0], 0.5);
  stepPlan{4} = newStep('l', [0.9, 0.2, 0.0, 0.0, 0.0, 0], 0.5);
end
