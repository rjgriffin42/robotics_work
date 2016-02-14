function timeVector = computeStepPlanTimeVector(stepPlan, dt)
  timeEnd = 0;
  for i = 1:length(stepPlan)
    timeEnd = timeEnd + stepPlan{i}.duration;
  end
  timeVector = [0:dt:timeEnd];
end
