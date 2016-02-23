function footstepPlan = computeFootstepPlan(stepPlan, doubleSupportRatio, dt)

  planLength = length(stepPlan);
  timeVector = computeStepPlanTimeVector(stepPlan, dt);

  % copy existing step plan
  footstepPlan.stepPlan = stepPlan;
  
  time = 0;
  for i = 1:planLength
      doubleSupportTime(i) = time;
      singleSupportTime(i) = time + doubleSupportRatio * stepPlan{i}.duration;
      heelStrikeTime(i) = time + stepPlan{i}.duration;
      
      time = time + stepPlan{i}.duration;
  end

  for i = 1:length(timeVector)
    for j = 1:planLength
      if (doubleSupportTime(j) >= timeVector(i))
        doubleSupportIndices(j) = i;
      end
      
      if (singleSupportTime(j) >= timeVector(i))
        singleSupportIndices(j) = i;
      end

      if (heelStrikeTime(j) >= timeVector(i))
        heelStrikeIndices(j) = i;
      end
    end
  end

  footstepPlan.timeVector = timeVector;
  footstepPlan.doubleSupportTime = doubleSupportTime;
  footstepPlan.singleSupportTime = singleSupportTime;
  footstepPlan.heelStrikeTime = heelStrikeTime;
  footstepPlan.doubleSupportIndices = doubleSupportIndices;
  footstepPlan.singleSupportIndices = singleSupportIndices;
  footstepPlan.heelStrikeIndices = heelStrikeIndices;
end
