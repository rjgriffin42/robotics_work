function footstepPlan = computeFootstepPlan(stepPlan, doubleSupportRatio, dt, initialDoubleSupportDuration)

  planLength = length(stepPlan);
  timeVector = computeStepPlanTimeVector(stepPlan, dt);
  
  doubleSupportTime(1) = 0;
  singleSupportTime(1) = initialDoubleSupportDuration;
  heelStrikeTime(1) = (1 - doubleSupportRatio) * stepPlan{1}.duration + initialDoubleSupportDuration;
  stepPlan{1}.duration = initialDoubleSupportDuration + (1 - doubleSupportRatio) * stepPlan{1}.duration;
  time = heelStrikeTime(1);

  for i = 2:planLength
      doubleSupportTime(i) = time;
      singleSupportTime(i) = time + doubleSupportRatio * stepPlan{i}.duration;
      heelStrikeTime(i) = time + stepPlan{i}.duration;
      
      time = time + stepPlan{i}.duration;
  end

  timeVector = 0:dt:time;

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

  footstepPlan.stepPlan = stepPlan;
  footstepPlan.timeVector = timeVector;
  footstepPlan.doubleSupportTime = doubleSupportTime;
  footstepPlan.singleSupportTime = singleSupportTime;
  footstepPlan.heelStrikeTime = heelStrikeTime;
  footstepPlan.doubleSupportIndices = doubleSupportIndices;
  footstepPlan.singleSupportIndices = singleSupportIndices;
  footstepPlan.heelStrikeIndices = heelStrikeIndices;
end
