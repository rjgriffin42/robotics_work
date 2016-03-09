%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based off of the work presented in
%    Johannes Englsberger, Twan Koolen, Sylvain Bertrand, Jerry Pratt, Christian
%      Ott, and Alin Albu-Shcaffer, "Trajectory generation for continuous leg 
%      forces during double support and heel-to-toe shift based on divergent
%      component of motion." 2014 IEEE/RSJ International Conference on 
%      Intelligent Robots and Systems. 2014.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcmTrajectory, dcmDotTrajectory, vrpTrajectory] = planCDSClosedFormDCM(footstepPlan,...
    doubleSupportRatio, t_alpha, omega0, plannerDT, timeVector)

  stepPlan = footstepPlan.stepPlan;
  alphaPlan = 0.5;

  planLength = length(stepPlan);
  for i = 1:planLength
    vrpKnots{i+1} = stepPlan{i}.pose(1:3);
  end
  vrpKnots{1} = [0 0.1 0];
  % add in first one for the initial foot position
  dcmInitial{planLength+1} = vrpKnots{planLength+1};

  % compute knot points for DCM trajectory single support
  for i = planLength:-1:1
    duration{i} = stepPlan{i}.duration;
    dcmEnd{i} = dcmInitial{i+1};

    doubleSupportDuration{i} = duration{i} * doubleSupportRatio;
    singleSupportDuration{i} = duration{i} - doubleSupportDuration{i};

    alpha = exp(omega0 * duration{i});
    
    dcmInitial{i} = (dcmEnd{i} - vrpKnots{i}) / alpha + vrpKnots{i};
  end

  % compute points for continuous double support trajectory
  tDeltaIni{1} = 0;
  tDeltaEnd{1} = (1 - t_alpha) * doubleSupportDuration{1};
  dcmIniDS{1} = dcmInitial{1};
  dcmEoDS{1} = vrpKnots{1} + exp(omega0*tDeltaEnd{1}) * (dcmInitial{1} -vrpKnots{1});
  dcmDotIniDS{1} = 1 / omega0 * (dcmIniDS{1} - vrpKnots{1});
  dcmDotEoDS{1} = 1 / omega0 * (dcmEoDS{1} - vrpKnots{1});

  for i = 2:planLength
    tDeltaIni{i} = t_alpha * doubleSupportDuration{i};
    tDeltaEnd{i} = (1 - t_alpha) * doubleSupportDuration{i};

    dcmIniDS{i} = vrpKnots{i-1} + exp(-omega0 * tDeltaIni{i}) * (dcmInitial{i} - vrpKnots{i-1});
    dcmEoDS{i} = vrpKnots{i} + exp(omega0 * tDeltaEnd{i}) * (dcmInitial{i} - vrpKnots{i});

    dcmDotIniDS{i} = omega0 * (dcmIniDS{i} - vrpKnots{i-1});
    dcmDotEoDS{i} = omega0 * (dcmEoDS{i} - vrpKnots{i});
  end
  tDeltaIni{planLength+1} = t_alpha * doubleSupportDuration{1};
  tDeltaEnd{planLength+1} = 0;

  dcmIniDS{planLength+1} = vrpKnots{planLength} + ...
      exp(-omega0 * tDeltaIni{planLength+1}) * (dcmInitial{planLength+1} - vrpKnots{planLength});
  dcmEoDS{planLength+1} = vrpKnots{planLength+1};

  dcmDotIniDS{planLength+1} = omega0 * (dcmIniDS{planLength+1} - vrpKnots{planLength});
  dcmDotEoDS{planLength+1} = [0 0 0];

  % compute whole DCM Trajectory
  index = 0;
  for stepIndex = 1:planLength
    if stepIndex == 1;
      t_ds = 0:plannerDT:tDeltaEnd{1};
      %t_ds = 0:plannerDT:doubleSupportDuration{stepIndex};
    else
      t_ds = plannerDT:plannerDT:doubleSupportDuration{stepIndex};
    end
    for j = 1:3
      tempVector = computeCubicSplineTrajectory([dcmIniDS{stepIndex}(j) dcmDotIniDS{stepIndex}(j)], ...
        [dcmEoDS{stepIndex}(j) dcmDotEoDS{stepIndex}(j)], t_ds);
      startIndex = index + 1;
      endIndex = index + length(t_ds);
      dcmTrajectory(startIndex:endIndex,j) = tempVector(:,1);
      dcmDotTrajectory(startIndex:endIndex,j) = tempVector(:,2);
    end

    index = index + length(t_ds);

    t_ss = plannerDT:plannerDT:singleSupportDuration{stepIndex};
    for j = 1:length(t_ss)
      dcmTrajectory(index+j, 1:3) = vrpKnots{stepIndex} + exp(omega0 * t_ss(j)) * (dcmEoDS{stepIndex} - vrpKnots{stepIndex});
      dcmDotTrajectory(index+j, 1:3) = omega0 * (dcmTrajectory(index+j, 1:3) - vrpKnots{stepIndex});
    end

    index = index + length(t_ss);
  end

  t_ds = plannerDT:plannerDT:tDeltaIni{planLength+1};
  for j = 1:3
    tempVector = computeCubicSplineTrajectory([dcmIniDS{planLength+1}(j) dcmDotIniDS{planLength+1}(j)], ...
      [dcmEoDS{planLength+1}(j) dcmDotEoDS{planLength+1}(j)], t_ds);
    startIndex = index + 1;
    endIndex = index + length(t_ds);
    dcmTrajectory(startIndex:endIndex,j) = tempVector(:,1);
    dcmDotTrajectory(startIndex:endIndex,j) = tempVector(:,2);
  end

  vrpTrajectory = dcmTrajectory - 1 / omega0 * dcmDotTrajectory;
end
