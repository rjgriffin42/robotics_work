%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based off of the work presented in
%    Johannes Englsberger, Twan Koolen, Sylvain Bertrand, Jerry Pratt, Christian
%      Ott, and Alin Albu-Shcaffer, "Trajectory generation for continuous leg 
%      forces during double support and heel-to-toe shift based on divergent
%      component of motion." 2014 IEEE/RSJ International Conference on 
%      Intelligent Robots and Systems. 2014.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcmTrajectory, vrpTrajectory] = planCDSClosedFormDCM(footstepPlan,...
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

  %compute knot points for DCM trajectory single support
  for i = planLength:-1:1
    duration{i} = stepPlan{i}.duration;
    dcmEnd{i} = dcmInitial{i+1};

    doubleSupportDuration{i} = duration{i} * doubleSupportRatio;
    singleSupportDuration{i} = duration{i} - doubleSupportDuration{i};

    alpha = exp(omega0 * singleSupportDuration{i});
    
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
    tDeltaIni{i} = -t_alpha * doubleSupportDuration{i};
    tDeltaEnd{i} = (1 - t_alpha) * doubleSupportDuration{i};

    dcmIniDS{i} = vrpKnots{i-1}+exp(omega0*tDeltaIni{i})*(dcmInitial{i}-vrpKnots{i-1});
    dcmEoDS{i} = vrpKnots{i}+exp(omega0*tDeltaEnd{i})*(dcmInitial{i}-vrpKnots{i});
    dcmDotIniDS{i} = 1 / omega0 * (dcmIniDS{i} - vrpKnots{i-1});
    dcmDotEoDS{i} = 1 / omega0 * (dcmEoDS{i} - vrpKnots{i});
  end

  for i = 1:planLength
    disp(i)
    disp(dcmIniDS{i});
    disp(dcmDotIniDS{i});
    disp(dcmEoDS{i});
    disp(dcmDotEoDS{i});
  end

  % compute whole DCM Trajectory
  stepIndex = 1;
  time = doubleSupportDuration{stepIndex};
  baseTime = 0;
  inDoubleSupport = true;

  index = 0;
  for i = 1:planLength
    t_ds = 0:plannerDT:doubleSupportDuration{i};
    disp(doubleSupportDuration{i});
    disp(i)

    for j = 1:3
      p = computeCubicSplineTrajectory([dcmIniDS{i}(j) dcmDotIniDS{i}(j)],...
                                       [dcmEoDS{i}(j) dcmDotEoDS{i}(j)], t_ds);
      dcmTrajectory(index+1:(length(t_ds)+index), j) = p(:, 1);
    end
    index = index + length(t_ds);

    t_ss = tDeltaEnd{i}:plannerDT:(singleSupportDuration{i}+tDeltaIni{i+1});
    for j = 1:length(t_ss)
      dcmTrajectory(index+j,1:3) = vrpKnots{i} + exp(omega0 * t_ss(j)) - vrpKnots{i};
    end
    index = index + length(t_ss);
  end

  %{
  for i = 1:length(timeVector)
    t = timeVector(i);

    if inDoubleSupport == true
      cubic_spline = computeCubicSplineTrajectory

      %dcmTrajectory(i, 1:3) = dcmInitial{stepIndex};
      vrpTrajectory(i, 1:3) = dcmTrajectory(i,1:3)-1/omega0*dcmDotTrajectory(i,1:3);
    else
      dcmTrajectory(i,1:3) = vrpKnots{stepIndex} + ...
        exp(omega0 * (t - baseTime)) * ...
        (dcmInitial{stepIndex} - vrpKnots{stepIndex});
      vrpTrajectory(i,1:3)=vrpKnots{stepIndex};
    end

    if t > time
      if inDoubleSupport == true
        time = time + singleSupportDuration{stepIndex};
        baseTime = baseTime + doubleSupportDuration{stepIndex};
        inDoubleSupport = false;
      else
        time = time + doubleSupportDuration{stepIndex};
        baseTime = baseTime + singleSupportDuration{stepIndex};
        inDoubleSupport = true;
        stepIndex = stepIndex + 1;
      end
    end
  end
  %}
end
