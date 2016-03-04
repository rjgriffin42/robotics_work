%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based off of the work presented in
%    Johannes Englsberger, Twan Koolen, Sylvain Bertrand, Jerry Pratt, Christian
%      Ott, and Alin Albu-Shcaffer, "Trajectory generation for continuous leg 
%      forces during double support and heel-to-toe shift based on divergent
%      component of motion." 2014 IEEE/RSJ International Conference on 
%      Intelligent Robots and Systems. 2014.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcmTrajectory, vrpTrajectory] = planClosedFormDCM(footstepPlan,...
    doubleSupportRatio, omega0, timeVector)

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

  % compute whole DCM Trajectory
  stepIndex = 1;
  time = doubleSupportDuration{stepIndex};
  baseTime = 0;
  inDoubleSupport = true;

  for i = 1:length(timeVector)
    t = timeVector(i);

    if inDoubleSupport == true
      dcmTrajectory(i, 1:3) = dcmInitial{stepIndex};
      vrpTrajectory(i, 1:3) = vrpKnots{stepIndex};
    else
      duration = stepPlan{stepIndex}.duration;

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
end
