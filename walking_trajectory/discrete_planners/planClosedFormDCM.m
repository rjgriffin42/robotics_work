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

  for i = planLength:-1:1
    duration = stepPlan{i}.duration;
    dcmEnd{i} = dcmInitial{i+1};

    doubleSupportDuration = duration * doubleSupportRatio;
    singleSupportDuration = duration - doubleSupportDuration;

    alpha = exp(omega0 * singleSupportDuration);
    
    disp('dcmEnd')
    disp(dcmEnd{i})
    disp('vrp')
    disp(vrpKnots{i})
    disp('alpha')
    disp(alpha)

    dcmInitial{i} = (dcmEnd{i} - vrpKnots{i}) / alpha + vrpKnots{i};
    disp('dcmInitial')
    disp(dcmInitial{i})
  end

  %{
  for i = 1:planLength
    disp('vrpKnots')
    disp(vrpKnots{i})
    disp('dcmInitial')
    disp(dcmInitial{i})
  end
  %}

  % compute bang-bang DCM Trajectory
  stepIndex = 1;
  time = stepPlan{stepIndex}.duration;
  baseTime = 0;
  for i = 1:length(timeVector)
    t = timeVector(i);
    duration = stepPlan{stepIndex}.duration;

    dcmTrajectory(i,1:3) = vrpKnots{stepIndex} + ...
      exp(omega0 * (t - baseTime)) * ...
      (dcmInitial{stepIndex} - vrpKnots{stepIndex});
    vrpTrajectory(i,1:3)=vrpKnots{stepIndex};

    if t > time
      time = time + duration;
      baseTime = baseTime + duration;
      stepIndex = stepIndex + 1;
    end
  end
end