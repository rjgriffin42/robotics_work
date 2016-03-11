%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based off of the work presented in
%    Johannes Englsberger, Twan Koolen, Sylvain Bertrand, Jerry Pratt, Christian
%      Ott, and Alin Albu-Shcaffer, "Trajectory generation for continuous leg 
%      forces during double support and heel-to-toe shift based on divergent
%      component of motion." 2014 IEEE/RSJ International Conference on 
%      Intelligent Robots and Systems. 2014.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcmTrajectory, vrpTrajectory] = planClosedFormDCM(footstepPlan, dcm0, omega0, plannerDT)

  stepPlan = footstepPlan.stepPlan;
  alphaPlan = 0.5;

  planLength = length(stepPlan);
  for i = 1:planLength
    vrpKnots{i+1} = stepPlan{i}.pose(1:3);
  end
  vrpKnots{1} = dcm0;

  % add in first one for the initial foot position
  dcmFinal = 0.5 * (stepPlan{planLength}.pose(1:3) + stepPlan{planLength-1}.pose(1:3));
  dcmInitial{planLength+1} = dcmFinal;

  %compute knot points for DCM trajectory single support
  for i = planLength:-1:1
    duration{i} = stepPlan{i}.duration;
    dcmEnd{i} = dcmInitial{i+1};

    alpha = exp(omega0 * duration{i});
    
    dcmInitial{i} = (dcmEnd{i} - vrpKnots{i}) / alpha + vrpKnots{i};
  end

  % compute whole DCM Trajectory
  baseIndex = 0;
  for stepIndex = 1:planLength
    if stepIndex == 1
      t = 0:plannerDT:stepPlan{stepIndex}.duration;
    else
      t = plannerDT:plannerDT:stepPlan{stepIndex}.duration;
    end

    for j = 1:length(t)
      dcmTrajectory(baseIndex+j,1:3) = vrpKnots{stepIndex} + exp(omega0 * t(j)) * (dcmInitial{stepIndex} - vrpKnots{stepIndex});
      dcmDotTrajectory(baseIndex+j, 1:3) = omega0 * (dcmTrajectory(baseIndex+j, 1:3) - vrpKnots{stepIndex});
      vrpTrajectory(baseIndex+j,1:3) = vrpKnots{stepIndex};
    end
    baseIndex = baseIndex + length(t);
  end
end
