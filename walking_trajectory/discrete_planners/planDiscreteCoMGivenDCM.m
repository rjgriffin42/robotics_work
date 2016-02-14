function [comTrajectory, comDotTrajectory] = planCoMGivenDCM(dcmTrajectory, ...
    omegaTrajectory, comInitial, timeVector)

  comTrajectory = zeros(length(timeVector), 3);
  comDotTrajectory = zeros(length(timeVector), 3);
  comTrajectory(1, :) = comInitial;
  comDotTrajectory(1, :) = omegaTrajectory(1) * (dcmTrajectory(1, :) - ...
      comTrajectory(1, :));

  for i = 1:length(timeVector)-1
    dt = timeVector(i+1) - timeVector(i);
    k1 = omegaTrajectory(i) * (dcmTrajectory(i, :) - comTrajectory(i, :));
    k2 = omegaTrajectory(i+1) * (dcmTrajectory(i+1, :) - comTrajectory(i, :) -...
        dt * k1);
    comDotTrajectory(i+1, :) = k1/2 + k2/2;
    comTrajectory(i+1, :) = comTrajectory(i, :) + dt * comDotTrajectory(i+1, :);
  end
end
