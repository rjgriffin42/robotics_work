function [s2, k2, alpha] = computeLQRTimeVaryingLinearTerm(S11, S12, S13, B2, Q, B, ...
    D, R1, coefficients, timeKnots, plannerDT)

  R1inv = inv(R1);
  temp = (S12 + Q * D) * R1inv;
  A2 = [zeros(2) temp; -eye(2) S13' * R1inv];
  A2inv = [inv(temp) * S13' * R1inv -eye(2); inv(temp) zeros(2)];

  n = length(coefficients);
  k = size(coefficients{1},2);

  % compute coefficients
  for j = n:-1:1
  	beta{j}(:,k) = -A2inv * B2 * coefficients{j}(:,k);
    gamma{j}(:,k) = R1inv * D * Q * coefficients{j}(:,k) -...
        1/2 * R1inv * B' * beta{j}(:,k);

    for i = (k-1):-1:1
    	beta{j}(:,i) = A2inv * ((i+1) * beta{j}(:,i+1) - ...
            B2 * coefficients{j}(:,i));
        gamma{j}(:,i) = R1inv * D * Q * coefficients{j}(:,i) - ...
            1/2 * R1inv * B' * beta{j}(:,i);
    end
   
    beta_sum = 0;
    for index = 1:k
        beta_sum = beta_sum + beta{j}(:,index) * ...
            (timeKnots{j+1} - timeKnots{j})^(index-1);
    end
    if j == n
        alpha{j} = -exp(A2 * (timeKnots{j+1} - timeKnots{j})) \ beta_sum;    
    else
        alpha{j} = exp(A2 * (timeKnots{j+1} - timeKnots{j})) \ ...
            (alpha{j+1} + beta{j+1}(:,1) - beta_sum);
    end
  end
  
  % compute trajectory solution
  rowIndex = 0;
  for j = 1:n
      t = timeKnots{j}:plannerDT:(timeKnots{j+1}-plannerDT);
      beta_sum = 0;
      gamma_sum = 0;
      for i = 1:length(t)
          for index = 1:k
              beta_sum = beta_sum + beta{j}(:,index) * (t(i) - timeKnots{j})^(index-1);
              gamma_sum = gamma_sum + gamma{j}(:,index) * (t(i) - timeKnots{j})^(index-1);
          end
          s2(:,rowIndex+i) = exp(A2 * (t(i) - timeKnots{j})) * alpha{j} + beta_sum;
          k2 = -1/2 * R1inv * B' * exp(A2 * (t(i) - timeKnots{j})) * alpha{j} + gamma_sum;
      end
      rowIndex = rowIndex + length(t);
  end
end