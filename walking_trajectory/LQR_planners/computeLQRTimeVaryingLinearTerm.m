function [s2, k2, alpha] = computeLQRTimeVaryingLinearTerm(S11, S12, S13, A2, B2, Q, B, ...
    D, R1, coefficients, timeKnots, plannerDT)

  R1inv = inv(R1);
  temp = (S12 + Q * D) * R1inv;
  A2inv = inv(A2);

  n = length(coefficients);
  k = size(coefficients{1},2);

  % compute coefficients
  for j = n:-1:1
  	beta(:,j,k) = -A2inv * B2 * coefficients{j}(:,k);
    gamma(:,j,k) = R1inv * D * Q * coefficients{j}(:,k) -...
        1/2 * R1inv * B' * beta(:,j,k);

    for i = (k-1):-1:1
    	beta(:,j,i) = A2inv * (i * beta(:,j,i+1) - ...
            B2 * coefficients{j}(:,i));
        gamma(:,j,i) = R1inv * D * Q * coefficients{j}(:,i) - ...
            1/2 * R1inv * B' * beta(:,j,i);
    end
   
    beta_sum = 0;
    dt = timeKnots{j+1} - timeKnots{j}
    for index = 1:k
        beta_sum = beta_sum + beta{j}(:,index) * ...
            dt^(index-1);
    end
    if (j == n)
        s2dt = zeros(4,1)
    else
        s2dt = alpha(:,j+1) + beta(:,j+1,1);
    end
    alpha(:,j) = expm(A2 * dt) \ (s1dt - beta_sum);
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
