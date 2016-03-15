function [s2] = computeLQRTimeVaryingLinearTerm(S11, S12, S13, B2, Q, D, R1, c)
  R1inv = inv(R1);
  temp = (S12 + Q * D) * R1inv;
  A2 = [0 temp; -eye(2) S13' * inv(R1)];
  A2inv = [inv(temp) * S13' * R1inv -eye(2); inv(temp) 0];

  n = rows(c);
  k = cols(c);
  for j = n:-1:1
  	beta(j,:) = -inv(A2) * B * c(j,:);
    gamma(j,:) = R1inv * D * Q * c(j,:) - 1/2 * R1inv * B' * beta(j,:);

    for i = (k-1):-1:0
    	beta(j,i) = A2inv * ((i+1) * beta(j,i+1) - B2 * c(j,i));
        gamma(j,i) = R1inv * D * Q * c(j,i) - 1/2 * R1inv * B' * beta(j,i);
    end
  end
end