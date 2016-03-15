function c = computeContinuousCubicCoefficients(xInitial, xFinal, duration)
  c(1:2, 1) = xInitial(1:2, 1);
  c(1:2, 2) = xInitial(1:2, 2);
  c(1:2, 3) = 3 / duration^2 * (xFinal(1:2, 1) - xInitial(1:2, 1))...
      - 1 / (2 * duration) * (xFinal(1:2, 2) + 4 * xInitial(1:2, 2));
  c(1:2, 4) = 2 / duration^3 * (xInitial(1:2, 1) - xFinal(1:2, 1)) + ...
      1 / duration^2 * (xFinal(1:2, 2) + xInitial(1:2, 2));
end