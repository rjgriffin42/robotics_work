function cubicSplineCoefficients = computeCubicSplineCoefficients(xInitial, xFinal)
  coefficients = [1.0  0.0  0.0  0.0
                  0.0  1.0  0.0  0.0
                 -3.0 -2.0  3.0 -1.0
                  2.0  1.0 -2.0  1.0];
  boundaryConditions = [xInitial(1) xInitial(2) xFinal(1) xFinal(2)]';
  cubicSplineCoefficients = coefficients * boundaryConditions;
end
