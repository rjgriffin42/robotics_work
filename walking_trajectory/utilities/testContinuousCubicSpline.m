clear
clc

plannerDT = 0.005;

xKnots(1) = 0.5;
yKnots(1) = 0.0;

xKnots(2) = 0.3;
yKnots(2) = 0.2;

%{
xKnots(3) = 0.4;
yKnots(3) = 0.6;

xKnots(4) = 0.6;
yKnots(4) = 0.9;

xKnots(5) = 0.8;
yKnots(5) = 0.3;

xKnots(6) = 0.7;
yKnots(6) = 0.0;
%}

duration(1) = 2.0;
duration(2) = 1.5;
duration(3) = 1.0;
duration(4) = 2.5;
duration(5) = 0.7;

numberOfSegments = length(xKnots) - 1;

x_coefficients = computeContinuousCubicCoefficients(xKnots, duration);
y_coefficients = computeContinuousCubicCoefficients(yKnots, duration);

rowIndex = 0;
for i = 1:numberOfSegments
	if i == 1
		t = 0:plannerDT:duration(i);
	else
		t = plannerDT:plannerDT:duration(i);
	end

    for j = 1:length(t)
	    x(rowIndex+j) = x_coefficients{i}(1) + x_coefficients{i}(2) * t(j) + ...
	       x_coefficients{i}(3) * t(j)^2 + x_coefficients{i}(4) * t(j)^3;
	    y(rowIndex+j) = y_coefficients{i}(1) + y_coefficients{i}(2) * t(j) + ...
	       y_coefficients{i}(3) * t(j)^2 + y_coefficients{i}(4) * t(j)^3;
	end
	timeVector(rowIndex+1:rowIndex+length(t)) = t;

    rowIndex = rowIndex + length(t);
end

figure;
subplot(3,1,1)
plot(timeVector, x);
subplot(3,1,2)
plot(timeVector, y);
subplot(3,1,3)
plot(x, y)