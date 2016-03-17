clear
clc

plannerDT = 0.005;

xKnots{1} = [0.5; 0.0];
xKnots{2} = [0.3; 0.2];
%{
xKnots{3} = [0.4; 0.6];
xKnots{4} = [0.6; 0.9];
xKnots{5} = [0.8; 0.3];
xKnots{6} = [0.7; 0.0];
%}

duration(1) = 2.0;
duration(2) = 1.5;
duration(3) = 1.0;
duration(4) = 2.5;
duration(5) = 0.7;

numberOfSegments = length(xKnots) - 1;

coefficients = computeContinuousCubicCoefficients(xKnots, duration);

rowIndex = 0;
for i = 1:numberOfSegments
	if i == 1
		t = 0:plannerDT:duration(i);
	else
		t = plannerDT:plannerDT:duration(i);
	end

    for j = 1:length(t)
    	x(:,rowIndex+j) = coefficients{i}(:,1) + coefficients{i}(:,2) * t(j) + ...
           coefficients{i}(:,3) * t(j)^2 + coefficients{i}(:,4) * t(j)^3;
	end
	timeVector(rowIndex+1:rowIndex+length(t)) = t;

    rowIndex = rowIndex + length(t);
end

figure;
subplot(3,1,1)
plot(timeVector, x(1,:));
subplot(3,1,2)
plot(timeVector, x(2,:));
subplot(3,1,3)
plot(x(1,:), x(2,:))