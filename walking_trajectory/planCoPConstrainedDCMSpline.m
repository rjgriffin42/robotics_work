include
clear all;

plannerDT = 0.005;
g = 9.81;

left = 0;
right = 1;

swingTime = 0.5;
numberOfSteps = 5;
startSide = right;

% define nominal foot locations
leftBackFoot = [0.25, -0.5];
leftFrontFoot = [0.25, 0.5];
rightFrontFoot = [-0.25, 0.5];
rightBackFoot = [-0.25, -0.5];

slope_left = (leftFrontFoot - rightBackFoot);
slopeStart_left = rightBackFoot;
slope_right = (rightFrontFoot - leftBackFoot);
slopeStart_right = leftBackFoot;

% create group of switching steps
leftLead{1} = leftFrontFoot;
leftLead{2} = rightBackFoot;
rightLead{1} = rightFrontFoot;
rightLead{2} = leftBackFoot;

Q = 1e-3;     % tracking cost
D = 1e10;   % dynamics cost
R = 1e-6;   % acceleration cost

timeVector = 0:plannerDT:(numberOfSteps * swingTime);

% compute desired dcm trajectory
omega = sqrt(g / 0.5);
for i = 1:length(timeVector)
    dcmDesiredTrajectory(i,1:2) = [0 0];
end

% define dynamic parameters
dcmInitial = [-0.1, -0.05];
dcmDotInitial = [0, 0];
knotsPerSegment = 21;
totalKnots = numberOfSteps * (knotsPerSegment - 1) + 1;

currentSide = startSide;
segmentTime = swingTime / (knotsPerSegment - 1);
segmentIndex = 0;
for step = 1:numberOfSteps
    for segment = (segmentIndex+1):((knotsPerSegment-1)+segmentIndex)
        segmentDuration(segment) = segmentTime;
        
        if (currentSide == left)
            y_x(segment) = slope_left(1);
            y0_x(segment) = slopeStart_left(1);
            y_y(segment) = slope_left(2);
            y0_y(segment) = slopeStart_left(2);
        else
            y_x(segment) = slope_right(1);
            y0_x(segment) = slopeStart_right(1);
            y_y(segment) = slope_right(2);
            y0_y(segment) = slopeStart_right(2);
        end
    end
    
    segmentIndex = segmentIndex+(knotsPerSegment-1);
    if step == numberOfSteps
        
    end
    if (currentSide == left)
        currentSide = right;
    else
        currentSide = left;
    end
end



if (currentSide == left)
    y_x(totalKnots) = slope_left(1);
    y0_x(totalKnots) = slopeStart_left(1);
    y_y(totalKnots) = slope_left(2);
    y0_y(totalKnots) = slopeStart_left(2);
else
    y_x(totalKnots) = slope_right(1);
    y0_x(totalKnots) = slopeStart_right(1);
    y_y(totalKnots) = slope_right(2);
    y0_y(totalKnots) = slopeStart_right(2);
end
        
y_x_temp = y_x;
y_y_temp = y_y;
y_x = zeros(totalKnots); y_y = zeros(totalKnots);
for i = 1:length(y_x_temp)
    y_x(i,i) = y_x_temp(i);
    y_y(i,i) = y_y_temp(i);
end

% create projection matrices
[betaMatrices, betaDotMatrices, betaDotDotMatrices] = ...
    computeCubicProjection(segmentDuration);

% create discretized projection
rowIndex = 0;
for segment = 1:totalKnots-1
    if segment == 1
        t = 0:plannerDT:segmentDuration(segment);
    else
        t = plannerDT:plannerDT:segmentDuration(segment);
    end
    T = [t'.^0, t'.^1, t'.^2, t'.^3];
    Phi(rowIndex+1:rowIndex+size(T,1),:) = T * betaMatrices{segment};
    PhiDot(rowIndex+1:rowIndex+size(T,1),:) = T(:,1:3) * betaDotMatrices{segment};
    PhiDotDot(rowIndex+1:rowIndex+size(T,1),:) = T(:,1:2) * betaDotDotMatrices{segment};
    rowIndex = rowIndex + size(T,1);
end


% add constraints on dcm initial and final value
CE = zeros(4, 3*totalKnots);
CE(1,1) = 1;
CE(2,totalKnots) = 1;
CE(3,totalKnots+1) = 1;
CE(4,2*totalKnots) = 1;
ce(1) = dcmInitial(1);
ce(2) = dcmDesiredTrajectory(end,1);
ce(3) = dcmInitial(2);
ce(4) = dcmDesiredTrajectory(end,2);

% add constraints on CoP position
for i = 1:totalKnots
    CI((2*i-1),2*totalKnots+i) = 1;
    CI((2*i),2*totalKnots+i) = -1;
    ci((2*i-1),:) = 1;
    ci((2*i),:) = 0;
end

G00 = Q * Phi' * Phi + D * (PhiDot - omega * Phi)' * (PhiDot - omega * Phi) + R * PhiDotDot' * PhiDotDot;
G11 = G00;
G01 = zeros(totalKnots);
G02 = D * (PhiDot - omega * Phi)' * (omega * Phi * y_x);
G12 = D * (PhiDot - omega * Phi)' * (omega * Phi * y_y);
G22 = omega^2 * D * ((Phi * y_x)' * (Phi * y_x) + (Phi * y_y)' * (Phi * y_y));

G = [G00 G01 G02;...
     G01' G11 G12;...
     G02' G12' G22];

g0 = 2 * (D * omega * (Phi * y0_x')' * (PhiDot - omega * Phi) - Q * dcmDesiredTrajectory(:,1)' * Phi);
g1 = 2 * (D * omega * (Phi * y0_y')' * (PhiDot - omega * Phi) - Q * dcmDesiredTrajectory(:,2)' * Phi);
g2 = 2 * omega^2 * D * ((Phi * y0_x')' * (Phi * y_x) + (Phi * y0_y')' * (Phi * y_y));
g = [g0 g1 g2];

u = quadprog(G, g', CI, ci', CE, -ce');

dcmKnots(:,1) = u(1:totalKnots);
dcmKnots(:,2) = u(totalKnots+1:2*totalKnots);

alpha = u(2*totalKnots+1:end);
copKnots(:,1) = y0_x' + y_x * alpha;
copKnots(:,2) = y0_y' + y_y * alpha;

dcmTrajectory(:,1:2) = Phi * dcmKnots;
dcmDotTrajectory(:,1:2) = PhiDot * dcmKnots;

copTrajectory(:,1:2) = Phi * copKnots;

for i = 1:length(dcmTrajectory)
    vrpTrajectory(i,1:2) = dcmTrajectory(i,:) - dcmDotTrajectory(i,:) / omega;
end

figure;
subplot(3,1,1)
plot(timeVector, dcmTrajectory(:,1), timeVector, copTrajectory(:,1), timeVector, vrpTrajectory(:,1))
subplot(3,1,2)
plot(timeVector, dcmTrajectory(:,2), timeVector, copTrajectory(:,2), timeVector, vrpTrajectory(:,2))
subplot(3,1,3)
plot(dcmTrajectory(:,1), dcmTrajectory(:,2), copTrajectory(:,1), copTrajectory(:,2), vrpTrajectory(:,1), vrpTrajectory(:,2))
