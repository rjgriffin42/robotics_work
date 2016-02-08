% Author: Robert Griffin (rjgriffin42@gmail.com)
% Date: 2/7/14

include;
clear all;

% define constants
dt = 0.005;
mass = 75; % kg

% define gait parameters
ds_ratio = 0.2; % double support ratio
zcom_nominal = 0.85; % nominal com height;
to_ratio = 0.7;
hs_ratio = 0.0;

% define LQR parameters

% define initial conditions
l_foot_pose0 = [0, 0.1, 0, 0, 0, 0];
r_foot_pose0 = [0,-0.1, 0, 0, 0, 0];
com0 = [0, 0, 0.85];
comdot0 = [0, 0, 0];
comdotdot0 = [0, 0, 0];
cop0 = [com0(1) com0(2) (l_foot_pose(3) + r_foot_pose(3)];
omega0 = 1/sqrt(9.81/com0(3));
dcm0 = com0 + omega0*comdot0;
dcmdot0 = zeros(1,3);

% define step plan
step_plan = forward_step_plan();

% compute discrete time vector
t = compute_step_plan_time_vector(step_plan, dt);

% compute inertial step plan
step_plan = transform_step_plan_to_inertial_frame(l_foot_pose0, r_foot_pose0, ...
  step_plan);


