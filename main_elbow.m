% Example script for solving for translation and rotational alignment, and using these
% paramenters to run mekf-acc, rts-acc, mekf-dof, and rts-dof for the elbow joint 
% must specify ua.acc, ua.gyr, fa.acc, fa.gyr (accelerometer and gyro
% measurements for the upper arm and forearm IMU). gyro in rad/s, acc in m/s^2
% Assums right-handed coordinate frame with 
% +x points distally parallel with the bone, +y point forward in the 'I-pose' if
% IMU attached to the right side. 

% Author: Howard Chen, PhD
% Affiliation: University of Alabama in Huntsville

addpath(genpath('./filters')); 
addpath(genpath('./alignment'));
addpath(genpath('./kinematics')); 

% parameters
freq = 120; 
gyr_noise = 0.005; 
con_acc_mekf_acc = .01;
con_acc_rts_acc = 0.02; 
con_acc_mekf_dof = 0.05;  
con_dof_mekf_dof = 0.05;  
con_acc_rts_dof = 0.01;  
con_dof_rts_dof = 0.02;  

%% solve for translation alignment
[rUA2,rFA] = lev_calc(ua.acc,ua.gyr,fa.acc,fa.gyr,1/freq);

%% solve for rotational alignment
[~,~,ua.quat_s,fa.quat_s] = mekf_acc_s(ua,fa,freq,gyr_noise, con_acc_rts_acc, rUA2,rFA); 
[q1_imu, q2_imu,~, ~,er] = laidig_elbow_align(ua.quat_s,ua.gyr(3:end-2,:),fa.quat_s,fa.gyr(3:end-2,:),100);

%% run filter
% multiplicative kalman filter with linear acceleration constraint
el_mekf_acc = mekf_acc(ua,fa,freq,gyr_noise, con_acc_mekf_acc, rUA2,rFA,q1_imu,q2_imu); 
eul_el_mekf_acc = rad2deg(quatToEuler(el_mekf_acc)); 

% multiplicative kalman soother with linear acceleration constraint
el_rts_acc = mekf_acc_s(ua,fa,freq,gyr_noise, con_acc_rts_acc, rUA2,rFA,q1_imu,q2_imu); 
eul_el_rts_acc = rad2deg(quatToEuler(el_rts_acc));

% multiplicative kalman filter with linear acceleration and dof constraint
el_mekf_dof = mekf_elbow_acc(ua,fa,freq,gyr_noise, con_acc_mekf_dof, con_dof_mekf_dof, rUA2,rFA,q1_imu,q2_imu); 
eul_el_mekf_dof = rad2deg(quatToEuler(el_mekf_dof));

% multiplicative kalman smoother with linear acceleration and dof constraint
el_rts_dof = mekf_elbow_acc_s(ua,fa,freq,gyr_noise, con_acc_rts_dof, con_dof_rts_dof, rUA2,rFA,q1_imu,q2_imu); 
eul_el_rts_dof = rad2deg(quatToEuler(el_rts_dof));

%% Plot data
subplot(3,1,1);
plot([eul_el_mekf_acc(:,1),eul_el_rts_acc(:,1),eul_el_mekf_dof(:,1),eul_el_rts_dof(:,1)]); 
ylabel('flexion');
legend('mekf-acc','rts-acc','mekf-dof','rts-dof'); 

subplot(3,1,2);
plot([eul_el_mekf_acc(:,2),eul_el_rts_acc(:,2),eul_el_mekf_dof(:,2),eul_el_rts_dof(:,2)]); 
ylabel('abduction'); 

subplot(3,1,3);
plot([eul_el_mekf_acc(:,3),eul_el_rts_acc(:,3),eul_el_mekf_dof(:,3),eul_el_rts_dof(:,3)]); 
ylabel('rotation');

