% Example script for solving for translation and rotational alignment, and using these
% paramenters to run mekf-acc, rts-acc, mekf-dof, and rts-dof for the wrist joint 
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
con_acc_rts_dof = 0.04;  
con_dof_rts_dof = 0.08;  

%% solve for translational alignment parameters
[rFA2,rHA] = lev_calc(fa.acc,fa.gyr,ha.acc,ha.gyr,1/freq);

%% solve for rotational alignment parameters
[~,~,fa.quat_s,ha.quat_s] = mekf_acc_s(fa,ha,freq,gyr_noise, con_acc_rts_acc, rFA2,rHA); 
[q1_imu,q2_imu,er] = wristAlign(fa.quat_s,ha.quat_s,50);

%% run filter
% multiplicative kalman filter with linear acceleration constraint
wr_mekf_acc = mekf_acc(fa,ha,freq,gyr_noise, con_acc_mekf_acc, rFA2,rHA,q1_imu,q2_imu); 
eul_wr_mekf_acc = rad2deg(quatToXYZ(wr_mekf_acc)); 

% multiplicative kalman soother with linear acceleration constraint
wr_rts_acc = mekf_acc_s(fa,ha,freq,gyr_noise, con_acc_rts_acc, rFA2,rHA,q1_imu,q2_imu); 
eul_wr_rts_acc = rad2deg(quatToXYZ(wr_rts_acc));

% multiplicative kalman filter with linear acceleration and dof constraint
wr_mekf_dof = mekf_wrist_acc(fa,ha,freq,gyr_noise, con_acc_mekf_dof, con_dof_mekf_dof, rFA2,rHA,q1_imu,q2_imu); 
eul_wr_mekf_dof = rad2deg(quatToXYZ(wr_mekf_dof));

% multiplicative kalman smoother with linear acceleration and dof constraint
wr_rts_dof = mekf_wrist_acc_s(fa,ha,freq,gyr_noise, con_acc_rts_dof, con_dof_rts_dof, rFA2,rHA,q1_imu,q2_imu); 
eul_wr_rts_dof = rad2deg(quatToXYZ(wr_rts_dof));

%% Plotting
subplot(3,1,1);
plot([eul_wr_mekf_acc(:,1),eul_wr_rts_acc(:,1),eul_wr_mekf_dof(:,1),eul_wr_rts_dof(:,1)]); 
ylabel('rotation');
legend('mekf-acc','rts-acc','mekf-dof','rts-dof'); 

subplot(3,1,2);
plot([eul_wr_mekf_acc(:,2),eul_wr_rts_acc(:,2),eul_wr_mekf_dof(:,2),eul_wr_rts_dof(:,2)]); 
ylabel('flexion'); 

subplot(3,1,3);
plot([eul_wr_mekf_acc(:,3),eul_wr_rts_acc(:,3),eul_wr_mekf_dof(:,3),eul_wr_rts_dof(:,3)]); 
ylabel('deviation');