# imu-constraint-matlab-el-wr
Code pertaining to the manuscript 'Drift-free joint angle calculation using inertial measurement units without magnetometers: an exploration of sensor fusion methods for the elbow and wrist' 

# contains code for the following:
* `lev_calc.m` IMU translational alignment
* `laidig_elbow_align.m` IMU Rotational alignment for elbow
* `wristAlign.m` IMU Rotational alignment for wrist
* `mekf_acc.m` Multiplicative Extended Kalman Filter with linear acceleration constraint
* `mekf_acc_s.m` Multiplicative Rauch-Tung-Striebel Smoother with linear acceleration constraint
* `mekf_elbow_acc.m` Multiplicative Extended Kalman Filter with linear acceleration and elbow rotational constraint
* `mekf_elbow_acc_s.m` Multiplicative Rauch-Tung-Striebel Smoother with linear acceleration constraint and elbow rotational constraint
* `mekf_wrist_acc.m` Multiplicative Extended Kalman Filter with linear acceleration and elbow rotational constraint
* `mekf_wrist_acc_s.m` Multiplicative Rauch-Tung-Striebel Smoother with linear acceleration constraint and elbow rotational constraint
* `main_elbow.m` Example run script for elbow
* `main_wrist.m` Example run script for wrist