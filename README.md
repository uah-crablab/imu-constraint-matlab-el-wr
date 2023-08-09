# imu-constraint-matlab-el-wr
Code pertaining to the Journal Paper: 

Chen, H.; Schall, M.C., Jr.; Martin, S.M.; Fethke, N.B. Drift-Free Joint Angle Calculation Using Inertial Measurement Units without Magnetometers: An Exploration of Sensor Fusion Methods for the Elbow and Wrist. Sensors 2023, 23, 7053. https://doi.org/10.3390/s23167053

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
