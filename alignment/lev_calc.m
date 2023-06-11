function [r_prox,r_dist] = lev_calc(acc_prox, gyr_prox,acc_dist,gyr_dist,dT)
% Author: Howard Chen
% Affiliation: University of Alabama in Huntsville
% 
% implementation of the following translational alignment from:
% Seel, T., Raisch, J., & Schauer, T. (2014). IMU-based joint angle 
% measurement for gait analysis. Sensors, 14(4), 6891-6909.
%
% inputs: 
% acc_prox: acceleration measurement of IMU on proximal segment (m/s^2);
% gyr_prox: gyroscope measurement of IMU on proximal segment (rad/s);
% acc_dist: acceleration measurement of IMU on distal segment (m/s^2);
% gyr_dist: gyroscope measurement of IMU on dist segment (rad/s);
% dT: samplig period (seconds)
%
% outputs:
% r_prox: distance from joint center to IMU on proximal segment (meters)
% r_dost: distance from joint center to IMU on distal segment (meter)

max_iter = 20; 
gyr_dot1 = omDot(gyr_prox,dT); 
KK1 = K(gyr_prox,gyr_dot1); 

gyr_dot2 = omDot(gyr_dist,dT); 
KK2 = K(gyr_dist,gyr_dot2); 

n = size(gyr_dot1,1); 
e = zeros(n,1); 
J = zeros(n,6); 

r =zeros(6,1);

for ii = 1:max_iter
    for i = 1:n
        a1 = acc_prox(i+2,:)'- KK1(i*3-2:i*3,:)*r(1:3);
        a2 = acc_dist(i+2,:)'- KK2(i*3-2:i*3,:)*r(4:6);
        e1 = sqrt(a1(1)^2 + a1(2)^2 + a1(3)^2);
        e2 = sqrt(a2(1)^2 + a2(2)^2 + a2(3)^2);

        e(i) = e1 - e2;    
        J(i,1:3) = -(transpose(a1)*KK1(i*3-2:i*3,:))/e1; 
        J(i,4:6) =  (transpose(a2)*KK2(i*3-2:i*3,:))/e2; 

    end
    r = r - pinv(J)*e;
    er = e'*e
end

r_prox = r(1:3);
r_dist = r(4:6);

end

function out = omDot(gyro,dT)
% calculation of numerical derivative of gyroscope
n = length(gyro); 
out = zeros(n-4,3); 
for i = 3:n-2
    out(i-2,:) = (gyro(i-2,:)-8*gyro(i-1,:)+8*gyro(i+1,:)-gyro(i+2,:))./(12*dT); 
end

end

function out = K(gyr,gyr_dot)
% lever arm equation in matrix form 
gyr(1:2,:) = []; %remove first two indices since there is no correspondance to omega dot
n = size(gyr_dot,1); 
out = zeros(3*n,3); 
for i = 1:n
    out(i*3-2:i*3,:) = [-gyr(i,2)^2-gyr(i,3)^2, gyr(i,1)*gyr(i,2)-gyr_dot(i,3), gyr_dot(i,2)+gyr(i,1)*gyr(i,3);...
                         gyr_dot(i,3)+gyr(i,1)*gyr(i,2), -gyr(i,1)^2-gyr(i,3)^2, gyr(i,2)*gyr(i,3)-gyr_dot(i,1);...
                         gyr(i,1)*gyr(i,3)-gyr_dot(i,2) gyr_dot(i,1)+gyr(i,2)*gyr(i,3) -gyr(i,1)^2-gyr(i,2)^2];

end
end
