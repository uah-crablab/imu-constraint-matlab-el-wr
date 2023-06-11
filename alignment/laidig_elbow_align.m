function [q1, q2, j1, j2,err] = laidig_elbow_align(ua_quat, ua_gyr, fa_quat, fa_gyr,max_iter)
% Author: Howard Chen
% Affiliation: University of Alabama in Huntsville
% 
% implementation of the following:
% Laidig, D., Müller, P., & Seel, T. (2017). Automatic anatomical calibration 
% for IMU-based elbow angle measurement in disturbed magnetic fields. 
% Current directions in biomedical engineering, 3(2), 167-170.
%
% Inputs:
% ua_quat: IMU orientation of the sensor attached to upper arm
% ua_gyr: gyro measurements of the IMU attached to upper arm
% fa_quat: IMU orientation of the sensor attached to forearm
% fa_gyr: gyro measurements of the IMU attached to forearm
% max_itre: maximum interation of the solver 
%
% output: 
% q1: rotation of the upper arm IMU
% q2: rotation of the forearm IMU
% j1: upper arm axis solved by the algorithm
% j2: forearm axis solved by the algorithm

x =  .001.*ones(4,1); 
w_rel = (quatRotate(ua_quat,ua_gyr) - quatRotate(fa_quat,fa_gyr))'; 

use_J1_1 = true; 
use_J2_1 = true; 

JJ = J11(ua_quat,ua_gyr,fa_quat,fa_gyr,x); 
[ee, j1, j2] = e11(ua_quat,fa_quat,w_rel,x); 

for j = 1:max_iter
    if abs(sin(x(1))) < 0.5
        use_J1_1 = ~use_J1_1; 
    end
    
    if abs(sin(x(3))) < 0.5
        use_J2_1 = ~use_J2_1; 
    end
    
    if use_J1_1 == true && use_J2_1 == true 
        JJ = J11(ua_quat,ua_gyr,fa_quat,fa_gyr,x); 
        [ee, j1, j2] = e11(ua_quat,fa_quat,w_rel,x);
    end
    
    if use_J1_1 == true && use_J2_1 == false
        JJ = J12(ua_quat,ua_gyr,fa_quat,fa_gyr,x); 
        [ee, j1, j2] = e12(ua_quat,fa_quat,w_rel,x); 
    end
    
    if use_J1_1 == false && use_J2_1 == true
        JJ = J21(ua_quat,ua_gyr,fa_quat,fa_gyr,x); 
        [ee, j1, j2] = e21(ua_quat,fa_quat,w_rel,x); 
    end
    
        
    if use_J1_1 == false && use_J2_1 == false
        JJ = J22(ua_quat,ua_gyr,fa_quat,fa_gyr,x); 
        [ee, j1, j2] = e22(ua_quat,fa_quat,w_rel,x); 
    end
    
    err = ee'*ee
    
    x_prev = x; 
    x = x - pinv(JJ)*ee; 
    
    if j ~=1 && all(abs((x-x_prev)./x) < .001)
        disp('solution converged')
        break
    end
    
    if j == max_iter
        disp('max iter has been reached'); 
    end
end

e3 = [0 0 1]; 

alpha1 = acos(dot(e3,j1)); 
u1 = cross(e3, j1); 
alpha2 = acos(dot(e3,j2)); 
u2 = cross(e3, j2); 

q1 = [cos(alpha1/2) sin(alpha1/2).*u1./vecnorm(u1)]; 
q2 = [cos(alpha2/2) sin(alpha2/2).*u2./vecnorm(u2)]; 

q_el = quatMultiply(quatConj(quatMultiply(ua_quat,q1)),quatMultiply(fa_quat,q2)); 

eul = quatToEuler(q_el);

% q1 = quatMultiply(q1,quatZ(eul(zero_ind,1))); 
% q2 = quatMultiply(q2,quatX(-eul(zero_ind,3))); 

q1 = quatMultiply(q1,quatZ(min(eul(:,1)))); 
q2 = quatMultiply(q2,quatX(-median(eul(:,3)))); 

end

function [out, j1_1, j2_1] = e11(ua_quat, fa_quat, wrel, x)

    j1_1 = [sin(x(1))*cos(x(2)); sin(x(1))*sin(x(2)); cos(x(1))]; 
    j2_1 = [sin(x(3))*cos(x(4)); sin(x(3))*sin(x(4)); cos(x(3))]; 
    
    jn = cross(quatRotate(ua_quat,j1_1),quatRotate(fa_quat,j2_1));
    out = dot(wrel, jn)'; 

end

function [out, j1_1, j2_2] = e12(ua_quat, fa_quat, wrel, x)

    j1_1 = [sin(x(1))*cos(x(2)); sin(x(1))*sin(x(2)); cos(x(1))]; 
    j2_2 = [cos(x(3)); sin(x(3))*sin(x(4)); sin(x(3))*cos(x(4))]; 
    
    jn = cross(quatRotate(ua_quat,j1_1),quatRotate(fa_quat,j2_2));
    out = dot(wrel, jn)'; 

end

function [out, j1_2, j2_1] = e21(ua_quat, fa_quat, wrel, x)

    j1_2 = [cos(x(1)); sin(x(1))*sin(x(2)); sin(x(1))*cos(x(2))]; 
    j2_1 = [sin(x(3))*cos(x(4)); sin(x(3))*sin(x(4)); cos(x(3))]; 
    
    jn = cross(quatRotate(ua_quat,j1_2),quatRotate(fa_quat,j2_1));
    out = dot(wrel, jn)'; 

end

function [out, j1_2, j2_2] = e22(ua_quat, fa_quat, wrel, x)

    j1_2 = [cos(x(1)); sin(x(1))*sin(x(2)); sin(x(1))*cos(x(2))]; 
    j2_2 = [cos(x(3)); sin(x(3))*sin(x(4)); sin(x(3))*cos(x(4))]; 

    jn = cross(quatRotate(ua_quat,j1_2),quatRotate(fa_quat,j2_2));
    out = dot(wrel, jn)'; 

end


function out = J11(ua_quat, ua_gyr,fa_quat,fa_gyr,x) 
    n = size(ua_quat,1); 
    out = zeros(n,4);
    x1 = x(1); 
    x2 = x(2);
    x3 = x(3);
    x4 = x(4); 
    
    for i = 1:n 
    X1 = ua_quat(i,1); 
    X2 = ua_quat(i,2); 
    X3 = ua_quat(i,3); 
    X4 = ua_quat(i,4); 
    gx = ua_gyr(i,1); 
    gy = ua_gyr(i,2); 
    gz = ua_gyr(i,3); 

    XX1 = fa_quat(i,1); 
    XX2 = fa_quat(i,2); 
    XX3 = fa_quat(i,3); 
    XX4 = fa_quat(i,4); 
    ggx = fa_gyr(i,1); 
    ggy = fa_gyr(i,2); 
    ggz = fa_gyr(i,3); 

    out(i,:) = [((sin(x1)*(2*X1*X2 - 2*X3*X4) + cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) + (sin(x1)*(2*X1*X3 + 2*X2*X4) - cos(x1)*cos(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)) - ((sin(x1)*(2*X1*X2 - 2*X3*X4) + cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)) + (sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*cos(x2)*(2*X1*X3 - 2*X2*X4) - cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*cos(x2)*(2*X1*X3 - 2*X2*X4) - cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) - (sin(x1)*(2*X1*X3 + 2*X2*X4) - cos(x1)*cos(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)), ((cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4) + sin(x1)*sin(x2)*(2*X1*X3 - 2*X2*X4))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) + (sin(x1)*sin(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) + ((cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4) + sin(x1)*sin(x2)*(2*X1*X3 - 2*X2*X4))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)) - (cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) + ((sin(x1)*sin(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)) + (cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)), ((sin(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)) + (sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX3 - 2*XX2*XX4) - cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) + (sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - cos(x3)*cos(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)) + ((sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX3 - 2*XX2*XX4) - cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) - (sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - cos(x3)*cos(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)), - ((cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4) + sin(x3)*sin(x4)*(2*XX1*XX3 - 2*XX2*XX4))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) + (sin(x3)*sin(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4) + sin(x3)*sin(x4)*(2*XX1*XX3 - 2*XX2*XX4))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)) - (cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x3)*sin(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)) + (cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2))]; 
    
    end
end

function out = J12(ua_quat, ua_gyr,fa_quat,fa_gyr,x) 
    n = size(ua_quat,1); 
    out = zeros(n,4);
    x1 = x(1); 
    x2 = x(2);
    x3 = x(3);
    x4 = x(4); 
    
    for i = 1:n 
    X1 = ua_quat(i,1); 
    X2 = ua_quat(i,2); 
    X3 = ua_quat(i,3); 
    X4 = ua_quat(i,4); 
    gx = ua_gyr(i,1); 
    gy = ua_gyr(i,2); 
    gz = ua_gyr(i,3); 

    XX1 = fa_quat(i,1); 
    XX2 = fa_quat(i,2); 
    XX3 = fa_quat(i,3); 
    XX4 = fa_quat(i,4); 
    ggx = fa_gyr(i,1); 
    ggy = fa_gyr(i,2); 
    ggz = fa_gyr(i,3); 

    out(i,:) = [((sin(x1)*(2*X1*X2 - 2*X3*X4) + cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) + (sin(x1)*(2*X1*X3 + 2*X2*X4) - cos(x1)*cos(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)) - ((sin(x1)*(2*X1*X2 - 2*X3*X4) + cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)) + (sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*cos(x2)*(2*X1*X3 - 2*X2*X4) - cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*cos(x2)*(2*X1*X3 - 2*X2*X4) - cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) - (sin(x1)*(2*X1*X3 + 2*X2*X4) - cos(x1)*cos(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)), ((cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4) + sin(x1)*sin(x2)*(2*X1*X3 - 2*X2*X4))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) + (sin(x1)*sin(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) + ((cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4) + sin(x1)*sin(x2)*(2*X1*X3 - 2*X2*X4))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)) - (cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) + ((sin(x1)*sin(x2)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)) + (cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 + 2*X2*X3))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)), - ((sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + cos(x3)*cos(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) + (sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - cos(x3)*cos(x4)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + cos(x3)*cos(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)) + (sin(x3)*(2*XX1*XX4 + 2*XX2*XX3) - cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - cos(x3)*cos(x4)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)) - (sin(x3)*(2*XX1*XX4 + 2*XX2*XX3) - cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)), ((cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + sin(x3)*sin(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)) + (sin(x3)*sin(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3) + sin(x3)*sin(x4)*(2*XX1*XX3 + 2*XX2*XX4))*(cos(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)) - (sin(x3)*sin(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3) + sin(x3)*sin(x4)*(2*XX1*XX3 + 2*XX2*XX4))*(sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x1)*(2*X1*X2 - 2*X3*X4) + cos(x2)*sin(x1)*(2*X1*X4 + 2*X2*X3)) + (cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + sin(x3)*sin(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x1)*(2*X1*X3 + 2*X2*X4) + cos(x2)*sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2))]; 
    end
end
function out = J21(ua_quat, ua_gyr,fa_quat,fa_gyr,x) 
    n = size(ua_quat,1); 
    out = zeros(n,4);
    x1 = x(1); 
    x2 = x(2);
    x3 = x(3);
    x4 = x(4); 
    
    for i = 1:n 
    X1 = ua_quat(i,1); 
    X2 = ua_quat(i,2); 
    X3 = ua_quat(i,3); 
    X4 = ua_quat(i,4); 
    gx = ua_gyr(i,1); 
    gy = ua_gyr(i,2); 
    gz = ua_gyr(i,3); 

    XX1 = fa_quat(i,1); 
    XX2 = fa_quat(i,2); 
    XX3 = fa_quat(i,3); 
    XX4 = fa_quat(i,4); 
    ggx = fa_gyr(i,1); 
    ggy = fa_gyr(i,2); 
    ggz = fa_gyr(i,3); 

    out(i,:) = [((sin(x1)*(2*X1*X3 - 2*X2*X4) + cos(x1)*cos(x2)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) + (sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - cos(x1)*cos(x2)*(2*X1*X3 + 2*X2*X4) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) + ((sin(x1)*(2*X1*X3 - 2*X2*X4) + cos(x1)*cos(x2)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)) + (sin(x1)*(2*X1*X4 + 2*X2*X3) - cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) + ((sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - cos(x1)*cos(x2)*(2*X1*X3 + 2*X2*X4) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)) - (sin(x1)*(2*X1*X4 + 2*X2*X3) - cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)), ((cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3) + sin(x1)*sin(x2)*(2*X1*X3 + 2*X2*X4))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)) + (cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) + sin(x1)*sin(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)) + ((cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3) + sin(x1)*sin(x2)*(2*X1*X3 + 2*X2*X4))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)) - (sin(x1)*sin(x2)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x4)*sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) + sin(x1)*sin(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)) + (sin(x1)*sin(x2)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4))*(sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x4)*sin(x3)*(2*XX1*XX4 + 2*XX2*XX3)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)), ((sin(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)) + (sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX3 - 2*XX2*XX4) - cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x3)*(2*XX1*XX2 - 2*XX3*XX4) + cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) + (sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - cos(x3)*cos(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)) + ((sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX3 - 2*XX2*XX4) - cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) - (sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - cos(x3)*cos(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)), - ((cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4) + sin(x3)*sin(x4)*(2*XX1*XX3 - 2*XX2*XX4))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) + (sin(x3)*sin(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4) + sin(x3)*sin(x4)*(2*XX1*XX3 - 2*XX2*XX4))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)) - (cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x3)*sin(x4)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)) + (cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - sin(x3)*sin(x4)*(2*XX1*XX4 + 2*XX2*XX3))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2))]; 
    
    end
end

function out = J22(ua_quat, ua_gyr,fa_quat,fa_gyr,x) 
    n = size(ua_quat,1); 
    out = zeros(n,4);
    x1 = x(1); 
    x2 = x(2);
    x3 = x(3);
    x4 = x(4); 
    
    for i = 1:n 
    X1 = ua_quat(i,1); 
    X2 = ua_quat(i,2); 
    X3 = ua_quat(i,3); 
    X4 = ua_quat(i,4); 
    gx = ua_gyr(i,1); 
    gy = ua_gyr(i,2); 
    gz = ua_gyr(i,3); 

    XX1 = fa_quat(i,1); 
    XX2 = fa_quat(i,2); 
    XX3 = fa_quat(i,3); 
    XX4 = fa_quat(i,4); 
    ggx = fa_gyr(i,1); 
    ggy = fa_gyr(i,2); 
    ggz = fa_gyr(i,3); 

    out(i,:) = [((sin(x1)*(2*X1*X3 - 2*X2*X4) + cos(x1)*cos(x2)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)) + (sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - cos(x1)*cos(x2)*(2*X1*X3 + 2*X2*X4) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) + ((sin(x1)*(2*X1*X3 - 2*X2*X4) + cos(x1)*cos(x2)*(X1^2 - X2^2 - X3^2 + X4^2) + cos(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)) + (sin(x1)*(2*X1*X4 + 2*X2*X3) - cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) + ((sin(x1)*(X1^2 + X2^2 - X3^2 - X4^2) - cos(x1)*cos(x2)*(2*X1*X3 + 2*X2*X4) + cos(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)) - (sin(x1)*(2*X1*X4 + 2*X2*X3) - cos(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) + cos(x1)*cos(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)), ((cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3) + sin(x1)*sin(x2)*(2*X1*X3 + 2*X2*X4))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)) + (cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) + sin(x1)*sin(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)) + ((cos(x2)*sin(x1)*(2*X1*X4 - 2*X2*X3) + sin(x1)*sin(x2)*(2*X1*X3 + 2*X2*X4))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)) - (sin(x1)*sin(x2)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) + cos(x4)*sin(x3)*(2*XX1*XX3 + 2*XX2*XX4) - sin(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((cos(x2)*sin(x1)*(X1^2 - X2^2 + X3^2 - X4^2) + sin(x1)*sin(x2)*(2*X1*X2 - 2*X3*X4))*(cos(x4)*sin(x3)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x3)*(2*XX1*XX3 - 2*XX2*XX4) + sin(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4)) + (sin(x1)*sin(x2)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x2)*sin(x1)*(2*X1*X2 + 2*X3*X4))*(cos(x3)*(2*XX1*XX4 + 2*XX2*XX3) + sin(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 - 2*XX3*XX4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)), - ((sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + cos(x3)*cos(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)) + (sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - cos(x3)*cos(x4)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((sin(x3)*(2*XX1*XX3 - 2*XX2*XX4) + cos(x3)*cos(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) + cos(x3)*sin(x4)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)) + (sin(x3)*(2*XX1*XX4 + 2*XX2*XX3) - cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((sin(x3)*(XX1^2 + XX2^2 - XX3^2 - XX4^2) - cos(x3)*cos(x4)*(2*XX1*XX3 + 2*XX2*XX4) + cos(x3)*sin(x4)*(2*XX1*XX4 - 2*XX2*XX3))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)) - (sin(x3)*(2*XX1*XX4 + 2*XX2*XX3) - cos(x3)*sin(x4)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + cos(x3)*cos(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2)), ((cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + sin(x3)*sin(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)) + (sin(x3)*sin(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)))*(gy*(2*X1*X4 - 2*X2*X3) - gz*(2*X1*X3 + 2*X2*X4) - ggy*(2*XX1*XX4 - 2*XX2*XX3) + ggz*(2*XX1*XX3 + 2*XX2*XX4) - gx*(X1^2 + X2^2 - X3^2 - X4^2) + ggx*(XX1^2 + XX2^2 - XX3^2 - XX4^2)) - ((cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3) + sin(x3)*sin(x4)*(2*XX1*XX3 + 2*XX2*XX4))*(cos(x2)*sin(x1)*(X1^2 - X2^2 - X3^2 + X4^2) - cos(x1)*(2*X1*X3 - 2*X2*X4) + sin(x1)*sin(x2)*(2*X1*X2 + 2*X3*X4)) - (sin(x3)*sin(x4)*(XX1^2 - XX2^2 - XX3^2 + XX4^2) - cos(x4)*sin(x3)*(2*XX1*XX2 + 2*XX3*XX4))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X4 + 2*X2*X3) - gz*(2*X1*X2 - 2*X3*X4) - ggx*(2*XX1*XX4 + 2*XX2*XX3) + ggz*(2*XX1*XX2 - 2*XX3*XX4) + gy*(X1^2 - X2^2 + X3^2 - X4^2) - ggy*(XX1^2 - XX2^2 + XX3^2 - XX4^2)) - ((cos(x4)*sin(x3)*(2*XX1*XX4 - 2*XX2*XX3) + sin(x3)*sin(x4)*(2*XX1*XX3 + 2*XX2*XX4))*(cos(x1)*(2*X1*X4 + 2*X2*X3) + sin(x1)*sin(x2)*(X1^2 - X2^2 + X3^2 - X4^2) - cos(x2)*sin(x1)*(2*X1*X2 - 2*X3*X4)) + (cos(x4)*sin(x3)*(XX1^2 - XX2^2 + XX3^2 - XX4^2) + sin(x3)*sin(x4)*(2*XX1*XX2 - 2*XX3*XX4))*(cos(x1)*(X1^2 + X2^2 - X3^2 - X4^2) + cos(x2)*sin(x1)*(2*X1*X3 + 2*X2*X4) - sin(x1)*sin(x2)*(2*X1*X4 - 2*X2*X3)))*(gx*(2*X1*X3 - 2*X2*X4) - gy*(2*X1*X2 + 2*X3*X4) - ggx*(2*XX1*XX3 - 2*XX2*XX4) + ggy*(2*XX1*XX2 + 2*XX3*XX4) - gz*(X1^2 - X2^2 - X3^2 + X4^2) + ggz*(XX1^2 - XX2^2 - XX3^2 + XX4^2))]; 
    
    end
end


function out = quatRotate(quat,vecs)
    transepose_flag = false; 
    if size(vecs,1) == 3
       vecs = vecs';
       transepose_flag = true; 
    end
    
    n = size(vecs,1); 
    out = quatMultiply(quatMultiply(quat,[zeros(n,1) vecs]),quatConj(quat)); 
    out(:,1) = []; 
    
    if transepose_flag == true
        out = out'; 
    end
end
% 
% % %% math
% I = eye(3); O = zeros(3); 
% syms XX1 XX2 XX3 XX4 X1 X2 X3 X4 x1 x2 x3 x4 gx gy gz ggx ggy ggz
% X = [X1;X2;X3;X4];
% XX = [XX1;XX2;XX3;XX4];
% g = [gx;gy;gz];
% gg = [ggx; ggy; ggz]; 
% C = (X(1)^2-transpose(X(2:4))*X(2:4)).*I + 2*X(2:4)*transpose(X(2:4)) + 2*X(1)*skew(X(2:4));
% CC = (XX(1)^2-transpose(XX(2:4))*XX(2:4)).*I + 2*XX(2:4)*transpose(XX(2:4)) + 2*XX(1)*skew(XX(2:4));
% syms c1 c2 c3 c4 c5 c6 c7 c8 c9
% syms cc1 cc2 cc3 cc4 cc5 cc6 cc7 cc8 cc9
% 
% % C = [c1 c2 c3; c4 c5 c6; c7 c8 c9];
% % CC = [cc1 cc2 cc3; cc4 cc5 cc6; cc7 cc8 cc9];
% 
% j1_1 = [sin(x1)*cos(x2); sin(x1)*sin(x2); cos(x1)]; 
% j2_1 = [sin(x3)*cos(x4); sin(x3)*sin(x4); cos(x3)]; 
% 
% j1_2 = [cos(x1); sin(x1)*sin(x2); sin(x1)*cos(x2)]; 
% j2_2 = [cos(x3); sin(x3)*sin(x4); sin(x3)*cos(x4)]; 
% 
% e11 = transpose(C*g-CC*gg)* cross(C*j1_1,CC*j2_1); 
% J11 = jacobian(e11,[x1 x2 x3 x4]); 
% 
% e12 = transpose(C*g-CC*gg)* cross(C*j1_1,CC*j2_2); 
% J12 = jacobian(e12,[x1 x2 x3 x4]); 
% 
% e21 = transpose(C*g-CC*gg)* cross(C*j1_2,CC*j2_1); 
% J21 = jacobian(e21,[x1 x2 x3 x4]); 
% 
% e22 = transpose(C*g-CC*gg)* cross(C*j1_2,CC*j2_2); 
% J22 = jacobian(e22,[x1 x2 x3 x4]); 
% 
