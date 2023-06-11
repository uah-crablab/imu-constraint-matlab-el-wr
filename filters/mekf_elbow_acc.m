function [q,uaQuat,faQuat] = mekf_elbow_acc(ua,fa,freq,gyrNoise, conNoise, dofNoise, rUA,rFA,q1,q2)
%Multiplicative extended kalman filter with linear acceleration and dof
%constraint for elbow joint 
%   Author: Howard Chen, PhD  
%   University of Alabama in Huntsville

% Inputs: 
% ua.gyr- gyroscope measurement of upper arm segment (rad/s)
% ua.acc- accelerometer measurement of upper arm segment (m/s^2)
% fa.gyr- gyroscope measurement of forearm segment (rad/s)
% fa.acc- accelerometer measurement of forearm segment (m/s^2)
% freq: sampling frequency (Hz)
% gyrNoise: gyroscope noise (stdev, rad/s)
% conNoise: linear constraint noise
% dofNoise: rotational constraint noise
% rUA: distance from joint center to upper arm IMU
% rFA: distance from joint center to forearm IMU
% q1: rotation of upper arm IMU to align with joint center
% q2: rotation of forearm IMU to align with joint ceenter 
%
% Outputs:
% q_el: elbow orientation from filter (hamiltonian quaternion)
% uaQuat: orienation of the upper arm segment from filter (hamiltonian
% quaternion)
% faQuat: orientation of the forearm segment from filter (hamiltonian
% quaternion)

dT = 1/freq; 

lUA = quatToDCM(q1)'*rUA;
ua.gyr = ua.gyr*quatToDCM(q1); 
ua.acc = ua.acc*quatToDCM(q1); 

lFA = quatToDCM(q2)'*rFA; 
fa.gyr = fa.gyr*quatToDCM(q2); 
fa.acc = fa.acc*quatToDCM(q2); 

gyr_dot1 = omDot(ua.gyr,dT); 
KK1 = K(ua.gyr,gyr_dot1); 

gyr_dot2 = omDot(fa.gyr,dT); 
KK2 = K(fa.gyr,gyr_dot2); 

n = size(gyr_dot1,1); 
I = eye(3);
Q = eye(6).*gyrNoise.^2; 

R = [I.*conNoise^2 zeros(3,1);...
     zeros(1,3) dofNoise.^2]; %

uaQuat = zeros(n,4); 
faQuat = zeros(n,4);

O = zeros(3); 
G = eye(6).*dT; 

P = Q; 
X = [1;0;0;0;1;0;0;0]; 

g = [0;0;9.81]; 
for i = 1:n
    % Nominal state (4.54a)
    gyr1 = ua.gyr(i+2,:); 
    gyr2 = fa.gyr(i+2,:); 
    
    AA = (eye(4) + 0.5.*[0 -gyr1; gyr1' -skew(gyr1)].*dT);
    BB = (eye(4) + 0.5.*[0 -gyr2; gyr2' -skew(gyr2)].*dT); 

    X = [AA zeros(4); zeros(4) BB]*X; 

    % quaterion to DCM (body to earth)
    C1 = (X(1)^2-X(2:4)'*X(2:4)).*eye(3) + 2*X(2:4)*X(2:4)' + 2*X(1)*skew(X(2:4));
    C2 = (X(5)^2-X(6:8)'*X(6:8)).*eye(3) + 2*X(6:8)*X(6:8)' + 2*X(5)*skew(X(6:8));

    F = [I-skew(gyr1).*dT O; O I-skew(gyr2).*dT];
    P = F*P*F'+G*Q*G'; 

    a1 = ua.acc(i+2,:)'-KK1(i*3-2:i*3,:)*lUA;
    a2 = fa.acc(i+2,:)'-KK2(i*3-2:i*3,:)*lFA;

    H = [C1*skew(a1) -C2*skew(a2)]; 
    
    J1 = transpose(C1)*C2; 
    J2 = transpose(C2)*C1; 
    
    H(4,:) = [-[0 1 0]*J1*[1;0;0], [1 0 0]*J1*[1;0;0], 0, 0, -[0 0 1]*J2*[0;0;1], [0 1 0]*J2*[0;0;1]]; 

    k =  P*H'*(H*P*H'+R)^-1;
    x = k*([(C1*a1 - C2*a2);-[0 0 1]*transpose(C1)*C2*[1; 0; 0]]); 

    P=(eye(6)-k*H)*P;

    X(1:4) = quatNormalize(quatMultiply(X(1:4)',[1 0.5.*x(1:3)']))'; 
    X(5:8) = quatNormalize(quatMultiply(X(5:8)',[1 0.5.*x(4:6)']))'; 

    uaQuat(i,:) = X(1:4)'; 
    faQuat(i,:) = X(5:8)'; 
    
end
    q = quatMultiply(quatConj(uaQuat),faQuat); 
end


function [ out ] = skew( a )
    out = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
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
