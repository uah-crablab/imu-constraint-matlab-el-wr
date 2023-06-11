function [q_wr,faQuat,haQuat] = mekf_wrist_acc(fa,ha,freq,gyrNoise, conNoise, dofNoise, rUA,rFA,q1,q2)
%Multiplicative extended kalman filter with linear acceleration and dof
%constraint for wrist joint 
%   Author: Howard Chen, PhD  
%   University of Alabama in Huntsville

% Inputs: 
% fa.gyr- gyroscope measurement of forearm segment (rad/s)
% fa.acc- accelerometer measurement of forearm segment (m/s^2)
% ha.gyr- gyroscope measurement of hand segment (rad/s)
% ha.acc- accelerometer measurement of hand segment (m/s^2)
% freq: sampling frequency (Hz)
% gyrNoise: gyroscope noise (stdev, rad/s)
% conNoise: linear constraint noise
% dofNoise: rotational constraint noise
% rFA: distance from joint center to forearm IMU
% rHA: distance from joint center to hand IMU
% q1: rotation of forearm IMU to align with joint center
% q2: rotation of hand IMU to align with joint ceenter 
%
% Outputs:
% q_wr: wrist orientation from filter (hamiltonian quaternion)
% faQuat: orienation of the forearm segment from filter (hamiltonian
% quaternion)
% haQuat: orientation of the hand segment from filter (hamiltonian
% quaternion)

dT = 1/freq; 

lUA = quatToDCM(q1)'*rUA;
fa.gyr = fa.gyr*quatToDCM(q1); 
fa.acc = fa.acc*quatToDCM(q1); 

lFA = quatToDCM(q2)'*rFA; 
ha.gyr = ha.gyr*quatToDCM(q2); 
ha.acc = ha.acc*quatToDCM(q2); 

gyr_dot1 = omDot(fa.gyr,dT); 
KK1 = K(fa.gyr,gyr_dot1); 

gyr_dot2 = omDot(ha.gyr,dT); 
KK2 = K(ha.gyr,gyr_dot2); 


n = size(gyr_dot1,1); 
I = eye(3); O = zeros(3); 
Q = eye(6).*gyrNoise.^2; 

R = [I.*conNoise^2 zeros(3,1);...
     zeros(1,3) dofNoise.^2]; %

faQuat = zeros(n,4); 
haQuat = zeros(n,4);

O = zeros(3); 
G = eye(6).*dT; 

P = Q; 
X = [1;0;0;0;1;0;0;0]; 

for i = 1:n
    % Nominal state (4.54a)
    gyr1 = fa.gyr(i+2,:); 
    gyr2 = ha.gyr(i+2,:); 
    
    AA = (eye(4) + 0.5.*[0 -gyr1; gyr1' -skew(gyr1)].*dT);
    BB = (eye(4) + 0.5.*[0 -gyr2; gyr2' -skew(gyr2)].*dT); 

    X = [AA zeros(4); zeros(4) BB]*X; 

    % quaterion to DCM (body to earth)
    C1 = (X(1)^2-X(2:4)'*X(2:4)).*eye(3) + 2*X(2:4)*X(2:4)' + 2*X(1)*skew(X(2:4));
    C2 = (X(5)^2-X(6:8)'*X(6:8)).*eye(3) + 2*X(6:8)*X(6:8)' + 2*X(5)*skew(X(6:8));

    F = [I-skew(gyr1).*dT O; O I-skew(gyr2).*dT];

    P = F*P*F'+G*Q*G'; 

    a1 = fa.acc(i+2,:)'-KK1(i*3-2:i*3,:)*lUA;
    a2 = ha.acc(i+2,:)'-KK2(i*3-2:i*3,:)*lFA;

    H = [C1*skew(a1) -C2*skew(a2)]; 
    
    J1 = transpose(C1)*C2; 
    J2 = transpose(C2)*C1; 
    
    H(4,:) = [[0 0 1]*J1*[0;0;1], 0, -[1 0 0]*J1*[0;0;1], -[0 1 0]*J2*[0;1;0], [1 0 0]*J2*[0;1;0], 0]; 


    k =  P*H'*(H*P*H'+R)^-1;
    x = k*([(C1*a1 - C2*a2);- [0 1 0]*transpose(C1)*C2*[0; 0; 1]]); % wrist


    P=(eye(6)-k*H)*P;

    X(1:4) = quatNormalize(quatMultiply(X(1:4)',[1 0.5.*x(1:3)']))'; 
    X(5:8) = quatNormalize(quatMultiply(X(5:8)',[1 0.5.*x(4:6)']))'; 

    faQuat(i,:) = X(1:4)'; 
    haQuat(i,:) = X(5:8)'; 
    
end
    q_wr = quatMultiply(quatConj(faQuat),haQuat); 
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
