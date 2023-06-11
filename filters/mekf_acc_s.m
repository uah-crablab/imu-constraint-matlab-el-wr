function [q_s,q,proxQuat_s,distQuat_s] = mekf_acc_s(prox,dist,freq,gyrNoise, conNoise, rProx,rDist,q1,q2)
%Multiplicative rauch-tung-striebel smoother with linear acceleration
%constraint
%   Author: Howard Chen, PhD  
%   University of Alabama in Huntsville

% Inputs: 
% prox.gyr: gyroscope measurement of proximal segment (rad/s)
% prox.acc: accelerometer measurement of proximal segment (m/s^2)
% dist.gyr: gyroscope measurement of distal segment (rad/s)
% dist.acc: accelerometer measurement of distal segment (m/s^2)
% freq: sampling frequency (Hz)
% gyrNoise: gyroscope noise (stdev, rad/s)
% conNoise: constraint noise
% rProx: distance from joint center to proximal IMU
% rDist: distance from joint center to distal IMU
% q1: rotation of proximal IMU to align with joint center
% q2: rotation of distal IMU to align with joint ceenter 
%
% Outputs:
% q_s- relative orientation of the smoother (hamiltonian quaternion)
% q- relative orientation from the filter (hamiltonian quaternion)
% proxQuat_s- orienation of the proximal segment from smoother (hamiltonian
% quaternion)
% distQuat_s- orientation of the distal segment from smoother (hamiltonian
% quaternion)

if nargin == 7
    q1 = [1 0 0 0];
    q2 = [1 0 0 0]; 
end
dT = 1/freq; 

lUA = quatToDCM(q1)'*rProx;
prox.gyr = prox.gyr*quatToDCM(q1); 
prox.acc = prox.acc*quatToDCM(q1); 

lFA = quatToDCM(q2)'*rDist; 
dist.gyr = dist.gyr*quatToDCM(q2); 
dist.acc = dist.acc*quatToDCM(q2); 

gyr_dot1 = omDot(prox.gyr,dT); 
KK1 = K(prox.gyr,gyr_dot1); 

gyr_dot2 = omDot(dist.gyr,dT); 
KK2 = K(dist.gyr,gyr_dot2); 

n = size(gyr_dot1,1); 
I = eye(3); 
Q = eye(6).*gyrNoise.^2; 

R = I.*conNoise^2; %

uaQuat = zeros(n,4); 
uaQuat(1,:) = [1 0 0 0];
faQuat = zeros(n,4);
faQuat(1,:) = [1 0 0 0];

O = zeros(3); 
G = eye(6).*dT; 

P = zeros(6,6,n);
P(:,:,1) = Q; 
P_k = P; 
X = [ones(1,n);zeros(3,n);ones(1,n);zeros(3,n)];
x = zeros(6,n); 
for i = 2:n
    % Nominal state (4.54a)
    gyr1 = prox.gyr(i+2,:); 
    gyr2 = dist.gyr(i+2,:); 
    
    AA = (eye(4) + 0.5.*[0 -gyr1; gyr1' -skew(gyr1)].*dT);
    BB = (eye(4) + 0.5.*[0 -gyr2; gyr2' -skew(gyr2)].*dT); 
    
    X(:,i) = [AA zeros(4); zeros(4) BB]*X(:,i-1); 

    % quaterion to DCM (body to earth)
    C1 = (X(1,i)^2-X(2:4,i)'*X(2:4,i)).*eye(3) + 2*X(2:4,i)*X(2:4,i)' + 2*X(1,i)*skew(X(2:4,i));
    C2 = (X(5,i)^2-X(6:8,i)'*X(6:8,i)).*eye(3) + 2*X(6:8,i)*X(6:8,i)' + 2*X(5,i)*skew(X(6:8,i));

    F = [I-skew(gyr1).*dT O; O I-skew(gyr2).*dT];
    P(:,:,i) = F*P_k(:,:,i-1)*F'+G*Q*G'; 

    a1 = prox.acc(i+2,:)'-KK1(i*3-2:i*3,:)*lUA;
    a2 = dist.acc(i+2,:)'-KK2(i*3-2:i*3,:)*lFA;

    H = [C1*skew(a1) -C2*skew(a2)];  

    k =  P(:,:,i)*H'*(H*P(:,:,i)*H'+R)^-1;
    x(:,i) = k*((C1*a1 - C2*a2)); 

    P_k(:,:,i)=(eye(6)-k*H)*P(:,:,i) ;

    X(1:4,i) = quatNormalize(quatMultiply(X(1:4,i)',[1 0.5.*x(1:3,i)']))'; 
    X(5:8,i) = quatNormalize(quatMultiply(X(5:8,i)',[1 0.5.*x(4:6,i)']))'; 
    uaQuat(i,:) = X(1:4,i)'; 
    faQuat(i,:) = X(5:8,i)'; 
    
end
Xs = [ones(1,n);zeros(3,n);ones(1,n);zeros(3,n)];
xs = x(:,n); 
Xs(:,n) = X(:,n); 
proxQuat_s = uaQuat;
distQuat_s = faQuat; 
Ps = P(:,:,n); 
for i = n-1:-1:1
    gyr1 = prox.gyr(i+2,:); 
    gyr2 = dist.gyr(i+2,:);
    
    F = [I-skew(gyr1).*dT O; O I-skew(gyr2).*dT];
    
    C = P_k(:,:,i)*F'*P(:,:,i+1)^-1;
 
    xs = x(:,i) + C*xs;
    Ps = P_k(:,:,i) + C*(Ps-P(:,:,i+1))*C'; 
    
    % Correct nominal state
    Xs(1:4,i) = (eye(4) + 0.5.*[0 -xs(1:3)'; xs(1:3) -skew(xs(1:3))])*X(1:4,i+1);
    Xs(1:4,i) = Xs(1:4,i)./norm(Xs(1:4,i)); 
    
    Xs(1:4,i) = quatNormalize(quatMultiply(X(1:4,i)',[1 0.5.*xs(1:3)']))'; 
    Xs(5:8,i) = quatNormalize(quatMultiply(X(5:8,i)',[1 0.5.*xs(4:6)']))'; 
    proxQuat_s(i,:) = Xs(1:4,i)'; 
    distQuat_s(i,:) = Xs(5:8,i)'; 
    
end
    q = quatMultiply(quatConj(uaQuat),faQuat); 
    q_s = quatMultiply(quatConj(proxQuat_s),distQuat_s); 

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
