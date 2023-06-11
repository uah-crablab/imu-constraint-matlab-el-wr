function [R] = quatToDCM(q)
%quatToDCM converts from quaternion vector to DCM
% input: quaternion vector (nx4) first element is real
% output: Rotation matrix (3x3xn)
%
% Author: Howard Chen

    R = zeros(3,3,size(q,1));
    
    R(1,1,:)= q(:,1).^2+q(:,2).^2-q(:,3).^2-q(:,4).^2;
    R(1,2,:)=-2.*q(:,1).*q(:,4)+2.*q(:,2).*q(:,3);
    R(1,3,:)= 2.*q(:,1).*q(:,3)+2.*q(:,2).*q(:,4);
    R(2,1,:)= 2.*q(:,1).*q(:,4)+2.*q(:,2).*q(:,3);
    R(2,2,:)= q(:,1).^2-q(:,2).^2+q(:,3).^2-q(:,4).^2;
    R(2,3,:)=-2.*q(:,1).*q(:,2)+2.*q(:,3).*q(:,4);
    R(3,1,:)=-2.*q(:,1).*q(:,3)+2.*q(:,2).*q(:,4);
    R(3,2,:)=2.*q(:,1).*q(:,2)+2.*q(:,3).*q(:,4);
    R(3,3,:)=q(:,1).^2-q(:,2).^2-q(:,3).^2+q(:,4).^2;
end


