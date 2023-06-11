function [off1,off2,er] = wristAlign(faQuat,haQuat,max_iter)

n = size(faQuat,1);
ee = zeros(n,1); 
JJ = zeros(n,6);
ii = 0; 

off1 = [1;0;0;0]';
off2 = [1;0;0;0]'; 
x = [0 0 0 0 0 0]'; 

for jj = 1:max_iter
    q_fa = quatMultiply(faQuat,off1);
    q_ha = quatMultiply(haQuat,off2); 
    for i = 1:n
        ua = quatToDCM(q_fa(i,:));
        fa = quatToDCM(q_ha(i,:));
           
        ee(i) = [0 1 0]*transpose(ua)*fa*[0; 0; 1]; 

        J1 = transpose(ua)*fa; 
        J2 = transpose(fa)*ua; 
        JJ(i,:) = [[0 0 1]*J1*[0;0;1], 0, -[1 0 0]*J1*[0;0;1], -[0 1 0]*J2*[0;1;0], [1 0 0]*J2*[0;1;0], 0]; 

    end
    
    x_prev = x; 
    x = -.0001*JJ'*ee; 
    off1 = quatMultiply(off1,[1 0.5.*x(1:3)']);
    off2 = quatMultiply(off2,[1 0.5.*x(4:6)']);
    
    er = ee'*ee
    ii = ii+1;

    if jj ~=1 && all(abs((x-x_prev)./x) < .001)
        disp('solution converged')
        break
    end

    if jj == max_iter
        disp('max iter has been reached'); 
    end
    
end

end

% function [ out ] = skew( a )
% %SKEW Skew symmetry matrix
% %   skew(a)*b = cross(a,b)
%   out = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0]; 
% end
% %%
% syms R11 R12 R13 R21 R22 R23 R31 R32 R33
% R = [R11 R12 R13; R21 R22 R23; R31 R32 R33]; 
% 
% syms RR11 RR12 RR13 RR21 RR22 RR23 RR31 RR32 RR33
% RR = [RR11 RR12 RR13; RR21 RR22 RR23; RR31 RR32 RR33]; 
% syms x1 x2 x3 x4; 
% x = [x1 x2 x3 x4]; 
% ua =  R*rotZ(x(1))*rotY(x(2));
% fa =  RR*rotZ(x(3))*rotY(x(4));
% e = [0 1 0]*transpose(ua)*fa*[0; 0; 1]
