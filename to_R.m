
%% Rodrigue's formula 
function [R] = to_R(w)

    omega = [0 -w(3) w(2);  
             w(3) 0  -w(1);
             -w(2) w(1) 0]; 
    theta = norm(w); 

    R =  eye(3)  +  (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*omega^2;
    
end