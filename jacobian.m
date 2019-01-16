% This function returns the Jacobian of the projection of 3d points (in world 
% coordinate system) to the image pixels.

function [dm_dp] = jacobian(A, w, w_3d, c_3d)
    
    % w_3d = 3xN
    % c_3d = 3xN
    % dm_dp = 2Nx6
    
    image_3d = A*c_3d;
    dm_dp = zeros(2*size(w_3d,2), 6);
    R = to_R(w); %calling rodrigues formula 
    omega = skew_symmetric(w); 
    I = eye(3);
    dR_dr1 = (1/norm(w)^2)*(w(1)*omega + ...
        skew_symmetric(cross(w,(I - R)*I(:,1))))*R ;
    dR_dr2 = (1/norm(w)^2)*(w(2)*omega + ...
        skew_symmetric(cross(w,(I - R)*I(:,2))))*R ;
    dR_dr3 = (1/norm(w)^2)*(w(3)*omega + ...
        skew_symmetric(cross(w,(I - R)*I(:,3))))*R ;

    for i=1:size(w_3d,2)
        dM_dp = [dR_dr1*w_3d(:,i) dR_dr2*w_3d(:,i) dR_dr3*w_3d(:,i) eye(3)];
        dmt_dM = A;
        
        dm_dmt = [1.0/image_3d(3,i) 0  -image_3d(1,i)/image_3d(3,i)^2;
                   0       1.0/image_3d(3,i) -image_3d(2,i)/image_3d(3,i)^2];
 
        dm_dp(2*(i-1)+1:2*(i-1)+2, :) = dm_dmt * dmt_dM * dM_dp;
    end

end

function vx = skew_symmetric(w)
    vx =[0 -w(3) w(2);  
        w(3) 0  -w(1);
        -w(2) w(1) 0];
end