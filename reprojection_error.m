% This function calculates the reprojection error between projected 3d
% world points (using A, R, T matrices) and actual 2d image points

% c_2d = 2xN
% w_3d = 3xN
% v = 2Nx1

function [v] = reprojection_error(A, R, T , c_2d , w_3d)
    
    v = zeros(2*size(w_3d,2), 1);
    for i=1:size(w_3d,2)
        proj_2d = R*(w_3d(:,i)) + T;
        proj_2d = A * proj_2d;
        proj_2d = proj_2d(1:2, :) ./ repmat(proj_2d(3, :), [2 1]);
        v(2*(i-1)+1:2*(i-1)+2) = proj_2d - c_2d(:,i);
    end
    
end