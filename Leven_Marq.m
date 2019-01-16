% This function applies the Levenberg Marquardt algorithm provided in the pdf 
% for non-linear optimization of the pose returned by ransac function. The
% function implements both with and without M-estimator

% matches : contains the indexes of matched descriptors of sift features
%                 database and sift features extracted from an image

% f_image : locations of all sift features extracted from an image

% sift_locations : 3d locations of all sift features of box extracted from
%                 the init_texture folder.

% Intrinsic MAtrix(3x3) : intrinsic camera matrix
% R_best_from_ransac(3x3) : best rotation matrix returned by ransac
% T_best_from_ransac =(3x1): best translation vector returned by ransac
% use_weights   : boolean variable indicating whether to use tukey
%                 estimator's weight matrix or not


function [p] = Leven_Marq(intrinsic_Matrix, R_best_from_ransac, T_best_from_ransac, matches, f_image , sift_locations, use_weights)
    
    threshold = 0.0001;
    iter = 1000;

    R = R_best_from_ransac;
    T = T_best_from_ransac;
    
    % rotationMatrixToVector returns a vector which represents the axis
    % corresponds to the rotation axis and magnitude corresponds to the
    % rotation angle in radians. 
    
    %The function uses the Rodrigues rotation matrix implementation
    w = rotationMatrixToVector(R); 
    i=0;
    p = [w' ; T];
    u = threshold + 1;
    lambda = 0.001;
    c_2d = f_image(1:2, matches(1,:));
    w_3d = sift_locations(1:3, matches(2,:));

    while (i <= iter && u > threshold)
        % J = 2Nx6
        % e = 2Nx1
        % new_e = 2Nx1

        % calculating Jacobian and reprojection error using matched 2d-3d 
        % correspondences
        R = rotationVectorToMatrix(p(1:3));
        T = p(4:6);
        c_3d = R*w_3d + T;
        J = jacobian(intrinsic_Matrix, p(1:3), w_3d, c_3d);
        e = reprojection_error(intrinsic_Matrix, rotationVectorToMatrix(p(1:3)), p(4:6) , c_2d, w_3d);
        % W = thresholding(e);

        % use_weights is a boolean variable which tell the function whether
        % to use weight matrix from tukey estimator or not
        if use_weights == 0
            update = - (J'*J + lambda*eye(size(J,2)))\(J'*e);
        else
            update = - (J'*W*J + lambda*eye(size(J,2)))\(J'*W*e);
        end

        % calculating reprojection error between projected 3d world point
        % and matched corresponding 2d image points
        p_temp = p + update;
        new_e = reprojection_error(intrinsic_Matrix, rotationVectorToMatrix(p_temp(1:3)), p_temp(4:6) , c_2d, w_3d);

        % update of parameters based on the comparison between old reprojection
        % error and new reprojection error
        if(norm(new_e) > norm(e))
           lambda = 10*lambda ;
        else
           disp(i)
           fprintf('Reprojection error: %f \n',norm(new_e));
           fprintf('Update : %f \n',norm(update));
           lambda = lambda*0.1;
           p = p_temp;
        end

        u = norm(update);
        i = i+1;
    end

end 
