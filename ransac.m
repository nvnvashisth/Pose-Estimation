% This function returns the best rotation matrix, best translation vector and 
% inliers calculated by using estimateWorldCameraPose function with
% randomly selecting points for some number of iterations.

% matches       : contains the indexes of matched descriptors of sift features
%                 3d and sift features extracted from an image

% f_image       : locations sift features extracted from an image

% sift_locations: 3d locations of all sift features of box extracted from
%                 the init_texture folder.

% num_points    : number of points given to estimateWorldCameraPose
%                 function in matlab

% iter          : maximum the number of iterations
% threshold     : the threshold of the distances between actual points and 
%                 the projected point

% inlierRatio   : the minimum ratio of number of inliers to the total
%                 mached features


function [bestInlierNum, bestInlierRatio  , best_R , best_T, best_inlierIdx ] = ransac(sift_feature_matches, f_image, sift_locations, threshold, num_points, iterations ,inlierRatio, cameraparams, intrinsic_matrix)
     
     bestInlierNum = 0;  %Number of points out of num_points which gave the least reprojection error 
     bestInlierRatio = 0; %Number of best inliers/total number of points used by PNP
     best_R = [];
     best_T = [];
     
     for i=1:iterations
         % Randomly select 'num_points' number of points and call
         % estimateWorldCameraPose, these points are used by pnp to
         % estimate camera pose 
         idx = randperm(size(sift_feature_matches,2),num_points);
         sample_index = sift_feature_matches(:,idx);
         sample_2d = f_image(1:2,sample_index(1,:));
         sample_3d = sift_locations(1:3,sample_index(2,:));
         data = struct('imagePoints',sample_2d', 'worldPoints',sample_3d','params',cameraparams);
         %% test without reprojectionerror tbd
         
         try
             [R,T] = estimateWorldCameraPose(data.imagePoints,data.worldPoints,data.params, 'MaxReprojectionError', 100);
         catch
             %disp('Error found');
             continue;
         end
         
         % reproject the 3d matched points on the image and calculate the
         % reprojection error
         t = -R*T';
         world_to_camera = R*sift_locations(:, sift_feature_matches(2,:)) + repmat(t, [1 size(sift_feature_matches, 2)]);
         proj_2d = intrinsic_matrix * world_to_camera  ;
         proj_2d_homogeneus = proj_2d(1:2, :) ./ repmat(proj_2d(3, :), [2 1]);
         v = proj_2d_homogeneus - f_image(1:2, sift_feature_matches(1,:));
         vec = [];
         for j=1:size(sift_feature_matches,2)
             vec(:,j) = norm(v(:,j));
         end
         
         % Compute the inliers with distances smaller than the threshold
         inlierIdx = find(abs(vec)<=threshold);
         inlierNum = length(inlierIdx);
         
         % Update the number of inliers and fitting model if better model is found
         if (inlierNum>=round(inlierRatio*size(sift_feature_matches, 2))) && (inlierNum>bestInlierNum)
             bestInlierNum = inlierNum;
             best_R = R;
             best_T = t;
             bestInlierRatio = bestInlierNum/size(sift_feature_matches, 2);
             fprintf('Best Inlier Number: %d \n',bestInlierNum);
             best_inlierIdx = inlierIdx;
         end
    end 
end
