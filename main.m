% This program assumes, 
% Vl_SIFT libraries
% Results of Exercise 1 are stored under ../data/data/images/init_texture_result
% Results of Exercise 2 are stored under ../data/data/images/detection_result
% Results of Exercise 3 are stored under ../data/data/images/tracking_result

%% Exercise 1.0 - Model preparation

PATH='../data/data/model/teabox.ply';
[vertex,faces]=read_ply(PATH);

Adj=[1,1,1,1,0,0,1,1;1,1,1,1,0,1,1,1;1,1,1,1,0,1,1,0;1,1,1,1,1,1,1,0;1,1,1,1,1,1,0,0;1,1,1,1,1,1,0,1;1,1,1,1,1,0,0,1;1,1,1,1,1,0,1,1];
IntrinsicMatrix = [2960.37845 0 0; 0 2960.37845 0; 1841.68855 1235.23369 1];
normals = [0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; -1 0 0; -1 0 0; 0 0 -1;0 0 -1; 0 -1 0; 0 -1 0 ];
cameraParams = cameraParameters('IntrinsicMatrix',IntrinsicMatrix);

DSC_9743=[1376,1020;2239,1066;2310,1114;1347,1133;0,0;0,0;2278,1589;1376,1614];
DSC_9744=[158,933;2202,1150;1922,1232;1402,987;0,0;2176,1619;1913,1734;1420,1406];
DSC_9745 = [1885,857;1937,1148;1537,1141;1585,854;0,0;1918,1650;1550,1656;0,0];
DSC_9746=[2314,974;1709,1195;1462,1095;2074,909;2287,1393;1714,1687;1477,1564;0,0];
DSC_9747=[2293,1124;1307,1117;1357,996;2248,1007;2258,1604;1341,1597;0,0;0,0];
DSC_9748=[1747,1178;1319,938;1583,883;2055,1102;1759,1670;1350,1344;0,0;2044,1575];
DSC_9749 = [1601,1185;1649,911;1938,905;1983,1186;1615,1679;0,0;0,0;1969,1677];
DSC_9750 = [1461,1145;2053,968;2295,1028;1709,1236;1471,1608;0,0;2268,1447;1707,1727];
Training_images=[DSC_9743;DSC_9744;DSC_9745;DSC_9746;DSC_9747;DSC_9748;DSC_9749;DSC_9750];

%Task 1.1 : Get Rotation and Translation matrices
for i=1:8
    image=Training_images((i-1)*8+1:8*i,:);
    imagePoints=image(Adj(i,:)>0,:);
    worldPoints=vertex(Adj(i,:)>0,:);
    disp('Camera pose estimation for image')
    disp(i)
    [worldOrientation,worldLocation] = estimateWorldCameraPose(imagePoints,worldPoints,cameraParams);
    pcshow(worldPoints,'VerticalAxis','Y','VerticalAxisDir','down','MarkerSize',30);
    hold on
    plotCamera('Size',0.05,'Orientation',worldOrientation,'Location',worldLocation);
    rotation(:,:,i) = worldOrientation;
    translation(i,:) = -rotation(:,:,i)*worldLocation';

end
%
hold off
disp('Saving images with corners plotted on them')
%% Task 1.1 : saving results in init_texture_result folder
% plotting the 3d model corners to every image

path = '../data/data/images/init_texture/DSC_97';
path_result ='../data/data/images/init_texture_result/DSC_97';
for j = 1:8
    filename =strcat(string(42 + j),'.jpg');
    filepath =char(strcat(path, filename));
    im_2d = imread(filepath);
    worldPoints=vertex(Adj(j,:)>0,:);
    c_3d = IntrinsicMatrix'*(rotation(:,:,j)*worldPoints' + translation(j,:)');
    c_2d = c_3d(1:2,:)./repmat(c_3d(3,:),[2,1]);
    imshow(im_2d);
    hold on;
    scatter(c_2d(1,:)', c_2d(2,:)', 15,'r');
    filename_output =strcat(string(42 + j),'.fig');
    output_path = char(strcat(path_result, filename_output));
    savefig(output_path);
    hold off;
end
% Task 1.2: Sift Features detector and feature descriptor extraction
%Mapping Sift feature onto 3d model of box

path = '../data/data/images/init_texture/DSC_97';
sift_descriptors = [];
sift_locations = [];
origin = [0,0,0];
if exist('sift_data.mat','file')==2
    disp('Loading sift data .....')
    load('sift_data.mat');
else
    for j=1:8
        disp('Training SIFT features for image ---> ' )
        disp(j);
        filename =strcat(string(42 + j),'.jpg');
        filepath =char(strcat(path, filename));
        RGB = imread(filepath);
        I = single(rgb2gray(RGB));
        [f_image,d_image] = vl_sift(I);
        for i = 1:size(f_image,2)
            x_3d_Cord = intersection_image_world(faces ,normals, origin, f_image(1:2,i)', IntrinsicMatrix' , rotation(:,:,j),translation(j,:),vertex);   
            if ( ~isempty(x_3d_Cord))
                if (size(x_3d_Cord)>=2)
                    disp(size(x_3d_Cord))
                end
                sift_locations = [sift_locations,  x_3d_Cord];
                sift_descriptors = [sift_descriptors , d_image(:,i)];
            end
        end
    end
    %save('Ex1_sift_data.mat', 'sift_locations', 'sift_descriptors');
end

%% Task 1.2 : saving results in init_texture_folder 
% plotting all camera poses and extracted sift features

pcshow(vertex, 'b', 'VerticalAxis','Z','VerticalAxisDir','up', ...
         'MarkerSize',10);
hold on 
pcshow(sift_locations', 'g', 'VerticalAxis','Z','VerticalAxisDir','up', ...
         'MarkerSize',2);
hold on 
for j =1:8
    t = -rotation(:,:,j)'*translation(j,:)';
    pcshow(t', 'r', 'VerticalAxis','Z','VerticalAxisDir','up', ...
         'MarkerSize',5);
    hold on 
    plotCamera('Size',0.02,'Orientation',rotation(:,:,j),'Location',...
        t', 'Color', [1 0 0]);
    hold on 
    
end
hold off
path_result ='../data/data/images/init_texture_result/'
filename_result = 'sift_locations_camera_pose.fig';
resultpath = char(strcat(path_result, filename_result));
savefig(resultpath);


%% Exercise 2: Pose Estimation with PNP: Detection

path = '../data/data/images/detection/DSC_97';
result_path_task2 = '../data/data/images/detection_result/DSC_97';

r_ransac_history = cell(1,24);
t_ransac_history = cell(1,24);
r_optimum_history = cell(1,24);
t_optimum_history = cell(1,24);

%% Iterating through images and calculating the reprojection error
%% and the best inliers

threshold = 15;
num_iters = 1000;
num_points =7;  %For more accurate Pnp evaluation, we are taking 7 instead of 4 points.
inlierRatio = 0;

for j = 1:24
    % extracting sift features from each image
    filename =strcat(string(50 + j),'.jpg');
    filepath =char(strcat(path, filename));
    image = imread(filepath);
    fprintf('******Detecting Image******** : %d \n', j);
    I = single(rgb2gray(image));
    disp('*******Extracting Features******');
    [f_image,d_image] = vl_sift(I);
    
    % matching the sift features for 2D-3D correspondance
    disp('******Matching Features********');
    [matches, scores] = vl_ubcmatch(d_image, sift_descriptors);
    disp('Best Matched descriptors (Ransac)....');
    
    % searching the inliers and implemented ransac_pnp algo from the
    % multiple view geometry psuedo code by Richard Hatley and Andrew
    % Zisserman
    [bestInNum, bestInRatio ,  best_R , best_T, best_inlierIndex ] = ransac(matches,f_image, sift_locations, threshold, num_points, num_iters ,inlierRatio, cameraParams, IntrinsicMatrix');
    best_location = -best_R'*best_T;
    disp('Plotting figure 1')
    
    % plotting the mached features and inliers
    %plot_result(best_R, best_location, vertex,f_image(1:2,:)', f_image(1:2,matches(1,:))',f_image(1:2,matches(1,best_inlierIndex))',filepath);
    plot_best_features(f_image(1:2,:)', f_image(1:2,matches(1,:))',f_image(1:2,matches(1,best_inlierIndex))',filepath);
    resultfilename =strcat(string(50 + j),'.fig');
    resultfilepath =char(strcat(result_path_task2, resultfilename));
    savefig(resultfilepath);
    
    hold off;
    close;
    
    % performing non linear optimization
    disp('Non Linear optimization ... ')
    [p] = Leven_Marq(IntrinsicMatrix', best_R , best_T, matches(:, best_inlierIndex),f_image , sift_locations,0);
    R_optimum = rotationVectorToMatrix(p(1:3));
    T_optimum = p(4:6);
    
    % saving the plots of 3d corners projected on the image using
    % parameters returned by ransac_pnp and Leven_Marq functions
    disp('Plotting figure 2')
    c_3d_ransac = IntrinsicMatrix'*(best_R*vertex' + best_T);
    c_2d_ransac = c_3d_ransac(1:2,:)./repmat(c_3d_ransac(3,:),[2,1]);
    c_3d_lm = IntrinsicMatrix'*(R_optimum*vertex' + T_optimum);
    c_2d_lm = c_3d_lm(1:2,:)./repmat(c_3d_lm(3,:),[2,1]);
    line_list=[1,5,4,2;
                2,1,3,6;
                3,4,7,2;
                4,1,3,8;
                5,1,8,6;
                6,2,5,7; 
                7,3,6,8;
                8,7,5,4];
    imshow(image);
    hold on;
    scatter(c_2d_ransac(1,:)', c_2d_ransac(2,:)', 15,'r','filled');
    for i =1:8
        line(c_2d_ransac(1,line_list(i,1:2))', c_2d_ransac(2,line_list(i,1:2))','Color','r', 'LineWidth',1.5);
        line(c_2d_ransac(1,line_list(i,[1,3]))', c_2d_ransac(2,line_list(i,[1,3]))','Color','r','LineWidth',1.5);
        line(c_2d_ransac(1,line_list(i,[1,4]))', c_2d_ransac(2,line_list(i,[1,4]))','Color','r','LineWidth',1.5);
    end
    
    
    hold on 
    scatter(c_2d_lm(1,:)', c_2d_lm(2,:)', 15,'g','filled');
    for i =1:8
        line(c_2d_lm(1,line_list(i,1:2))', c_2d_lm(2,line_list(i,1:2))','Color','g', 'LineWidth',1.5);
        line(c_2d_lm(1,line_list(i,[1,3]))', c_2d_lm(2,line_list(i,[1,3]))','Color','g','LineWidth',1.5);
        line(c_2d_lm(1,line_list(i,[1,4]))', c_2d_lm(2,line_list(i,[1,4]))','Color','g','LineWidth',1.5);
    end
    
    %plot(c_2d_lm(1,:), c_2d_lm(2,:), '-x');
    %plot(c_2d_lm(1,:), c_2d_lm(2,:), 'b*-', 'LineWidth', 2, 'MarkerSize', 15);
    hold on;
    resultfilename2 =strcat(string(50 + j),'non_linear_optimization.fig');
    resultfilepath2 =char(strcat(result_path_task2, resultfilename2));
    savefig(resultfilepath2);
    hold off;
    close;
    
    % saving params R & T
    r_ransac_history{j} = best_R;
    t_ransac_history{j} = best_T;
    
    r_optimum_history{j} = R_optimum;
    t_optimum_history{j} = T_optimum;

end

% %% Saving Ransac and LM parameters 
% save('Ex2_ransac_and_lm_params.mat', 'r_ransac_history', ...
%     't_ransac_history', 'r_optimum_history', 't_optimum_history');
% 
% %% Plotting and saving all camera poses
% load('EX2_ransac_and_lm_params.mat');

pcshow(vertex, 'b', 'VerticalAxis','Z','VerticalAxisDir','up', ...
     'MarkerSize',20);
hold on

for i=1:24
    t_ransac = -r_ransac_history{i}'*t_ransac_history{i};
    t_optimum = -r_optimum_history{i}'*t_optimum_history{i};
    
    plotCamera('Size',0.02,'Orientation',r_optimum_history{i},'Location',...
        t_optimum', 'Color', [0 1 0]);
    hold on
    disp(i);
end

% saving the plot
path_result ='../data/data/images/detection_result/';
filename_result = 'all_camera_poses.fig';
resultpath = char(strcat(path_result, filename_result));
savefig(resultpath);



%% Exercise 3 : Pose refinement with nonlinear optimization: Tracking 

%% Initializing path and history variables

path = '../data/data/images/tracking/DSC_9';
result_path_task4 = '../data/data/images/tracking_result/DSC_9';
R_history = cell(1,48);
T_history = cell(1,48);


%% R & T matrices of the first image

% The best_inlierindex from the ransac.m is used here for the tracking.

% The initial frame pose of the camera is estimated below 
threshold = 8;
num_iters = 1000;
num_points =7;
inlierRatio = 0;
for j = 1
    % extracting sift features from the image
    filename =strcat(string(774 + j),'.jpg');
    filepath =char(strcat(path, filename));
    image = imread(filepath);
    I = single(rgb2gray(image));
    [f_image,d_image] = vl_sift(I);

    % searching for inliers
    disp('Matching descriptors ....');
    [matches, scores] = vl_ubcmatch(d_image, sift_descriptors);
    disp('Best Matched descriptors (Ransac)....');
    [bestInNum, bestInRatio ,  best_R , best_T, best_inlierIndex ] = ransac(matches,f_image, sift_locations, threshold, num_points, num_iters ,inlierRatio, cameraParams, IntrinsicMatrix');
       
    % performing non linear optimization
    [p] = Leven_Marq(IntrinsicMatrix', best_R , best_T, matches(:, best_inlierIndex),f_image , sift_locations, 0);
    
    % updating variables
    R_optimum = rotationVectorToMatrix(p(1:3));
    T_optimum = p(4:6);
    R_history{j} = R_optimum;
    T_history{j} = T_optimum;
end

% plotting the projected 3d corners on the image
c_3d_lm = IntrinsicMatrix'*(R_optimum*vertex' + T_optimum);
c_2d_lm = c_3d_lm(1:2,:)./repmat(c_3d_lm(3,:),[2,1]);
imshow(image);
hold on 
scatter(c_2d_lm(1,:)', c_2d_lm(2,:)', 15,'g','filled');
resultfilename4 =strcat(string(774 + j),'trackig.fig');
resultfilepath4 =char(strcat(result_path_task4, resultfilename4));
savefig(resultfilepath4);
close;

% initilizing variables used for future iterations.
R_initial = R_optimum;
T_initial = T_optimum;

%% Iteration on remaining images
 
for j = 2:47
    % extracting the sift features
    filename =strcat(string(774 + j),'.jpg');
    filepath =char(strcat(path, filename));
    image = imread(filepath);
    I = single(rgb2gray(image));
    [f_image,d_image] = vl_sift(I);
    disp(j);

    % searching for inliers
    disp('Matching descriptors ....');
    [matches, scores] = vl_ubcmatch(d_image, sift_descriptors);
    disp('Best Matched descriptors (Ransac)....');
    [bestInNum, bestInRatio ,  best_R , best_T, best_inlierIndex ] = ransac(matches,f_image, sift_locations, threshold, num_points, num_iters ,inlierRatio, cameraParams, IntrinsicMatrix');
  
    % performing non linear optimization
    [p] = Leven_Marq(IntrinsicMatrix', R_initial , T_initial, matches(:,best_inlierIndex), f_image, sift_locations, 0);
    
    % updating variables
    R_optimum = rotationVectorToMatrix(p(1:3));
    T_optimum = p(4:6);
    R_history{j} = R_optimum;
    T_history{j} = T_optimum;
    R_initial = R_optimum;
    T_initial = T_optimum;
    
    % plotting the projected 3d corners on the image
    c_3d_lm = IntrinsicMatrix'*(R_optimum*vertex' + T_optimum);
    c_2d_lm = c_3d_lm(1:2,:)./repmat(c_3d_lm(3,:),[2,1]);
    imshow(image);
    hold on 
    scatter(c_2d_lm(1,:)', c_2d_lm(2,:)', 15,'g','filled');
    resultfilename4 =strcat(string(774 + j),'trackig.fig');
    resultfilepath4 =char(strcat(result_path_task4, resultfilename4));
    savefig(resultfilepath4);
    close;
end

%% Saving R & T parameters of tracking
save('Ex3_rot_&_trans_tracking_params.mat', 'R_history', ...
    'T_history');

% Plotting and saving all camera poses of tracking
load('Ex3_rot_&_trans_tracking_params.mat');


pcshow(vertex, 'b', 'VerticalAxis','Z','VerticalAxisDir','up', ...
     'MarkerSize',20);
hold on

for i=1:47
    t_tracking = -R_history{i}'*T_history{i};

    
    plotCamera('Size',0.02,'Orientation',R_history{i},'Location',...
        t_tracking', 'Color', [1 0 0]);
    hold on

end

% saving the plot
path_result ='../data/data/images/tracking_result/';
filename_result = 'all_tracking_camera_poses.fig';
resultpath = char(strcat(path_result, filename_result));
savefig(resultpath);

% THE END


