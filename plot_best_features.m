
function [] = plot_best_features(f_all, f_matched, f_ransac, image) 
    
    image_2d = imread(image);
    imshow(image_2d);
    hold on
    scatter(f_all(:,1), f_all(:,2), 5,'r'); %Plotting all descriptors of the images
    hold on
    
    scatter(f_matched(:,1), f_matched(:,2), 5, 'b'); %Plotting all matched descriptors
    hold on
    
    scatter(f_ransac(:,1), f_ransac(:,2), 5,'g'); %Plotting best features of the image
    
    hold off
end