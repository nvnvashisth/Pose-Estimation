% The function finds out an intersection between Camera coordinates and
% world coordinates. The function first finds out if there is an
% intersection(i.e. does the ray points to the box) and then finds out the
% the face of the which intersects with the ray.

%Input: model_faces of the model w.r.t world 
%       normals w.r.t world 
%       origin w.r.t camera
%       sift_location: disk center of the image features 
%       A, R, T are the Intrinsic, Rotation and Translation matrices
%       model_vertex w.r.t world 

%Output: Interesection points in the world corrdinate system


function [d_world_coord] = intersection_image_world(model_faces, normals, origin, sift_location, A, R , T, model_vertex)

    model_faces = model_faces + repmat(ones(size(model_faces,1),1), [1, size(model_faces,2)]);
    dir = A\[sift_location,1]';
    
    best_coordinates = [];
    
    for i=1:size(model_faces,1)
        
        face = model_faces(i,:);
        
        % extracting the world coordinates of vertex for each face 
        v1 = model_vertex(face(1),:); 
        v2 =model_vertex(face(2),:);
        v3 = model_vertex(face(3),:);
        
        % transforming the vertex from world coordinates to camera coordinates 
        v1_camera = R*v1' + T';         
        v2_camera = R*v2' + T';         
        v3_camera = R*v3' + T';
        
        % The TriangleRayIntersection is called by keeping everything in camera coordinates 
        [intersect, t, u, v, x_coordinate] = TriangleRayIntersection (...
                                 origin, dir , v1_camera, v2_camera, v3_camera); 
                         
       if (intersect == true)
           normal = normals(i,:);
           x_coordinate = R'*(x_coordinate' - T');
           camera_center = -R'*T';
           direction_camera_to_xcoordinate= (x_coordinate-camera_center );
           m_normal = direction_camera_to_xcoordinate/norm(direction_camera_to_xcoordinate);
           
           temp = dot(m_normal', normal);
           if temp < 0.0
               best_coordinates = x_coordinate;
           end
       end
    end
    
    d_world_coord = best_coordinates;
end 