function [gridx,gridy,rot_inv_mat] = build_homography_trans(grid_x,grid_y,scale,cali_matrix)
% function build_homography_zrot is to build all the homography of the
% camera motion of x,y-translation and z-rotation
% rot_mat: from blurry to sharp
% rot_inv_mat: from sharp to blurry
%% new version
[TY,TX] = meshgrid(grid_x,grid_y);

num = scale(1)*scale(2);
rot_inv_mat = zeros(num,9);
cali_inv = eye(3)/cali_matrix;

for i =1:num
    P = [1,0,TX(i);0,1,TY(i);0,0,1];
    H = cali_matrix*P*cali_inv;
    H_inv = inv(H);
    rot_inv_mat(i,:) = reshape(H_inv,[1,9]);
end
gridx = rot_inv_mat(:,7);
gridy = rot_inv_mat(:,8);
end