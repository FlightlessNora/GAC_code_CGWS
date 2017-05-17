function rot_inv_mat = build_homography_zrot(grid_z,num,cali_matrix)
% function build_homography_zrot is to build all the homography of the
% camera motion of x,y-translation and z-rotation
% rot_mat: from blurry to sharp
% rot_inv_mat: from sharp to blurry
%% new version
R = sqrt(grid_z.*grid_z);
RZ = grid_z./R;
R_sin = sin(R);
R_cos = cos(R);
R_cosc= 1-R_cos;
temp1 = R_cos;
temp2 = R_cos;
temp3 = R_cos+RZ.*RZ.*R_cosc;
temp8 = RZ.*R_sin;
temp9 = -RZ.*R_sin;
rot_inv_mat = zeros(num,9);
cali_inv = eye(3)/cali_matrix;

for i =1:num
    if(R(i)>0)
        P = zeros(3,3);
        P(1) = temp1(i);
        P(2) = temp8(i);
        P(3) = 0;
        P(4) = temp9(i);
        P(5) = temp2(i);
        P(6) = 0;
        P(7) = 0;
        P(8) = 0;
        P(9) = temp3(i);
    else
        P = [1,0,0;0,1,0;0,0,1];
    end
    H = cali_matrix*P*cali_inv;
    H_inv = inv(H);
    rot_inv_mat(i,:) = reshape(H_inv,[1,9]);
end
end