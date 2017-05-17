function rot_inv_mat = build_homography(grid_x,grid_y,grid_z,scale,cali_matrix)
% function build_homography_zrot is to build all the homography of the
% camera motion of x,y-translation and z-rotation
% rot_mat: from blurry to sharp
% rot_inv_mat: from sharp to blurry

%% origin version
% [Y,X,Z] = meshgrid(grid_x,grid_y,grid_z);
% num = scale(1)*scale(2)*scale(3);
% M = zeros(num,9);
% rot_mat = zeros(num,9);
% rot_inv_mat = zeros(num,9);
% 
% for i =1:num
%     P = [cos(Z(i)),sin(Z(i)),X(i);-sin(Z(i)),cos(Z(i)),Y(i);0,0,1];
%     H = cali_matrix*P/cali_matrix;
%     H_inv = inv(H);
%     rot_mat(i,:) = reshape(H,[1,9]);
%     rot_inv_mat(i,:) = reshape(H_inv,[1,9]);
% end

%% new version
[TY,TX,RZ] = meshgrid(grid_x,grid_y,grid_z);
R = sqrt(RZ.*RZ);
RZ = RZ./R;
R_sin = sin(R);
R_cos = cos(R);
R_cosc= 1-R_cos;
temp1 = R_cos;
temp2 = R_cos;
temp3 = R_cos+RZ.*RZ.*R_cosc;
temp8 = RZ.*R_sin;
temp9 = -RZ.*R_sin;

num = scale(1)*scale(2)*scale(3);
rot_inv_mat = zeros(num,9);
cali_inv = eye(3)/cali_matrix;

for i =1:num
    if(R(i)>0)
%     P = exp(sqrt(-1)*[0,-RZ(i),RY(i);RZ(i),0,-RX(i);-RY(i),RX(i),0])+ [0,0,TX(i);0,0,TY(i);0,0,1];
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
        P = P + [0,0,TX(i);0,0,TY(i);0,0,0];
    else
        P = [1,0,TX(i);0,1,TY(i);0,0,1];
    end
    H = cali_matrix*P*cali_inv;
%     H = H + [0,0,TX(i);0,0,TY(i);0,0,0];
    H_inv = inv(H);
    rot_mat(i,:) = reshape(H,[1,9]);
    rot_inv_mat(i,:) = reshape(H_inv,[1,9]);
end
end