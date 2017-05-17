function [trans_xy,homoz_rot,homoz_rot_w] = homoz_rot_bilinear(orig_x,orig_y,mtx,...
    mty,transx,transy,scale,inv_rotz_mat)
% function [trans_xy,homoz_rot,trans_xy_t,homoz_rot_w] = homoz_rot_bilinear(orig_x,orig_y,mtx,...
%     mty,transx,transy,scale,inv_rotz_mat)
dim_im_small = orig_x*orig_y;
[y,x] = meshgrid(1:orig_y,1:orig_x);

temp1 = zeros(9,dim_im_small);
temp2 = zeros(9,dim_im_small);
temp3 = zeros(9,dim_im_small);

X_row = x(:)';
Y_row = y(:)';
temp1(1,:) = X_row;
temp2(2,:) = X_row;
temp3(3,:) = X_row;
temp1(4,:) = Y_row;
temp2(5,:) = Y_row;
temp3(6,:) = Y_row;
temp1(7,:) = ones(1,dim_im_small);
temp2(8,:) = ones(1,dim_im_small);
temp3(9,:) = ones(1,dim_im_small);

% 这里为什么用inv_rotz_mat计算proj1,2,3，原因解释见gupta的纸质论文，我用钢笔写在
% 第一页下方
proj1 = inv_rotz_mat*temp1;
proj2 = inv_rotz_mat*temp2;
proj3 = inv_rotz_mat*temp3;

clear temp1;clear temp2;clear temp3;

location1 = proj1./proj3;
location1_floor = floor(location1);
location2 = proj2./proj3;
location2_floor = floor(location2);
dif1 = location1 - location1_floor;
dif2 = location2 - location2_floor;
% weight = [(1-dif1).*(1-dif2),dif1.*(1-dif2),dif2.*(1-dif1),dif1.*dif2];
% 上面的系数的第二种算法：
scal = (1+dif1./(1-dif1)).*(1+dif2./(1-dif2));
scal = 1./scal;
weight = [scal,scal.*dif1./(1-dif1),scal.*dif2./(1-dif2),scal.*dif1./(1-dif1).*dif2./(1-dif2)];
% 核的第一种取法
homoNewZ = (2*repmat(Y_row,[scale(3),1])-location2_floor-1)*orig_x...
    +(2*repmat(X_row,[scale(3),1])-location1_floor);
% 核的第二种取法
% homoNewZ = (location2_floor-1)*orig_x + location1_floor;
temp_loc = [homoNewZ,homoNewZ-1,homoNewZ-orig_x,homoNewZ-orig_x-1];
homoz_rot = cell(1,scale(3));
homoz_rot_w = cell(1,scale(3));
for i = 1:scale(3)
    homoz_rot{i} = temp_loc(i,:);
    homoz_rot_w{i} = weight(i,:);
end
% 注意这里应该与hz的trans xy的生成形式一致，否则classify_fiber and slice就没有
% 太多的意义了
trans_xy = cell(1,scale(1)*scale(2));
% trans_xy_t = cell(1,scale(1)*scale(2));
for i = 1:scale(1)*scale(2)
    dx = transx(i);
    dy = transy(i);
    kerl = [2*mtx+1,2*mty+1];
    ker = zeros(kerl);
    lx = mtx+1+dx;
    ly = mty+1+dy;
    lxf = floor(lx);
    lyf = floor(ly);
    dif1 = lx - lxf;
    dif2 = ly - lyf;
    scal = (1+dif1./(1-dif1)).*(1+dif2./(1-dif2));
    scal = 1./scal;
    weight = [scal,scal.*dif1./(1-dif1),scal.*dif2./(1-dif2),scal.*dif1./(1-dif1).*dif2./(1-dif2)];
    lx_loc = [lxf,lxf-1,lxf,lxf-1];
    ly_loc = [lyf,lyf,lyf-1,lyf-1];
    for j = 1:4
        ker(lx_loc(j),ly_loc(j)) = weight(j);
    end
%     核的第一种取法
    trans_xy{i} = rot90(ker,2);
%     核的第二种取法
%     trans_xy{i} = ker;
end
end