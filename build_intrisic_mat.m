% 这里设置所有相机平移运动所用到的参数
% 相机拍摄图片butchershop,bookshelf时的参数，从Gupta的网页上查到的。如果相机
% 设置又所变更，这里的参数需要调整。
% focal_length = 50mm % 36mm for houses.jpg
% sensor resolution: 40pixels/mm
% depth: 2.5m
% sensor size: 640X480 pixels
% ccd_size = [16,12]; % ccd_size = [36,36]; for hz's picture.
% gamma:represents the skew coefficient between the x and the y axis,and
% often is 0
function [K,cx,cy]=build_intrisic_mat(im_size,focal_length)
% function [f1_d,f2_d,K,cx,cy]=build_intrisic_mat(im_size,focal_length)

% gamma = 0;
max_width_of_image = max(im_size);
max_height_of_image = max_width_of_image;
depth = 1;
ccd_size = [36,36];
max_width_of_ccd = ccd_size(2);
max_height_of_ccd = ccd_size(1);
m_x = max_width_of_image/max_width_of_ccd;
m_y = max_height_of_image/max_height_of_ccd;
alpha_x = focal_length*m_x;
alpha_y = focal_length*m_y;
orig_x = im_size(1);
orig_y = im_size(2);
cx = orig_x/2;
cy = orig_y/2;
u_x = orig_y/2;
u_y = orig_x/2;
% 这里为什么横纵轴相反，原理见论文Camera Calibration with One-Dimensional Objects
% 文中4.1 Computer Simulations部分提到，u_0 = 320,v_0 = 240此时，图像分辨率为640*480
% 即图像大小为480*640
% C_intrinsic = [alpha_x gamma cx;0 alpha_y cy;0 0 1];
f1_d = alpha_x/depth;
f2_d = alpha_y/depth;
K = [f1_d,0,u_x;0,f2_d,u_y;0,0,1];
% clearvars -except f1_d f2_d focal_length