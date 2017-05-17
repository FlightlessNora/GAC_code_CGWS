% ���������������ƽ���˶����õ��Ĳ���
% �������ͼƬbutchershop,bookshelfʱ�Ĳ�������Gupta����ҳ�ϲ鵽�ġ�������
% �����������������Ĳ�����Ҫ������
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
% ����Ϊʲô�������෴��ԭ�������Camera Calibration with One-Dimensional Objects
% ����4.1 Computer Simulations�����ᵽ��u_0 = 320,v_0 = 240��ʱ��ͼ��ֱ���Ϊ640*480
% ��ͼ���СΪ480*640
% C_intrinsic = [alpha_x gamma cx;0 alpha_y cy;0 0 1];
f1_d = alpha_x/depth;
f2_d = alpha_y/depth;
K = [f1_d,0,u_x;0,f2_d,u_y;0,0,1];
% clearvars -except f1_d f2_d focal_length