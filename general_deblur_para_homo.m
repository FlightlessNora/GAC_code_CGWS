%% 真实图像deblur时需要提前设置的参数。
%% 注意这些需要在外部设置好的参数：
function [rotpadx,rotpady,maxx,maxy] = general_deblur_para_homo(grid_x,grid_y,grid_z,length_add,cx,cy)
maxz = max(abs(max(grid_z)),abs(min(grid_z)))*pi/180;
% maxz为z轴旋转的最大角度，正负向都包括，以此来判断需要加多宽的边界
maxx = ceil(max(abs(grid_x)))+1;
maxy = ceil(max(abs(grid_y)))+1;
l1 = cx-cy*tan(maxz/2);
l2 = tan(maxz)*l1;
k2 = cy-cx*tan(maxz/2);
k1 = tan(maxz)*k2;
% 加边界的宽度在演算纸上进行了推导，上下加的边宽为k_1*k_2/sqrt(k_1^2+k_2^2)
% k_1和k_2分别为2,4象限被截去的三角形的直角边
% 下面两行区别在于length_add是为了得到有边界条件的模糊图像
rotpadx = ceil(k1*k2/sqrt(k1^2+k2^2))+1+length_add(1);
% rotpadx = ceil(k1*k2/sqrt(k1^2+k2^2));
% 左右加的边宽为l_1*l_2/sqrt(l_1^2+l_2^2)
% l_1和l_2分别为1,3象限被截去的三角形的直角边
rotpady = ceil(l1*l2/sqrt(l1^2+l2^2))+1+length_add(2);
end