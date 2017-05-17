function [grid_x,grid_y,grid_z,scale] = make_kernel_grid(rx_lim,ry_lim,rz_lim,unit)
% This function build kernel grid for 3D camera motion
% grid_size: 1x3 matrix indicates how many grid for each dim
% scale: 1x3 matrix indicates the unit width for each dim


x_scale = rx_lim./unit(1);
y_scale = ry_lim./unit(2);
z_scale = rz_lim./unit(3);
grid_x = ceil(x_scale(1))*unit(1):unit(1):floor(x_scale(2))*unit(1);
grid_y = ceil(y_scale(1))*unit(2):unit(2):floor(y_scale(2))*unit(2);
grid_z = ceil(z_scale(1))*unit(3):unit(3):floor(z_scale(2))*unit(3);

scale = zeros(1,3);
scale(1) = abs(ceil(x_scale(1)))+floor(x_scale(2))+1;
scale(2) = abs(ceil(y_scale(1)))+floor(y_scale(2))+1;
scale(3) = abs(ceil(z_scale(1)))+floor(z_scale(2))+1;
end