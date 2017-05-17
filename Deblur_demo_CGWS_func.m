function Deblur_demo_CGWS_func(maxit,innerit,im_id,win_num)
pat = pwd;
path(path,strcat(pat,'\qpc'))
path(path,strcat(pat,'\Multivariate Gaussian Distribution'))
path(path,strcat(pat,'\Images'))
RP = strcat(pat,'\Results_inner_',num2str(innerit),'_out_',num2str(maxit),'_CGWS');
if ~exist(RP)
    mkdir(RP);
end
results_path = RP;
ksetmat = strcat('blurred_kset_',num2str(im_id),'_',num2str(win_num),'.mat');
load (ksetmat);
col_file_name = strcat('blurred',num2str(im_id),'.jpg');
name = strfind(col_file_name,'.');
imname = col_file_name(1:(name-1));
col_im = im2double(imread(col_file_name));
if(size(col_im,3)>1)
    im = rgb2gray(col_im);
else
    im = col_im;
end
im = im2double(im);
im_size = size(im);
dim_im = im_size(1)*im_size(2);
length_add = [0,0];
k_size = 21;
pose_ind_size = 60;
bc = 'circular';
focal_length = 35;
unit = [0.0013,0.0013,0.0012];
rx_lim = [-6*unit(1),6*unit(1)];
ry_lim = [-6*unit(2),6*unit(2)];
rz_lim = [-7*unit(3),7*unit(3)];

[grid_x,grid_y,grid_z,scale] = make_kernel_grid(rx_lim,ry_lim,rz_lim,unit);

[K,cx,cy]=build_intrisic_mat(im_size,focal_length);

[gridx,gridy] = build_homography_trans(grid_x,grid_y,scale,K);

inv_rot_mat = build_homography_zrot(grid_z,scale(3),K);
[rotpadx,rotpady,mtx,mty] = general_deblur_para_homo(gridx,gridy,grid_z,length_add,cx,cy);

[trans_xy,homoz_rot,homoz_rot_w] = homoz_rot_bilinear(im_size(1),im_size(2),mtx,mty,gridx,gridy,...
    scale,inv_rot_mat);
traj_inv_mat = build_homography(grid_x,grid_y,grid_z,scale,K);
proj_set = make_homo_projection_bi(traj_inv_mat,loc,'bilinear');
[pose_ind,w0] = pose_restrict_bi(k_size,k_size,k_set,proj_set,loc,pose_ind_size);

[pose_sub_indx,pose_sub_indy,pose_sub_indz] = ind2sub(scale,pose_ind);
pose_sub_ind = [pose_sub_indx,pose_sub_indy,pose_sub_indz];
[px_mat,py_mat,pxx_mat,pxy_mat,pyy_mat] = gen_partialmat(im_size(1),im_size(2));

%% main loop
sigma_s = 2;sigma_r = 0.5;
support_size = 5;
dt = 1;
thresh = truncate_thresh(im,k_size,2,.5);
sigma_r = sigma_r*0.93^0;
dt = dt*0.93^0;
thresh = thresh*0.93^0;

mode = 'Gaussian';
iter = 1;
err =1;
tol = 1e-6;
L = im;

alpha = 8;
fs_alpha = 0.2;
%%
tic
while (iter<=maxit)&&(err>tol)
    
    [Px,Py] = prediction(L,sigma_s,sigma_r,support_size,dt,thresh);
    sigma_r = sigma_r*0.93;
    dt = dt*0.93;
    thresh = thresh*0.93;
    
    w = pose_weight_estimatation(Px,Py,im,w0,pose_ind_size,trans_xy,homoz_rot,...
        homoz_rot_w,rotpadx,rotpady,bc,dim_im,im_size,scale,pose_sub_ind);
    
    [slicel,trans_slice,trans_slice_t,Pose_sub_ind] = classify_all_slice(w,fs_alpha,trans_xy,scale,pose_sub_ind);
    % [fiberl,slicel,trans_slice,trans_slice_t,Pose_sub_ind,w] = classify_fiber_slice(w,fs_alpha,trans_xy,scale,pose_sub_ind);
    L_old = L;
    %     parameters for CGWS
    L = spv_deconv_partial(L_old,im,alpha,px_mat,py_mat,pxx_mat,pxy_mat,...
        pyy_mat,slicel,im_size,rotpadx,rotpady,trans_slice,...
        homoz_rot,homoz_rot_w,trans_slice_t,bc,Pose_sub_ind);
    
    ind_used = pose_ind;
    [pose_ind,w0] = pose_perturbation(ind_used,w,grid_x,grid_y,grid_z,mode);
    [pose_sub_indx,pose_sub_indy,pose_sub_indz] = ind2sub(scale,pose_ind);
    pose_sub_ind = [pose_sub_indx,pose_sub_indy,pose_sub_indz];
    
    iter = iter +1;
    err = sqrt(sum(sum((L-L_old).^2)))/size(L,1)/size(L,2);
end
alpha_end = 1;
L = spv_deconv_partial(L,im,alpha_end,px_mat,py_mat,pxx_mat,pxy_mat,...
        pyy_mat,slicel,im_size,rotpadx,rotpady,trans_slice,homoz_rot,...
        homoz_rot_w,trans_slice_t,bc,Pose_sub_ind);
Elaspedtime1 = toc;
%final deconvolution
alpha_end = 1;
if size(col_im,3)>1
    L_end = zeros(size(col_im));
    for i = 1:size(col_im,3)
        L_end(:,:,i) = spv_deconv_partial(L,col_im(:,:,i),alpha_end,px_mat,py_mat,pxx_mat,pxy_mat,...
            pyy_mat,slicel,im_size,rotpadx,rotpady,trans_slice,homoz_rot,...
            homoz_rot_w,trans_slice_t,bc,Pose_sub_ind);
    end
    deblurred_matnameL_end = strcat(imname,'_L_CGWS');
    deblurred_jpgnameL_end = strcat(imname,'_L_CGWS.jpg');
    Elaspedtime2 = toc;
    figure,imshow(L_end,[])
    title(num2str(Elaspedtime2));
    saveas(gcf,fullfile(results_path,deblurred_matnameL_end))
    close(gcf)
    save(fullfile(results_path,deblurred_matnameL_end),'L_end')
    imwrite(L_end,fullfile(results_path,deblurred_jpgnameL_end))
end
deblurred_matnameL = strcat(imname,'_L_gray_CGWS');
deblurred_jpgnameL = strcat(imname,'_L_gray_CGWS.jpg');
if size(col_im,3)==1
    imwrite(L,fullfile(results_path,deblurred_jpgnameL))
end
figure,imshow(L,[])
title(num2str(Elaspedtime1));
saveas(gcf,fullfile(results_path,deblurred_matnameL))
close(gcf)
save(fullfile(results_path,deblurred_matnameL),'L')
deblurred_kermat = strcat(imname,'_kermat');
save(fullfile(results_path,deblurred_kermat),'L','col_im','results_path',...
    'traj_inv_mat','w','ind_used','L_end')
clc