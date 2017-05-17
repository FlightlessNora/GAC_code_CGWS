function Rim = KT_gen_conv_spv(slicel,im_size,rotpadx,rotpady,TPad_im,...
    trans_slice_t,homoz_rot,homoz_rot_w,bc,pose_sub_ind)
ls = length(slicel);
Rim = zeros(im_size);
for p = 1:ls
    slice = pose_sub_ind(slicel{p}(1),3);
    Rim = Rim+slice_trans_tr(im_size,rotpadx,rotpady,TPad_im,...
        trans_slice_t{p},homoz_rot{slice},homoz_rot_w{slice},bc,0);
end
Rim = Rim(:);
end