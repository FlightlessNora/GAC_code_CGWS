function R = slice_trans_tr(im_size,rotpadx,rotpady,Pad_im,trans_slice,homoz_rot,homoz_rot_w,bc,var)
% 把第n个slice内的平移变换计算加权和
if var == 1
    Tt = padarray(Pad_im,[0,rotpady]);
    Tt = padarray(Tt,[rotpadx,0],'circular');
    padx1 = Tt(1:rotpadx,1:end-1);
    padx2 = Tt(end-rotpadx+1:end,2:end);
    Tt(1:rotpadx,:) = [zeros(rotpadx,1),padx1];
    Tt(end-rotpadx+1:end,:) = [padx2,zeros(rotpadx,1)];
    if strcmp(bc,'zero')
        Tim = imfilter(Tt,trans_slice);
    else
        Tim = imfilter(Tt,trans_slice,bc);
    end
    Tim = Tim(3:end-2,3:end-2);
    R = MyOp_z_rot(Tim(:)',homoz_rot,homoz_rot_w);
    R = reshape(R,im_size);
else
    %% 先将该slice要处理的图像进行旋转
    Rr = MyOp_conj_trans_z_rot(Pad_im(:)',homoz_rot,homoz_rot_w);
    Rr = reshape(Rr,im_size);
    Rt = padarray(Rr,[0,rotpady]);
    Rt = padarray(Rt,[rotpadx,0],'circular');
    padx1 = Rt(1:rotpadx,1:end-1);
    padx2 = Rt(end-rotpadx+1:end,2:end);
    Rt(1:rotpadx,:) = [zeros(rotpadx,1),padx1];
    Rt(end-rotpadx+1:end,:) = [padx2,zeros(rotpadx,1)];
    %% 整个slice进行平移
    if strcmp(bc,'zero')
        R = imfilter(Rt,trans_slice);
    else
        R = imfilter(Rt,trans_slice,bc);
    end
    R = R(3:end-2,3:end-2);
end
end