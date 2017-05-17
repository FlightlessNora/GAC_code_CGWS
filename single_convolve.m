function RIM = single_convolve(im_size,Pad_im,trans_xy,homoz_rot,homoz_rot_w,bc,scale,fiberl,rotpad)
%%
t = sub2ind([scale(1),scale(2)],fiberl(1,1),fiberl(1,2));
trans = trans_xy{t};
rotpadx = rotpad(1);
rotpady = rotpad(2);
Rt = padarray(Pad_im,[0,rotpady]);
Rt = padarray(Rt,[rotpadx,0],'circular');
padx1 = Rt(1:rotpadx,1:end-1);
padx2 = Rt(end-rotpadx+1:end,2:end);
Rt(1:rotpadx,:) = [zeros(rotpadx,1),padx1];
Rt(end-rotpadx+1:end,:) = [padx2,zeros(rotpadx,1)];
if strcmp(bc,'zero')
    Rt = imfilter(Rt,trans);
else
    Rt = imfilter(Rt,trans,bc);
end
R_IM = Rt(3:end-2,3:end-2);
n = fiberl(1,3);
nxy = homoz_rot{n};
wei = homoz_rot_w{n};
Rr = MyOp_z_rot(R_IM(:)',nxy,wei);
RIM = reshape(Rr,im_size);
end