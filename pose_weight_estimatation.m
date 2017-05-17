function w = pose_weight_estimatation(Px,Py,B,w0,n_basis,trans_xy,homoz_rot,...
    homoz_rot_w,rotpadx,rotpady,bc,dim_im,im_size,scale,pose_sub_ind)
% The function is used to estimate camera motion given predicted gradient maps.
% Required input
% Px,Py: predicted gradient maps
% B: blurry image
% omega: weights for partial derivative
% beta: Tihkonov parameter
omega = [25,12.5]; % weight for derivative
[Bx,By] = gradient(B,1);
[Bxx,Bxy1] = gradient(Bx,1);
[Bxy2,Byy] = gradient(By,1);
Bxy = Bxy1 + Bxy2;
%% build data for optimization
A = [];
% K * \partial L
for i=1:n_basis
%     tic
    Lx_basis_temp = single_convolve(im_size,Px,trans_xy,homoz_rot,homoz_rot_w,bc,scale,pose_sub_ind(i,:),[rotpadx,rotpady]);
    Ly_basis_temp = single_convolve(im_size,Py,trans_xy,homoz_rot,homoz_rot_w,bc,scale,pose_sub_ind(i,:),[rotpadx,rotpady]);
    [Lxx,Lxy1] = gradient(Lx_basis_temp,1);
    [Lxy2,Lyy] = gradient(Ly_basis_temp,1);
    Lxy = (Lxy1+Lxy2).*0.5;
    Lx_basis_temp = Lx_basis_temp(:);
    Ly_basis_temp = Ly_basis_temp(:);
    Lxx_basis_temp = Lxx(:);
    Lxy_basis_temp = Lxy(:);
    Lyy_basis_temp = Lyy(:);
    temp = [Lx_basis_temp;Ly_basis_temp;Lxx_basis_temp;Lxy_basis_temp;Lyy_basis_temp];
    A = [A,temp];
end
omega = omega*2; 
temp = A(1:dim_im,:);
ATA = temp'*temp*omega(1);
ATb = temp'*Bx(:)*omega(1);
temp = A(dim_im+1:dim_im*2,:);
ATA = ATA+temp'*temp*omega(1);
ATb = ATb+temp'*By(:)*omega(1);
temp = A(2*dim_im+1:dim_im*3,:);
ATA = ATA+temp'*temp*omega(2);
ATb = ATb+temp'*Bxx(:)*omega(2);
temp = A(3*dim_im+1:dim_im*4,:);
ATA = ATA+temp'*temp*omega(2);
ATb = ATb+temp'*Bxy(:)*omega(2);
temp = A(4*dim_im+1:dim_im*5,:);
ATA = ATA+temp'*temp*omega(2);
ATb = ATb+temp'*Byy(:)*omega(2);

clear bx_temp by_temp bxx_temp bxy_temp byy_temp;
clear Lx_basis_temp Ly_basis_temp Lxx_basis_temp Lxy_basis_temp Lyy_basis_temp;
clear temp;

beta =1;
H = ATA + beta*eye(size(ATA,1));
f = -ATb;
L = ones(1,length(w0));
k = 1;
A = [];
b = [];
l = zeros(length(w0),1);
u = [];

[w,~,~] = qpip(H,f,L,k,A,b,l,u,0,0,0);

w(w<max(w(:))*0.01)=0;
w = w./sum(w(:));
end