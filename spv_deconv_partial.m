function [L] = spv_deconv_partial(L_old,B,alpha_A,px_mat,py_mat,pxx_mat,...
    pxy_mat,pyy_mat,slicel,im_size,rotpadx,rotpady,trans_slice,homoz_rot,homoz_rot_w,...
    trans_slice_t,bc,pose_sub_ind)
% f = (KL-B)'(KL-B) + alpha_A*(px_mat*L-Px)'(px_mat*L-Px) +alpha_A*(py_mat*L-Py)'(py_mat*L-Py)
% \partial f = 2*(K'K+alpha_A*px_mat'*px_mat+alpha_A*py_mat'*py_mat)*L - 2*(K'B+alpha_A*px_mat'*Px+alpha_A*py_mat'*Py)
% KTK:  K'K
% KTB:  K'B
% ---------- Argument defaults ----------
if ~exist('tol','var') || isempty(tol) tol = 1e-1; end;
if ~exist('maxit','var') || isempty(maxit) maxit = 40; end;
omega = [50,25,12.5];
px_mat_t = px_mat';
py_mat_t = py_mat';
pxx_mat_t = pxx_mat';
pxy_mat_t= pxy_mat';
pyy_mat_t = pyy_mat';

k = L_old(:);
b = KT_gen_convol_spv(omega,B(:),px_mat,px_mat_t,py_mat,py_mat_t,pxx_mat,...
    pxx_mat_t,pxy_mat,pxy_mat_t,pyy_mat,pyy_mat_t,slicel,im_size,...
    rotpadx,rotpady,trans_slice_t,homoz_rot,homoz_rot_w,bc,pose_sub_ind);% 1.580398 seconds.

% ---------- Initialize -------------------------% Current iterate

r = KTK_gen_convol_spv(omega,k,alpha_A,px_mat,px_mat_t,py_mat,py_mat_t,pxx_mat,pxx_mat_t,...
    pxy_mat,pxy_mat_t,pyy_mat,pyy_mat_t,slicel,im_size,rotpadx,rotpady,homoz_rot,homoz_rot_w,...
    trans_slice,trans_slice_t,bc,pose_sub_ind) - b;%2.614439 seconds.

p = -r;                        % Current conjugate
iter = 0;
% ---------- Begin iteration ------------------------
while ((norm(r) > tol) && (iter < maxit))
    rr = r' * r;
    Ap = KTK_gen_convol_spv(omega,p,alpha_A,px_mat,px_mat_t,py_mat,py_mat_t,pxx_mat,pxx_mat_t,...
        pxy_mat,pxy_mat_t,pyy_mat,pyy_mat_t,slicel,im_size,rotpadx,rotpady,homoz_rot,homoz_rot_w,...
        trans_slice,trans_slice_t,bc,pose_sub_ind);% 2.665888 seconds.
    pAp = p' * Ap;
    alpha = rr / (pAp+0.001);    % Compute the step length
    k = k + alpha * p;       % Update x
    r = r + alpha * Ap;   % Update residual
    rr_new = r' * r;
    rho = rr_new / (rr+0.001);  % Compute beta
    p = rho * p - r;        % Update cg vector 
    iter = iter +1;
end
k(k>1) =1;
k(k<0) =0;
L = col2im(k,[1,1],size(B));
% figure,imshow(L,[])
end