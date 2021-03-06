function [slice_set,trans_slice,trans_slice_t,pose_sub_ind] = classify_all_slice(w0,alpha,trans_xy,scale,pose_sub_ind)

w_0 = find(w0==0);
w0(w_0) = [];
pose_ind_size = size(w0,1);
pose_sub_ind(w_0,:) = [];

pose_sub_xy = pose_sub_ind(:,1:2);
pose_sub_z = pose_sub_ind(:,3);
[slicel,num_fiber] = find_same_ele(pose_sub_z,pose_ind_size);
[fiberl,num_slice] = find_same_ele(pose_sub_xy,pose_ind_size);
slice_weight = compute_weight(slicel,num_fiber,w0);
fiber_weight = compute_weight(fiberl,num_slice,w0);
error0 = sum(slice_weight);
error = error0;
n_slice = 0;
tol = 1e-10;
iter = 0;
alpha_max = alpha;
while(error/error0>tol)
    w_f = max(fiber_weight);
    w_s = max(slice_weight);
    if (w_s>(alpha*w_f))
        s_loc = find(slice_weight == w_s);
        if (length(s_loc)>1)
            s_loc = s_loc(1);
        end
        i_s = slicel{s_loc};
        n_slice = n_slice+1;
        slice_set{n_slice} = i_s;
        ml = sub2ind([scale(1),scale(2)],pose_sub_ind(i_s,1),pose_sub_ind(i_s,2));
        ker = 0;
        for m = 1:length(i_s)
            ker = ker+w0(i_s(m))*trans_xy{ml(m)};
        end
        trans_slice{n_slice} = ker;
        trans_slice_t{n_slice} = rot90(ker,2);
        
        slice_weight(s_loc) = 0;
        [loc_reduce_f,cleared_set,num,label] = clear_reduced_location(i_s,fiberl);
        if (num==1)
            fiber_weight(label) = fiber_weight(label)-sum(w0(loc_reduce_f{label}));
            fiberl{label}(cleared_set{label}) = [];
        elseif (num>1)
            for j = 1:num
                la = label(j);
                fiber_weight(la) = fiber_weight(la)-sum(w0(loc_reduce_f{la}));
                fiberl{la}(cleared_set{la}) = [];
            end
        end
        error = sum(slice_weight)+sum(fiber_weight);
    else
        rat = ((2*iter-2)/(2*iter-1))^0.125;
        alpha = max(alpha*(0.1*rat + 0.9), alpha_max);
    end
    iter = iter+1;
end
if (n_slice == 0)
    slice_set = [];
end
end
%%
function weight = compute_weight(label,num,w0)
weight = zeros(1,num);
for k = 1:num
    L = label{k};
    weight(k) = sum(w0(L));
end
end

function [value,cleared_set,num,label] = clear_reduced_location(chosen,...
        set_to_be_clear)
    L = length(set_to_be_clear);
    value = cell(1,L);
    cleared_set = cell(1,L);
    num = 0;
    label = [];
    for i = 1:L
        set = set_to_be_clear{i};
        [v,~,c] = intersect(chosen,set);
        if ~isempty(c)
            num = num+1;
            label(num) = i;
            value{i} = v;
            cleared_set{i} = c;
        end
    end
end
%%
function [same,num_same] = find_same_ele(traj,l)
if (size(traj,2)==1)
    [xx,xy] = sort(traj);
    xx1 = [xx;0];
    xx2 = [0;xx];
    xxnew = (xx1-xx2);
    xxnew = xxnew(2:l+1);
    [xy1,~] = find(xxnew~=0);
    num_same = length(xy1);
    same = cell(1,num_same);
    same_ind = 1;
    for i = 1:num_same
        ind = xy1(i);
        same{i} = xy(same_ind:ind);
        same_ind = ind+1;
    end
elseif (size(traj,2)==2)
    trajx = traj(:,1);
    trajy = traj(:,2);
    [xx,xy] = sort(trajx);
    xx1 = [xx;0];
    xx2 = [0;xx];
    xxnew = (xx1-xx2);
    xxnew = xxnew(2:l+1);
    [xy1,~] = find(xxnew~=0);
    length_samex = length(xy1);
    samex_ind = 1;
    num_same = 1;
    for i = 1:length_samex
        ind_x = xy1(i);
        samex = xy(samex_ind:ind_x);
        l = length(samex);
        [ty,tyy] = sort(trajy(samex));
        ty1 = [ty;0];
        ty2 = [0;ty];
        yynew = ty1-ty2;
        yynew = yynew(2:l+1);
        [yy1,~] = find(yynew~=0);
        length_samey = length(yy1);
        samey_ind = 1;
        if (length_samey==1)
            ind_y = yy1;
            same{num_same} = samex(tyy(samey_ind:ind_y));
            samey_ind = ind_y+1;
            num_same = num_same+1;
        else
            for j = 1:length_samey
                ind_y = yy1(j);
                same{num_same} = samex(tyy(samey_ind:ind_y));
                samey_ind = ind_y+1;
                num_same = num_same+1;
            end
        end
        samex_ind = ind_x+1;
    end
    num_same = num_same-1;
elseif (size(traj,2)>2)
    disp('the dimension of the vecotr must be a column vector or m*2 vector');
end
end