function [pose_ind,w0] = pose_restrict_bi(row_k,col_k,...
    kernel_set,proj_set,loc_set,bar)

% This function initialize the camera motion using local blur kernels and
% compute the pose basis and their weight

half_row = (row_k-1)/2;
half_col = (col_k-1)/2;
num_loc = size(loc_set,1);
prob_set = [];
thresh = 0.01; % to determine each point
%% determine which poses to use
for i =1:num_loc
    temp_proj = proj_set{1,i};
    temp_proj = temp_proj + repmat([half_row+1,half_col+1],[size(temp_proj,1),1]);
%     temp_k = rot90(kernel_set{1,i},2);
    temp_k = kernel_set{1,i};

    ind = (temp_proj(:,2)-1)*row_k+temp_proj(:,1);
    ind = [ind,ind+1,ind+col_k,ind+col_k+1];
    temp_value = zeros(size(ind));
    temp_value(temp_k(ind)>thresh)=1;
    prob_set = [prob_set,sum(temp_value,2)];
end

prob = sum(prob_set,2)/4/num_loc;% 这就是论文中的初始W
[~,temp2] = sort(prob,'descend');
pose_ind = temp2(1:bar);
w0 = prob(pose_ind);
w0 = w0./sum(w0);
end