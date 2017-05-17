function [pose_ind,w] = pose_perturbation(pose_ind_old,w_old,grid_x,grid_y,grid_z,mode)
% This function add perturbation to possible camera poses.

len_x = length(grid_x);
len_y = length(grid_y);
len_z = length(grid_z);

%% determine which poses to use
if(strcmp(mode,'Gaussian'))
    mu = [0,0,0];
    sigma = [2,0.5,0.5;0.5,2,0.5;0.5,0.5,3];
    sample = round(mgd(1000,3,mu,sigma));
    sample_count = 1;
    
    [temp_value,temp_sortind] = sort(w_old,'descend');
    num_cut = sum((w_old<temp_value(1)*0.05));% 0.05这个数越大，更新的越快
    fprintf('%d ',num_cut);
    
    ind_temp = pose_ind_old(temp_sortind(1:end-num_cut));
    pose_ind = ind_temp;
    [pose_sub1,pose_sub2,pose_sub3] = ind2sub([len_x,len_y,len_z],ind_temp);
    pose_sub = [pose_sub1,pose_sub2,pose_sub3];
    used_sub = pose_sub;
    nn = length(ind_temp);
    count = 0;
    while(count<num_cut)
        add_ind = randi([1,nn]);
        sample_sub = pose_sub(add_ind,:)+sample(sample_count,:);
        if((sample_sub(1)>0) && (sample_sub(1)<=len_x) && (sample_sub(2)>0) ...
                && (sample_sub(2)<=len_y) && (sample_sub(3)>0) && (sample_sub(3)<= len_z))
            
            temp = abs(used_sub-repmat(sample_sub,[size(used_sub,1),1]));
            temp = sum(temp,2);
            if(sum(temp==0)==0)
                ind_new = sub2ind([len_x,len_y,len_z],sample_sub(1),...
                    sample_sub(2),sample_sub(3));
                pose_ind = [pose_ind;ind_new];
                used_sub = [used_sub;sample_sub];
                count = count +1;
            end
        end
        sample_count = sample_count+1;
        if (sample_count>1000)
            sample = round(mgd(1000,3,mu,sigma));
            sample_count = 1;
        end
    end
end

w = zeros(size(w_old));
w(1:end-num_cut) = w_old(temp_sortind(1:end-num_cut));