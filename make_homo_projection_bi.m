function proj_set = make_homo_projection_bi(rot_inv_mat,loc,mode)
% calculate projection correspondence 
% proj_set: from sharp to blurry

num_loc = size(loc,1);
num = size(rot_inv_mat,1);
proj_set = cell(1,num_loc);

for i=1:num_loc
    temp = zeros(9,3);
    temp(1,1) = loc(i,1);
    temp(4,1) = loc(i,2);
    temp(7,1) = 1;
    temp(2,2) = loc(i,1);
    temp(5,2) = loc(i,2);
    temp(8,2) = 1;
    temp(3,3) = loc(i,1);
    temp(6,3) = loc(i,2);
    temp(9,3) = 1;
    proj_temp = rot_inv_mat*temp;
    denom = proj_temp(:,3)+1e-10;
    
    location = proj_temp(:,1:2)./repmat(denom,[1,2])-repmat([loc(i,1),loc(i,2)],[num,1]);
    location_floor = floor(location);
    if(strcmp(mode,'round'))
        proj_set{1,i} = round(location);
    elseif(strcmp(mode,'floor'))
        proj_set{1,i} = floor(location);
    elseif(strcmp(mode,'bilinear'))
        proj_set{1,i} = location_floor;
%         dif = location - location_floor;
%         scal = (1+dif(:,1)./(1-dif(:,1)+0.001)).*(1+dif(:,2)./(1-dif(:,2)+0.001));
%         scal = 1./scal;
%         weight = [weight;scal,scal.*dif(:,1)./(1-dif(:,1)),scal.*dif(:,2)./(1-dif(:,2)),scal.*dif(:,1)./(1-dif(:,1)).*dif(:,2)./(1-dif(:,2))];
    end
end
end