function [Px,Py] = prediction(I,sigma_s,sigma_r,support_size,dt,threshold)
% The function is used to generate image gradient maps of the latent image
% in the prediction step.
% 
% Required input
% I: input image
% sigma_s: spatial sigma
% sigma_r: range sigma
% support_size: support size of bilateral filtering
% dt: time step in a single evolution 
% 

%% bilateral filetr
sigma = [sigma_s,sigma_r];
w = floor(support_size/2);
L1 = bfilter2(I,w,sigma);
% figure;imshow(L1);

%% shock filter
iter = 30;
L2 = shock(L1,iter,dt,1,'org');
% figure;imshow(L2);

%% thresholding
[Px,Py] = gradient(L2,1);
magnitude = sqrt(Px.^2+Py.^2);
ind = (magnitude>threshold);
Px = Px.*ind;
Py = Py.*ind;
end