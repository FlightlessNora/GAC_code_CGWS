function [threshold,threshold2] = truncate_thresh(L,k_size,r,rr)
imsize = size(L);
[L_x,L_y] = gradient(L,1);
magnitude = sqrt(L_x.^2+L_y.^2);
direction = atan(L_y./(L_x+0.0001));

r0 = zeros(imsize);
r1 = r0;r2 = r0;r3 = r0;
r1(direction<-0.25*pi)=1;
r2(direction<0)=1;
r2 = r2-r1;
r3(direction<0.25*pi)=1;
r4 = 1- r3;
r3 = r3-r2;

r1 = r1.*magnitude;
r2 = r2.*magnitude;
r3 = r3.*magnitude;
r4 = r4.*magnitude;

num = r*k_size;
num2 = ceil(rr*k_size*sqrt(prod(imsize)));

value1 = sort(r1(:),'descend');
value2 = sort(r2(:),'descend');
value3 = sort(r3(:),'descend');
value4 = sort(r4(:),'descend');

threshold = min([value1(num),value2(num),value3(num),value4(num)]);
threshold2 = min([value1(num2),value2(num2),value3(num2),value4(num2)]);