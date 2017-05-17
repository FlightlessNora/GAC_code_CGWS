%% ��ʵͼ��deblurʱ��Ҫ��ǰ���õĲ�����
%% ע����Щ��Ҫ���ⲿ���úõĲ�����
function [rotpadx,rotpady,maxx,maxy] = general_deblur_para_homo(grid_x,grid_y,grid_z,length_add,cx,cy)
maxz = max(abs(max(grid_z)),abs(min(grid_z)))*pi/180;
% maxzΪz����ת�����Ƕȣ������򶼰������Դ����ж���Ҫ�Ӷ��ı߽�
maxx = ceil(max(abs(grid_x)))+1;
maxy = ceil(max(abs(grid_y)))+1;
l1 = cx-cy*tan(maxz/2);
l2 = tan(maxz)*l1;
k2 = cy-cx*tan(maxz/2);
k1 = tan(maxz)*k2;
% �ӱ߽�Ŀ��������ֽ�Ͻ������Ƶ������¼ӵı߿�Ϊk_1*k_2/sqrt(k_1^2+k_2^2)
% k_1��k_2�ֱ�Ϊ2,4���ޱ���ȥ�������ε�ֱ�Ǳ�
% ����������������length_add��Ϊ�˵õ��б߽�������ģ��ͼ��
rotpadx = ceil(k1*k2/sqrt(k1^2+k2^2))+1+length_add(1);
% rotpadx = ceil(k1*k2/sqrt(k1^2+k2^2));
% ���Ҽӵı߿�Ϊl_1*l_2/sqrt(l_1^2+l_2^2)
% l_1��l_2�ֱ�Ϊ1,3���ޱ���ȥ�������ε�ֱ�Ǳ�
rotpady = ceil(l1*l2/sqrt(l1^2+l2^2))+1+length_add(2);
end