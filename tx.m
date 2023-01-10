function [ X,sym_pos] = tx( Ntx,sym_num,Constellations )
%TX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

order = log2(length(Constellations));
len = Ntx * sym_num;
pos = randi([1,2^order],1,len);
% Mtest
% pos = randi([1,4],1,len);
% pos = pos.^2;
% pos(pos == 9) = 4;
%Estest
% pos = randi([1,4],1,len);
% pos = pos + 5;
% pos(pos == 8) = 10;
% pos(pos == 9) = 11;
%
sym_pos = reshape(pos,Ntx,[]);
X = Constellations(sym_pos);
end

