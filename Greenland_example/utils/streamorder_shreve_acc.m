function s = streamorder_shreve_acc(St,Fpp)
% s = shreve_streamorder(St,Fpp)

s = streampoi(St,'channelheads','ix');
facc=Fpp.Z(s);

s = streampoi(St,'channelheads','logical');
s = double(s);
s(s==1)=facc;

for r = 1:numel(St.ix)
    s(St.ixc(r)) = s(St.ix(r))+s(St.ixc(r));
end
end
