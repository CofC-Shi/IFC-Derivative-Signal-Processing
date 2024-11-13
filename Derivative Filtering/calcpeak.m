function [Map,Mip] = calcpeak(LL,LM,RM,sig)

Map = zeros(1,length(LL));
Mip = zeros(1,length(LL));

for k = 1:length(LL)
    if(LM(k) == 0 || RM(k) == 0)
        continue
    end
    Map(k) = sum(sig(LL(k):LM(k)));
    Mip(k) = sum(sig(LL(k):RM(k)));
end