function x = zeroz(sig)

%peaks
P2N = find(sig(1:end-1) > 0 & sig(2:end) < 0);
%mins
N2P = find(sig(1:end-1) < 0 & sig(2:end) > 0);
%disp(size(P2N))
%disp('gap')
%disp(size(N2P))
%combine
x = sort([P2N, N2P]);
end