% Double Gaussian reconstruction with adjusted signs for peaks

function reco = bigaus_recon(LM,RM,Map,Mip,sig)

reco = zeros(1,length(sig));

offset = mean(sig);
Map = Map - offset;
Mip = Mip + offset;

%index vector
x = linspace(1,length(sig),length(sig));

for k = 1:length(LM)

    %center of bead
    t = round((LM(k)+RM(k))/2);
    %pk to pk transient time
    o = RM(k) - LM(k);
    %width control
    d = 0.3 * o;
    

    %calculations
    ex_1 = Map(k)*exp(-((x-t+(o/2)).^2)/(2*d^2)); % Positive peak
    ex_2 = Mip(k)*exp(-((x-t-(o/2)).^2)/(2*d^2)); % Negative peak

    reco = reco + ex_1 + ex_2;

end