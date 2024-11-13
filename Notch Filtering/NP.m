% find the peaks and return timestamp (peak indices) and peak amplitude

function [x,y] = NP(TH,sig)
    
    P = min(sig);
    
    sigD = BDC(sig);
    
    %mins
    N2P = find(sigD(1:end-1) < 0 & sigD(2:end) > 0);
    
    %thresholding
    x = [];
    y = [];
    
    count = 0;
    for k = 1:length(N2P)
    
        %ensure safe distance
        if(k > 1 && count > 0)
            if(N2P(k) - x(count) < 30)
                if(sig(N2P(k)) < sig(x(count)))
                x(count) = N2P(k);
                y(count) = sig(N2P(k));
                end
                continue
            end
        end
    
        %if min is large enough, add to list
        if(sig(N2P(k)) < TH*P)
            count = count + 1;
            x(end+1) = N2P(k);
            y(end+1) = sig(N2P(k));
        end
    end
end



