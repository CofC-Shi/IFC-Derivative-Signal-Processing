function events = detect_ev(sig,plotz)
    threshold = 1e-6;
    events = [];
    eind = 1;
    dip = 0;
    event = 0;
    start_x = 0;

    for k = 1:length(sig)
        if(sig(k) > threshold && event == 0)
            event = 1;
            start_x = k;
        end
        if(sig(k) < -threshold && event == 1)
            dip = 1;
        end
        if(sig(k) > -threshold && (dip + event == 2))
            event = 0;
            dip = 0;
            end_x = k;
            events(eind) = start_x;
            events(eind+1) = end_x;
            %disp(['beads start ... end: ',num2str(start_x),' ... ',num2str(end_x)])
            eind = eind + 2;
        end      
    end


    if(plotz == 1)
        plot(sig,'Color','g','LineWidth',2)
        hold on;
        for k = 1:length(events)/2
        line([events(1+2*(k-1)),events(1+2*(k-1))], ylim, 'Color', 'r', 'LineWidth', 0.5);
        line([events(2*k),events(2*k)], ylim, 'Color', 'b', 'LineWidth', 0.5);
        end
    end
end
