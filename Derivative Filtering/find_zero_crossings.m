function zero_crossings = find_zero_crossings(deriv)
    zero_crossings = [];
    for i = 1:length(deriv)-1
        if deriv(i) * deriv(i+1) < 0
            % Sign change indicates zero crossing
            zero_crossings = [zero_crossings, i];
        end
    end
end
