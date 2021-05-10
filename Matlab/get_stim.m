function stim_out = get_stim( G, dt )
%GET_STIM Summary of this function goes here
%   Detailed explanation goes here

    alpha = 0.333;
    r = 20.0;
    c = 334e-6;
    Smin = r/alpha;
    coeff = zeros(numel(G), 1);
    for i = 1:(numel(G))
        coeff(i) = ( c / ((c + dt*(numel(G)-1) - dt*(i-1))^2.0) / Smin );
    end

    stim_out = zeros(numel(G), 1);
    for j = 1:(numel(G)-2)
%         fprintf('%d \n', j);
        ss = 0;
        for i = 1:(j+1)
%             fprintf('%d  %d  %d  %d  %d\n', j, (numel(G)-1), i, (j+1), numel(coeff)-j+i-1);
%             fprintf('-----\n');
            ss = ss + coeff( numel(coeff)-j+i-1 ) * (G(i+1)-G(i));
        end
        stim_out(j) = ss;
    end

end

