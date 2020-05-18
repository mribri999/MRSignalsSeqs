function [FC, M0S, M1S, t_ss, G_ss] = conventional_flowcomp(params)

    G_ss = horzcat(linspace(params.g_ss*1e-3,params.g_ss*1e-3,(ceil(params.p_ss*1e-3/params.dt))),linspace(params.g_ss*1e-3,params.g_ss*1e-3/10,10));
    t_ss = length(G_ss)*params.dt;
    t_vec = 0:params.dt:(length(G_ss)-1)*params.dt;
    
    M0S = cumsum((G_ss.*params.dt).*1e6,2);          % [mT/m x ms]
    M1S = cumsum((G_ss.*t_vec*params.dt).*1e9,2);    % [mT/m x ms^2]
    
    M0S_ = M0S(1,end);
    M1S_ = M1S(1,end);
    
    for r = 1e-03:1e-03:5
        r_ = (ceil(r/params.dt/1000));
        h = r*params.smax;
        M02 = (-h*r+sqrt((h*r)^2+2*(h*r*M0S_ + M0S_^2 + 2*h*M1S_)))/2;
        M01 = M02 + M0S_;
        w1 = M01/h + r;
        w2 = M02/h + r;
        w1_ = (ceil(w1/params.dt/1000));
        w2_ = (ceil(w2/params.dt/1000));

        if  (w1_-2*r_ <= 1) || (w2_-2*r_ <= 1)
            h1 = M01/(r_-w1_)*100;
            h2 = -M02/(r_-w2_)*100;
            break
        end
    end

    G = horzcat(linspace(0,h1,r_),linspace(h1,h1,w1_-2*r_),linspace(h1,h2,2*r_),linspace(h2,h2,w2_-2*r_),linspace(h2,0,r_));
    FC = horzcat(G_ss*1000,G);

end
