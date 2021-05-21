function [FE, DeltaM1] = conventional_flowencode(params)

GAM = 2*pi*42.57;         % 1/(s*T)
DeltaM1 = pi/(GAM*params.VENC);  % mT/m x ms^2

M0S = params.M0S(1,end);
M1S = params.M1S(1,end);

for r = 0.01:(0.5-0.01)/490:0.5
  h = r*params.smax;
  M02 = abs(-h*r+sqrt((h*r)^2+2*(h*r*M0S + M0S^2 + 2*h*(M1S+DeltaM1))))/2;
  M01 = M0S + M02;
  w1 = abs(M01)/h + r;
  w2 = abs(M02)/h + r;
  r_ = (ceil(r/params.dt/1000));
  w1_ = (ceil(w1/params.dt/1000));
  w2_ = (ceil(w2/params.dt/1000));
  
  if  (w1_-2*r_ <= 1) || (w2_-2*r_ <= 1)
    break
  end
  
  G = horzcat(linspace(0,-h,r_),linspace(-h,-h,w1_-2*r_),linspace(-h,0,r_),linspace(h/r_,h,r_-1),linspace(h,h,w2_-2*r_),linspace(h,0,r_));
  FE = horzcat(params.G_ss*1000,G);
  
end

return