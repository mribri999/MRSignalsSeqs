function h = msinc(N,ncyc)
%
%	

x = ((([0:N-1]+0.5)-(N/2-0.5))/N)*2*ncyc;
hw = hamming(N);

h = sinc(x(:)).*hw(:);
h = h';

