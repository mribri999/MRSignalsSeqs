%
%	Generate 1-channel k-space with complex Gaussian noise.
%
if (exist('nscale')~=1)
  nscale = 0.001;
end;

im = circ(256,100);
ksp = ifftshift(ifft2(ifftshift(im)));
for n=1:1000
  kspace(:,:,n) = ksp + i*nscale*randn(256,256) + nscale*randn(256,256); 
end;

