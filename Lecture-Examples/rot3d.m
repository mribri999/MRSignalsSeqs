% This function simply activates the camer toolbar and allows 3D rotation
%
% Uses: cameratoolbar('SetMode','orbit','SetCoordSys','none');
%
% DBE 12/29/02

function rot3d

  cameratoolbar('Show');
  cameratoolbar('SetMode','orbit');
  cameratoolbar('SetCoordSys','none');

return