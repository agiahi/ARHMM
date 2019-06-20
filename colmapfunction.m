function farbe = colmapfunction(x,z2,brightness)
% ReLU(x) = x*(x>0)
if nargin==1, z2=rgb2hsv([0,0,1]);brightness=1;end
if nargin==2, brightness=1;end

x = min(max(x,-1),1); % clip to [-1,1]
angl2=z2; % discrete angles from discrete color palette
absz2=x; % saturation prop to responsibility

farbe = hsv2rgb([angl2,absz2,brightness-.25*x]);