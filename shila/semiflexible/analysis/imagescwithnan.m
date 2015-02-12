function [h, hcb] = imagescwithnan(a, cm, nanclr, hidenan, clims)
% IMAGESC with NaNs assigning a specific color to NaNs
%   Usage example:
%       a = peaks; 
%       a(a < 0) = nan; 
%       imagescwithnan(a, hot, [0 1 1])               % [0 1 1] is cyan
%   or
%       imagescwithnan(a, hot, [0 1 1], true)         % hide NaN color from
%                                                     % the colorbar
%   or
%       imagescwithnan(a, hot, [0 1 1], false, [0 4]) % set CLims
%
% This is an improved version of the snippet written by yuk on
% StackOverflow (http://stackoverflow.com/a/8482139/862188)
% 
% Author:      yuk     (http://stackoverflow.com/users/163080/yuk)
% Modified by: Zertrin (http://stackoverflow.com/users/862188/zertrin)
% Last modification: 10 May 2013

%# find minimum and maximum
if exist('clims','var')
    amin = clims(1);
    amax = clims(2);
else
    amin = min(a(:));
    amax = max(a(:));
end
%# size of colormap
n = size(cm,1);
%# color step
dmap = (amax - amin)/n;

%# standard imagesc
if exist('clims','var')
    him = imagesc(a,clims);
else
    him = imagesc(a);
end
%# add nan color to colormap
colormap([nanclr; cm]);
%# changing color limits
caxis([amin-dmap amax]);
%# place a colorbar
hcb = colorbar;
%# change Y limit for colorbar to avoid showing NaN color
if exist('hidenan','var') && hidenan
    ylim(hcb,[amin amax])
end

if nargout > 0
    h = him;
end
