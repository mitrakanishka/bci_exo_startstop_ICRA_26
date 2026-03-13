% eeg_interp() - interpolate data channels
%
% Usage: EEGOUT = eeg_interp(EEG, badchans, method);
%
% Inputs: 
%     EEG      - EEGLAB dataset
%     badchans - [integer array] indices of channels to interpolate.
%                For instance, these channels might be bad.
%                [chanlocs structure] channel location structure containing
%                either locations of channels to interpolate or a full
%                channel structure (missing channels in the current 
%                dataset are interpolated).
% Output: 
%     EEGOUT   - data set with bad electrode data replaced by
%                interpolated data
%
% Author: Arnaud Delorme, CERCO, CNRS, Mai 2006-

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [eeg, interpFun] = eeg_interp(eeg, chanlocs, bad_elec)

goodlocs = chanlocs(~bad_elec);
xelec = [ goodlocs.X ];
yelec = [ goodlocs.Y ];
zelec = [ goodlocs.Z ];
rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
xelec = xelec./rad;
yelec = yelec./rad;
zelec = zelec./rad;
badlocs = chanlocs(bad_elec);
xbad = [ badlocs.X ];
ybad = [ badlocs.Y ];
zbad = [ badlocs.Z ];
rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
xbad = xbad./rad;
ybad = ybad./rad;
zbad = zbad./rad;

[C, G] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad);
interpFun = @(x) estimateBadChannels(C, G, x);
badchansdata = interpFun(eeg(:,~bad_elec));

eeg(:,bad_elec) = badchansdata;
end

function values = estimateBadChannels(C, G, eeg)
meanvalues = mean(eeg, 2);
eeg = bsxfun(@minus, eeg, meanvalues);
eeg = [eeg zeros(size(eeg,1),1)];
C = C * eeg';

% apply results
% -------------
values = zeros(size(eeg,1), size(G,1));
for j = 1:size(G,1)
    values(:,j) = sum(C .* repmat(G(j,:)', [1 size(C,2)]));        
end
values = bsxfun(@plus, values, meanvalues);

end

function [C, Gsph] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad)

Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
Gsph  = computeg(xbad,ybad,zbad,xelec,yelec,zelec);

% compute solution for parameters C
% ---------------------------------
C = pinv([Gelec;ones(1,length(Gelec))]);

end
% compute G function
% ------------------
function g = computeg(x,y,z,xelec,yelec,zelec)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +... 
                (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
                (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);

g = zeros(length(x(:)),length(xelec));
%dsafds
m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    L = legendre(n,EI);
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);    
end
