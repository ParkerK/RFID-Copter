% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% interp2 for complex Z values, interpolating magnitude/phase and not re/im (like interp2 does)
%
%
%
% ***** Copyright / License / Authors *****
% Copyright 2007, 2008, 2009, 2010, 2011 Daniel Arnitz
%   Signal Processing and Speech Communication Laboratory, Graz University of Technology, Austria
%   NXP Semiconductors Austria GmbH Styria, Gratkorn, Austria
% Copyright 2012 Daniel Arnitz
%   Reynolds Lab, Department of Electrical and Computer Engineering, Duke University, USA
%
% This file is part of the PARIS Simulation Framework.
%
% The PARIS Simulation Framework is free software: you can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
%
% The PARIS Simulation Framework is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with the PARIS Simulation 
% Framework. If not, see <http://www.gnu.org/licenses/>.
%
%
% ***** Behavior *****
% Like interp2; interpolates magnitude/phase instead of real/imaginary part. The phase is unwrapped 
% along the first dimension. (default is method is 'linear', default extrapval is NaN)
%
%
% ***** Interface definition *****
% see Matlab function INTERP2; supports the following syntaxes
%    Zi = interp2(X, Y, Z, Xi, Yi)
%    Zi = interp2(X, Y, Z, Xi, Yi, method)
%    Zi = interp2(X, Y, Z, Xi, Yi, method, extrapval)
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
%
%
% ***** Todo *****
% ? check for remaining phase jumps in unwrapped Z; more intelligent unwrapping ("best dimension")
%
% *******************************************************************************************************

function Zi = interp2c(varargin)

% demux parameters
%     obligatory/default
X      = varargin{1};% interp2(X, Y, Z, Xi, Yi)
Y      = varargin{2};
Z      = varargin{3};
Xi     = varargin{4};
Yi     = varargin{5};
method = 'linear';
extrap = NaN;
%     optional 
if nargin == 6 % interp2(X, Y, Z, Xi, Yi, method)
   method = varargin{6};
elseif nargin == 7 % interp2(X, Y, Z, Xi, Yi, method, 'extrap'), interp2(X, Y, Z, Xi, Yi, method, extrapval)
   method = varargin{6};
   extrap = varargin{7};
end

% interpolate magnitude and phase of y separately
Zi_mag = interp2(X,Y,          abs(Z) , Xi,Yi, method,   abs(extrap));
Zi_arg = interp2(X,Y, unwrap(angle(Z)), Xi,Yi, method, angle(extrap)); % unwrapping will not work under all circumstances

% reassemble
Zi = Zi_mag .* exp(complex(0,Zi_arg));

