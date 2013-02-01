% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% interp1 for complex Y values, interpolating magnitude/phase and not re/im (like interp1 does)
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
% Like interp1; interpolates magnitude/phase instead of real/imaginary part
% (default is method is 'linear', default extrapval is NaN)
%
%
% ***** Interface definition *****
% see Matlab function INTERP1; supports the following syntaxes
%    yi = interp1(x, y, xi)
%    yi = interp1(x, y, xi, method)
%    yi = interp1(x, y, xi, method, 'extrap')
%    yi = interp1(x, y, xi, method, extrapval)
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
%
% *******************************************************************************************************

function yi = interp1c(varargin)

% demux parameters
%     obligatory/default
x      = varargin{1};% interp1(x, y, xi)
y      = varargin{2};
xi     = varargin{3};
method = 'linear';
extrap = NaN;
%     optional 
if nargin == 4% interp1(x, y, xi, method)
   method = varargin{4};
elseif nargin == 5% interp1(x, y, xi, method, 'extrap'), interp1(x, y, xi, method, extrapval)
   method = varargin{4};
   extrap = varargin{5};
end

% interpolate magnitude and phase of y separately
if isnumeric(extrap) % extrap = extrapval
   yi_mag = interp1(x,          abs(y) , xi, method,   abs(extrap));
   yi_arg = interp1(x, unwrap(angle(y)), xi, method, angle(extrap));
else % extrap = 'extrap'
   yi_mag = interp1(x,          abs(y) , xi, method, extrap);
   yi_arg = interp1(x, unwrap(angle(y)), xi, method, extrap);
end

% reassemble
yi = yi_mag .* exp(complex(0,yi_arg));

