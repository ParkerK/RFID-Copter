% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - modified DTFT phase estimator
%
% references:
% [1] H.C. So, A comparative study of two discrete-time phase delay estimators, IEEE Transactions on
%     Instrumentation and Measurement, Vol. 54, No. 6, December 2005, page 2501-2504
%
% IMPORTANT: the sampling frequency has to be a multiple of 4*sinusoid-frequency
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
% version = mdtft()
%    Just returns the version number (string).
% [phi, delta_err] = mdtft(x1, x2, settings)
%    Returns the phase shift between x1 and x2 in [rad] (x1 and x2 being sinusoids; sampling frequency
%    SETTINGS.FS, sinusoid frequency SETTINGS.f0). The function also returns the 90? phase shift 
%    rounding error in [samples] (see [1]).
%
%
% ***** Interface definition *****
% function [phi, delta_err] = mdtft(x1, x2, settings)
%    x1, x2      input signals (-> "phase difference")
%    settings    struct containing estimator settings
%       .fs           sampling frequency
%       .f0           frequency of sinusoids
%       .delta_warn   (optional; default 1e-6) a warning is thrown if the maximum 90 deg phase shift  
%                     rounding error [samples] exceeds this limit 
%       .delta_warn   (optional; default 10*eps) an error is thrown if the maximum 90 deg phase shift  
%                     rounding error [samples] exceeds this limit 
%   
%    phi         estimated phase difference between the sinusoids in X1 and X2
%    delta_err   90 deg phase shift rounding error in [samples]
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

function [phi, delta_err] = mdtft(x1, x2, settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   phi = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings
%     maximum error in delta rounding before a warning / an error is thrown (only if not defined in settings)
internalsettings.delta_warn = 10*eps;
internalsettings.delta_err  = 1e-6;


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fs', 'f0'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% complete settings if necessary: error bound for delta rounding
if ~isfield(settings, 'delta_warn')
   settings.delta_warn = internalsettings.delta_warn;
end
if ~isfield(settings, 'delta_err')
   settings.delta_err  = internalsettings.delta_err;
end

% check lengths
if length(x1) ~= length(x2)
   critwarn('Length of input vectors x1 and x2 does not match. Truncating to identical length.');
   if length(x1) > length(x2)
      x1 = x1(1:length(x2));
   else
      x2 = x2(1:length(x1));
   end
end


% *******************************************************************************************************
% estimator according to [1]

% calculate minimum delta 
settings.delta = settings.fs / (4 * settings.f0);
%     error created by rounding
delta_err = settings.delta - round(settings.delta);
%     check if (approximately) integer
if abs(delta_err) > settings.delta_err
   err('fs is not a multiple of 4*f0 (error = %.3e samples)', delta_err);
elseif abs(delta_err) > settings.delta_warn
   critwarn('fs is not a multiple of 4*f0 (error = %.3e samples)', delta_err);
end

%     round (used as index...)
settings.delta = round(settings.delta);

% analytic forms x(k) + j x(k-delta)
r1 = complex( x1(1+settings.delta:end), x1(1:end-settings.delta) );
r2 = complex( x2(1+settings.delta:end), x2(1:end-settings.delta) );

% zero-pad to make sure f0 is exactly on a bin
n = round( ceil(length(r1)*settings.f0/settings.fs)*settings.fs/settings.f0 ); % round just to make sure
r1 = [r1(:); zeros(n-length(r1),1)];
r2 = [r2(:); zeros(n-length(r2),1)];

% calculate the DFT
R1 = goertzel(r1, round(length(r1)*settings.f0/settings.fs + 1));
R2 = goertzel(r2, round(length(r2)*settings.f0/settings.fs + 1));

% phase estimate
phi = angle( R1 * conj(R2) );

