% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% channel - noise (time-invariant)
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
% version = channel_large()
%    Just returns the version number (string).
% output = channel_large(input, settings)
%    Adds noise according to SETTINGS to INPUT and returns the result in OUTPUT. 
%    Note that SETTINGS has to be complete even if type is "no noise" (i.e., 'off'). 
%    The function is able to process multiple receivers at once: Each line of SETTINGS.N0 and 
%    SETTINGS.FRXS represents one receiver; pass only a single line to create identical settings 
%    for all receivers. The first  dimension (lines) of INPUT has to be "time", the second dimension 
%    (columns) has to be "receivers".
%       
%
%
% ***** Interface definition *****
% function output = channel_noise(input, settings)
%    input       input vector [time, receiver no.]
%    settings    struct containing all necessary parameters
%       .type       type of noise {'wgn' white Gaussian, otherwise: no noise}
%       .n0         noise density [dBm, per ? Hz], e.g. [-50, 1000] = -50 dBm/kHz (one line per receiver)
%       .frxs       receiver (reader or tag) sampling frequencyies in Hz (one line per receiver)
%       .fs         sampling frequency in Hz
%
%    output      output vector (same size as input vector)
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


function output = channel_noise(input, settings)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   output = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'type', 'n0', 'frxs', 'fs'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% check size of input, settings.no and settings.frxs
%     n0 and frxs (...problem for two receivers)
if size(settings.n0, 1) ~= size(settings.frxs, 1)
   err('Size of settings.n0 does not match size of settings.frxs (lines).')
end
%     compared to input if setup not identical for all receivers
if size(settings.n0, 1) ~= 1 && size(settings.n0, 1) ~= size(input, 2)
   err('Size of input does not match setup. Check size of input, settings.n0 and settings.frxs.');
end


% *******************************************************************************************************
% "the function"

% calculate noise power out of noise density [dBm, per ? Hz] (noise power matching at receiver)
settings.varn = 10.^(settings.n0(:,1)/10-3)./settings.n0(:,2) * settings.fs ./ settings.frxs;
settings.varn = settings.varn(:)';

% setup identical for all receivers (?)
if length(settings.varn) == 1 && size(input,2) > 1 % identical settings for all receivers
   settings.varn = settings.varn * ones(1,size(input,2));
end
  
% add noise
switch lower(settings.type)
   % white Gaussian noise
   case 'wgn'
      output = input + randn(size(input)) .* repmat(sqrt(settings.varn), size(input,1), 1);
      
   % no noise   
   otherwise
      output = input;
end

