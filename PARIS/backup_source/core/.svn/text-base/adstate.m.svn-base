% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% assembly / detuning state <-> assembly / detuning values
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
% [as, ds] = adstate('settings->state', adsm, cat, rat_shift, enr, fsr)
%    Returns assembly state (AS) and detuning state (DS) for assembly CAT/RAT_SHIFT and detuning ENR/FSR,
%    resp., based on the assembly/detuning state mapping in ADSM.
% [cat, rat_shift, enr, fsr] = adstate('state->settings', adsm, as, ds)
%    Returns assembly settings CAT/RAT_SHIFT and detuning settings ENR/FSR for assembly state AS and
%    detuning state DS, resp., based on the assembly/detuning state mapping in ADSM.
%
%
% ***** Interface definition *****
%             [as, ds] = adstate(mode, adsm, cat, rat, enr, fsr) for mode='settings->state'
% [cat, rat, enr, fsr] = adstate(mode, adsm, as, ds)                   for mode='state->settings'
%    cat    assembly capacity in F (assembly)
%    rat    percent shift of assembly resistance from optimum (assembly)
%    enr    enhancement/boost of resonance (detuning)
%    fsr    frequency shift of resonance in Hz (detuning)
%    as     assembly state (index in assembly settings list)
%    ds     detuning state (index in detuning settings list)
%    mode   direction of mapping {'settings->state', 'state->settings'}
%    adsm   assembly/detuning state mapping (struct)
%       .assembly
%          .v_cat      list of assembly capacities in F (index is detuning state)
%          .v_rats     list of assembly resistance shifts in percent  (index is detuning state)
%       .detuning
%          .v_res_en   list of resonance boost values (index is detuning state)
%          .v_fshift   list of resonance frequency shifts in Hz  (index is detuning state)
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
% ? select nearest for close matches
%
% *******************************************************************************************************

function [out1, out2, out3, out4] = adstate(varargin)

% *******************************************************************************************************
% mode / input parameters

if nargin < 1
   error('Not enough input arguments.');
end

mode = lower(varargin{1});
switch mode
   case 'settings->state'
      if nargin < 6
         error('Not enough input arguments.');
      end
      adsm = varargin{2}; % assembly/detuning state matching
      cat  = varargin{3}; % assembly capacity
      rat  = varargin{4}; % assembly resistance SHIFT
      enr  = varargin{5}; % enhance/boost resonance
      fsr  = varargin{6}; % frequency shift of resonance
   case 'state->settings'
      if nargin < 4
         error('Not enough input arguments.');
      end
      adsm = varargin{2}; % assembly/detuning state matching
      as   = varargin{3}; % assembly state
      ds   = varargin{4}; % detuning state
   otherwise
      err('Unsupported mode "%s".', mode);
end


% *******************************************************************************************************
% map settings to state

if strcmpi(mode, 'settings->state')
   out1 = find( approxequal(adsm.assembly.v_cat,      cat) & approxequal(adsm.assembly.v_rats,   rat));
   out2 = find( approxequal(adsm.detuning.v_res_en,   enr) & approxequal(adsm.detuning.v_fshift, fsr));
   out3 = NaN;
   out4 = NaN;
end


% *******************************************************************************************************
% map state to settings

if strcmpi(mode, 'state->settings')
   out1 = adsm.assembly.v_cat(as);
   out2 = adsm.assembly.v_rats(as);
   out3 = adsm.detuning.v_res_en(ds);
   out4 = adsm.detuning.v_fshift(ds);
end


