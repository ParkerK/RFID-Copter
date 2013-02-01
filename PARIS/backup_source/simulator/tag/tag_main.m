% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% tag - main
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
% version = tag_main()
%    Just returns the version number (string).
% settings = tag_main('initialize', settings)
%   complete and sanitize settings
% settings = tag_main('re-initialize', settings)
%   reset to unpowered state
% output = tag_main('rx', input, settings)
%   receive, check if powered
% output = tag_main('rx_data', input, settings)
%   demodulate and decode
% output = tag_main('set_epcstate', input, settings)
%   set EPCglobal Class-1 Gen-2 state (not a complete / fully correct impl. of the EPC state diagram!)
% settings = tag_main('query', settings)
%   set up return link
% settings = tag_main('clen_addoverhead', settings)
%   add modulation overhead to carrier length
% settings = tag_main('prep_rn16', settings)
%   prepare modulation of RN16
% output = tag_main('tx_data_active', input, settings)
%   modulate data, returns an empty vector if tag is inactive
% output = tag_main('tx_data_all', input, settings)
%   modulate data, also include an inactive tag (passive reflection)
%
%
% ***** Interface definition *****
% function output = tag_main(mode, input, settings)
%    cmd        string with command (mode of operation) ... see Behavior
%    input      input signal (for some commands)
%    settings   settings structs (one for each reader)
%
%    output     depending on command: settings struct, results struct, output signal
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


function output = tag_main(varargin)
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
% input parameters

% number of input arguments
switch nargin
   case 2 % send
      cmd      = varargin{1};
      settings = varargin{2};
   case 3 % receive
      cmd      = varargin{1};
      input    = varargin{2};
      settings = varargin{3};
   otherwise
      err('Wrong amount of input parameters.');
end


% *******************************************************************************************************
% switch commands
switch lower(cmd)
   
   % ****************************************************************************************************
   % MISC: complete and sanitize settings
   case 'initialize'
      % modulation (esp. group delay and t0)
      settings.modulation = tag_modulation([0], settings.modulation);
      % make sure tag is unpowered
      settings.state.powered = false;
      settings.state.epc     = '';
      if isfield(settings.state, 'linkinfo')
         settings.state = rmfield(settings.state, 'linkinfo');
      end
      % return modified settings
      output = settings;
      
   % ****************************************************************************************************
   % MISC: reset to unpowered state
   case 're-initialize'
      settings.state.powered = false;
      settings.state.epc     = '';
      if isfield(settings.state, 'linkinfo')
         settings.state = rmfield(settings.state, 'linkinfo');
      end
      % return modified settings
      output = settings;
      
   % ****************************************************************************************************
   % PASSIVE: receive, check if powered
   case 'rx'
      % calculate transmitted (non-reflected) signal
      [output.rxsignal, tagmod_internal] = tag_modulation(input, [], [], settings.modulation);
      % power levels VIA TAG_MODULATION UNTIL TAG_POWER IS IMPLEMENTED
      output.pav_hat = tagmod_internal.pav; % available power
      output.pin_hat = tagmod_internal.pin; % input power (not reflected)
      output.pic_hat = tagmod_internal.pic; % chip input power
      output.powered = tagmod_internal.func; % functional?
      output.vdda    = tagmod_internal.vdda; % power supply voltage
      %     modify vdda setting
      output.settings = settings;
      output.settings.demodulation.vdda = output.vdda;
   
   % ****************************************************************************************************
   % PASSIVE: demodulate and decode
   case 'rx_data'
      % create tag clock (sampling mode)
      settings.clock.mode = 'fs';
      settings.clock.length = length(input);
      clock = tag_clock(settings.clock);
      % demodulate and sample
      tag_demodulated = tag_demodulation(input, clock, settings.demodulation);
      % decode
      [output.decoded, output.linkinfo] = tag_decoding(tag_demodulated, settings.decoding);
      
   % ****************************************************************************************************
   % PASSIVE: set EPCglobal Class-1 Gen-2 state
   % ... WARNING: THIS IS NOT A COMPLETE/CORRECT IMPLEMENTATION OF THE STATE DIAGRAM
   case 'set_epcstate'
      % not powered
      if ~settings.state.powered
         settings.state.epc = '';
         output = settings;
         return
      end
      % partial implementation of EPC Class-1 Gen-2 state diagram
      switch lower(settings.state.epc)
         case '' % was not powered => ready
            settings.state.epc = 'ready';
         case 'ready'
            if strcmpi(settings.state.linkinfo.cmd, 'query') && settings.state.linkinfo.crc5_ok
               settings.state.epc = 'reply';
            end
         case 'reply'
         case 'arbitrate'
         otherwise
            err('Unsupported EPC state "%s" while settings new state.', lower(settings.state.epc));
      end
      % return modified settings
      output = settings;
      
   % ****************************************************************************************************
   % PASSIVE: set up return link
   case 'query'
      % set up return link if in reply or arbitrate state
      if strcmpi(settings.state.linkinfo.cmd, 'reply') || strcmpi(settings.state.linkinfo.cmd, 'arbitrate')
         % encoding
         settings.encoding.trext = settings.state.linkinfo.trext; % pilot tone on/off
         settings.encoding.m     = settings.state.linkinfo.m;     % (1:FM0, 2,4,8:Miller)
         % modulation
         settings.modulation.lf  = settings.state.linkinfo.dr / settings.state.linkinfo.trcal;
      end
      % return modified settings
      output = settings;
      
   % ****************************************************************************************************
   % ACTIVE (MISC): add modulation overhead to carrier length
   case 'clen_addoverhead'
      % add the overhead required by the tag implementation
      output = input + settings.modulation.grpdel_s + 2*settings.modulation.t0_s; % [samples]
      % make carrier length a multiple of largest filterbank size (required by tag_modulation)
      output = output + settings.modulation.nfft - mod(output, settings.modulation.nfft);
      
   % ****************************************************************************************************
   % ACTIVE: prepare modulation of RN16
   case 'prep_rn16'
      % setup data
      if strcmpi(settings.state.epc, 'reply') 
         % encode command and make alphabet binary {0,1} (tag_encoding returns {-1, 1})
         settings.encdata = (tag_encoding(hex2bit(settings.id.rn16), settings.encoding, false) + 1) / 2;
      else
         settings.encdata = [0]; % unmodulated backscatter
      end
      % determine lengths and sanitize settings
      settings.modulation = tag_modulation(settings.encdata, settings.modulation);
      % return modified settings
      output = settings;
      
%    % ****************************************************************************************************
%    % ACTIVE: prepare continuous modulation with LF
%    case 'prep_tone'
%       settings.encdata = se
%       % determine lengths and sanitize settings
%       settings.modulation = tag_modulation(settings.encdata, settings.modulation);
%       % return modified settings
%       output = settings;
      
   % ****************************************************************************************************
   % ACTIVE: modulate data, returns an empty vector if tag is inactive
   case 'tx_data_active'
      if strcmpi(settings.state.epc, 'reply')
         % create tag clock
         settings.clock.mode   = 'fclk';
         settings.clock.length = settings.modulation.length_fclk_s;
         tag_clk = tag_clock(settings.clock);
         % modulate
         [output.txsignal, tagmod_internal] = tag_modulation(input, settings.encdata, tag_clk, settings.modulation);
         % power levels VIA TAG_MODULATION UNTIL TAG_POWER IS IMPLEMENTED
         output.pav_hat = tagmod_internal.pav; % available power
         output.pin_hat = tagmod_internal.pin; % input power (not reflected)
         output.pic_hat = tagmod_internal.pic; % chip input power
         output.powered = tagmod_internal.func; % functional?
         output.vdda    = tagmod_internal.vdda; % power supply voltage
         %     modify vdda setting
         output.settings = settings;
         output.settings.demodulation.vdda = output.vdda;
         %     add some more values to output (for debugging purposes)
         output.rho  = tagmod_internal.rho;  % linear refl. coeff. model: center value
         output.drho = tagmod_internal.drho; % linear refl. coeff. model: difference value
      else
         % create fake values
         output.txsignal = [];
         output.pav_hat  = NaN;
         output.pin_hat  = NaN;
         output.pic_hat  = NaN;
         output.powered  = NaN;
         output.vdda     = NaN;
         output.rho      = NaN;
         output.drho     = NaN;
         output.settings = settings;
         output.settings.demodulation.vdda = NaN;
      end
      
   % ****************************************************************************************************
   % ACTIVE: modulate data, also include inactive tags (passive reflection)
   case 'tx_data_all'   
      % create tag clock
      settings.clock.mode   = 'fclk';
      settings.clock.length = settings.modulation.length_fclk_s;
      tag_clk = tag_clock(settings.clock);
      % modulate
      [output.txsignal, tagmod_internal] = tag_modulation(input, settings.encdata, tag_clk, settings.modulation);      
      % power levels VIA TAG_MODULATION UNTIL TAG_POWER IS IMPLEMENTED
      output.pav_hat = tagmod_internal.pav; % available power
      output.pin_hat = tagmod_internal.pin; % input power (not reflected)
      output.pic_hat = tagmod_internal.pic; % chip input power
      output.powered = tagmod_internal.func; % functional?
      output.vdda    = tagmod_internal.vdda; % power supply voltage
      %     modify vdda setting
      output.settings = settings;
      output.settings.demodulation.vdda = output.vdda;
      %     add some more values to output (for debugging purposes)
      output.rho  = tagmod_internal.rho;  % linear refl. coeff. model: center value
      output.drho = tagmod_internal.drho; % linear refl. coeff. model: difference value
     
   % ****************************************************************************************************  
   otherwise
      err('Unsupported command: "%s".', lower(cmd));
end
