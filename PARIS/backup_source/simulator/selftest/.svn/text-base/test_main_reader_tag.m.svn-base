% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% (quick) selftest: reader_main and tag_main
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
% version = test_main_reader_tag()
%    Just returns the version number (string).
% sumoferrors = test_main_reader_tag(sumoferrors)
%    Tests the functions reader_main and tag_main (as a team), displays the results in the command window
%    and returns the overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_main_reader_tag(sumoferrors);
%    sumoferrors    sum of errors found not including this fcn call
%
%    sumoferrors    sum of errors found including this fcn call (+1 for each erroneous tested
%                   functionality)
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
% - data check for T->R when decoding is implemented
%
% *******************************************************************************************************

function sumoferrors = test_main_reader_tag(sumoferrors)
version = 'beta 3.0';


% ******************************************************************************************************
% version system

% just return version number
if nargin == 0
   sumoferrors = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% tolerances
internalsettings.tol_t = 0.05; % rel. err. tolerance for timings (tari, rtcal, trcal, ...) 
internalsettings.tol_d = 0.05; % rel. err. tolerance for distance estimates

% minimum modulation depth before a note is printed in case of an error
internalsettings.min_modd = 0.01;


% *******************************************************************************************************
% initialization

% output
disp('   = reader_main and tag_main (top-level functionality) ***');

% load settings and expected results
data = loadmat('results_main_reader_tag', '     ');
settings     = data.settings;
partitioning = data.partitioning;


% *******************************************************************************************************
% run tests and evaluate

for p = 1 : length(partitioning.names)
   
   % once for each new test block
   %     output
   disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(p), partitioning.names{p}));
   %     reset errors and errortext
   errors    = 0;
   errortext = '';
   %     reset recorded values (per block)
   record.modd = [];
      
   % test-specific messages
   switch lower(partitioning.names{p})
      case 'data r->t'
         disp(sprintf('        checking sent data for equality and decoded timings (tari, rtcal, and trcal) +/-%g%%',...
            internalsettings.tol_t*100));
      case 'data t->r'
         disp(sprintf('        checking type of encoding and presence of pilot tone; checking LF +/-%g%% + nfft',...
            internalsettings.tol_t*100));
         disp('        Warning: no checks of sent data (reader_decoding not yet implemented)');
      case 'mfcw'
         disp(sprintf('        checking distance +/-%g%% (compensating expected tag phase shift)',...
            internalsettings.tol_d*100));
      otherwise
         err('Unrecognized test "%s". Check partitioning.', partitioning.names{p});
   end   
   
   for j = 1 + partitioning.indices(p) : partitioning.indices(p+1)
            
      % **************************************************
      % INITIALIZATION
      
      settings{j}.tag.modulation.charfile = 'tagchar_modulator_synth01';
      
      
      % stored result values
      results{j}.modd = NaN; % tag modulation depth (only if tag is transmitter)
      
      % initialize reader and tag (complete settings)
      % ... IMPORTANT: No sanitation (modification) of parameters here!
      settings{j}.reader = reader_main('initialize', settings{j}.reader);      
      settings{j}.tag = tag_main('initialize', settings{j}.tag);
      
      % if the reader has to send in this test: prepare command and set neede carrier length
      if strcmpi(partitioning.names{p}, 'data R->T')
         settings{j}.reader = reader_main('prep_query', settings{j}.reader);
         % add overhead required by tag modulation to get carrier length
         settings{j}.reader.max_clen =...
            tag_main('clen_addoverhead', settings{j}.reader.modulation.length_s, settings{j}.tag) / settings{j}.rand.fs;
      end
      
      % if the tag has to send in this test: configure state and uplink, set needed carrier lengths
      if strcmpi(partitioning.names{p}, 'data T->R') || strcmpi(partitioning.names{p}, 'MFCW')
            % state
            settings{j}.tag.state.powered = 1; % powered ...
            settings{j}.tag.state.epc     = 'reply'; % ... and ready to reply
            % incomplete state (only necessary fields) for uplink (always configure to reader 1)
            settings{j}.tag.state.linkinfo.cmd   = 'reply';
            settings{j}.tag.state.linkinfo.tari  =  settings{j}.reader.modulation.tari;
            settings{j}.tag.state.linkinfo.rtcal =  settings{j}.reader.modulation.rtcal;
            settings{j}.tag.state.linkinfo.trcal =  settings{j}.reader.modulation.trcal;
            settings{j}.tag.state.linkinfo.m     =  settings{j}.reader.command.m;
            settings{j}.tag.state.linkinfo.dr    =  settings{j}.reader.command.dr;
            % switch on pilot tones for ranging tests and set a "zero bit RN16" for speed
            if strcmpi(partitioning.names{p}, 'MFCW')
               settings{j}.tag.state.linkinfo.trext = true;
               settings{j}.tag.id.rn16              = '';
            else
               settings{j}.tag.state.linkinfo.trext = false;
            end
         % set these values, ...
         settings{j}.tag = tag_main('query', settings{j}.tag);
         % ... encode and ...
         settings{j}.tag = tag_main('prep_rn16', settings{j}.tag);
         % ... get needed carrier length
         settings{j}.reader.max_clen = settings{j}.tag.modulation.length;
      end        
      
      % set the needed carrier length found above
      settings{j}.reader = reader_main('set_maxclen', settings{j}.reader);
      
       
      % **************************************************
      % TEST: query command
      if strcmpi(partitioning.names{p}, 'data R->T')
         % create modulated carrier
         reader_carrier = reader_main('tx_data', settings{j}.reader);
         % channel R -> T
         tag_carrier = channel_main({reader_carrier}, settings{j}.channel_rt);
         % demodulate, decode, determine state
         %     re-initialize tag (set to unpowered state)
         settings{j}.tag = tag_main('re-initialize', settings{j}.tag);
         %     receive, check if powered
         temp = tag_main('rx', tag_carrier{:}, settings{j}.tag);
         settings{j}.tag = temp.settings; % modifications to vdda
         tag_rxsignal = temp.rxsignal;
         settings{j}.tag.state.powered = temp.powered;
         %     set EPC state accordingly ('' or 'ready')
         settings{j}.tag = tag_main('set_epcstate', settings{j}.tag);
         %     demodulate and decode
         temp = tag_main('rx_data', tag_rxsignal, settings{j}.tag);
         tag_decoded = temp.decoded;
         settings{j}.tag.state.linkinfo = temp.linkinfo;
         %     determine state of tag
         settings{j}.tag = tag_main('set_epcstate', settings{j}.tag);
         %     setup return link (if possible)
         settings{j}.tag = tag_main('query', settings{j}.tag);
         
         % checks
         [errors, errortext] = check_result(j,...
            settings{j}.tag.state.linkinfo.tari, settings{j}.reader.modulation.tari,...
            errors, errortext, 'tari', 'relerr', internalsettings.tol_t);
         [errors, errortext] = check_result(j,...
            settings{j}.tag.state.linkinfo.rtcal, settings{j}.reader.modulation.rtcal,...
            errors, errortext, 'rtcal', 'relerr', internalsettings.tol_t);
         [errors, errortext] = check_result(j,...
            settings{j}.tag.state.linkinfo.trcal, settings{j}.reader.modulation.trcal,...
            errors, errortext, 'trcal', 'relerr', internalsettings.tol_t);
         [errors, errortext] = check_result(j, tag_decoded, settings{j}.reader.data,...
            errors, errortext, 'data', 'equal');
      end
      
      
      % **************************************************
      % TEST: ack reply
      % ... this is not fully multi-reader-multi-tag compliant; create settings accordingly
      if strcmpi(partitioning.names{p}, 'data T->R')
         % create unmodulated carriers
         reader_carrier = reader_main('tx_carrier', settings{j}.reader);
         % channels R -> T
         tag_carrier = channel_main({reader_carrier}, settings{j}.channel_rt);
         % set state and modulate
         tx_temp = tag_main('tx_data_active', tag_carrier{:}, settings{j}.tag);
         settings{j}.tag = tx_temp.settings; % modifications to vdda
         tag_modcarrier  = tx_temp.txsignal; % signal
         % channels T -> R
         reader_modcarrier = channel_main({tag_modcarrier}, settings{j}.channel_tr);
         % demodulation and decoding
         reader_baseband = reader_main('rx', reader_modcarrier{:}, settings{j}.reader);
         % cut leadin/leadout t0 (filter transients) and decode (TEMPORARY: use linktiming_tag to decode timings)
         t0_s = round(settings{j}.reader.demodulation.frs * settings{j}.tag.modulation.t0 * 0.75);
         timings = linktiming_tag(abs(reader_baseband(1+t0_s:end-t0_s)));
         
         % check modulation depth to detect unlucky setup hitting zero modulation depth
         mod = abs( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1),...
            tx_temp.rho + tx_temp.drho, settings{j}.tag.modulation.fc, 'linear') );
         umd = abs( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1),...
            tx_temp.rho - tx_temp.drho, settings{j}.tag.modulation.fc, 'linear') );
         record.modd = [record.modd; (mod(:) - umd(:)) ./ (mod(:) + umd(:))];                    

         % checks 
         [errors, errortext] = check_result(j,...
            timings.tlf_s, settings{j}.reader.demodulation.frs/settings{j}.reader.modulation.lf,...
            errors, errortext, 'lf', 'mse',...
            (internalsettings.tol_t*settings{j}.reader.demodulation.frs/settings{j}.reader.modulation.lf...
            + settings{j}.tag.modulation.nfft*settings{j}.reader.demodulation.frs/settings{j}.rand.fs)^2);
         [errors, errortext] = check_result(j, timings.enc_m, settings{j}.reader.command.m,...
            errors, errortext, 'm', 'equal');
         [errors, errortext] = check_result(j, timings.enc_trext, 0,... % no pilot for these tests
            errors, errortext, 'trext', 'equal');
      end
      
      
      % **************************************************
      % TEST: MFCW
      % ... this is not fully multi-reader-multi-tag compliant; create settings accordingly
      if strcmpi(partitioning.names{p}, 'MFCW')
         % create unmodulated carriers (including MFCW secondary carriers)
         %     setup MFCW system
         settings{j}.reader = reader_main('prep_mfcw', settings{j}.reader);
         %     set length to maximum carrier length
         settings{j}.reader = reader_main('set_maxclen', settings{j}.reader);
         %     create unmodulated MFCW carriers
         reader_carrier = reader_main('tx_mfcw', settings{j}.reader);
         % channels R -> T
         tag_carrier = channel_main({reader_carrier}, settings{j}.channel_rt);
         % set state and modulate
         tx_temp = tag_main('tx_data_active', tag_carrier{:}, settings{j}.tag);
         settings{j}.tag = tx_temp.settings; % modifications to vdda
         tag_modcarrier  = tx_temp.txsignal; % signal
         % channels T -> R
         reader_modcarrier = channel_main({tag_modcarrier}, settings{j}.channel_tr);
         % demodulation and distance estimation
         reader_baseband = reader_main('rx', reader_modcarrier{:}, settings{j}.reader);
         dist_est = reader_main('mfcw_est', reader_baseband, settings{j}.reader);
         
         % check modulation depths to detect unlucky setup hitting zero modulation depth
         fi = settings{j}.tag.modulation.fc + [0, settings{j}.reader.ranging.mfcw_addseccarriers.fi];
         mod = abs( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1),...
            tx_temp.rho + tx_temp.drho, fi, 'linear') );
         umd = abs( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1),...
            tx_temp.rho - tx_temp.drho, fi, 'linear') );
         record.modd = [record.modd; (mod(:) - umd(:)) ./ (mod(:) + umd(:))];              

         % get theoretical values
         settings{j}.reader.mfcw_theorysim.t_max           = diff(settings{j}.reader.ranging.mfcw_addseccarriers.range_s) / settings{j}.rand.fs;
         settings{j}.reader.mfcw_theorysim.tag.tagchar_mod = settings{j}.tag.modulation.charfile;
         settings{j}.reader.mfcw_theorysim.dist            = settings{j}.channel_rt{1,1}.largescale.dist;
         dist_exp = mfcw_theorysim(settings{j}.reader.mfcw_theorysim);
         dist_exp = [dist_exp.f2(2), dist_exp.f3(1), dist_exp.f4(1), dist_exp.f5(1), dist_exp.f6(1)]; % serialize; 2-6FCW, no model
         
         % checks
         for nc = 1 : settings{j}.reader.ranging.mfcw_calcdist.nc
            [errors, errortext] = check_result(j, dist_est.dist{nc+1}, dist_exp(nc),...
               errors, errortext, sprintf('%iFCW',nc+1), 'relerr', internalsettings.tol_d);
         end
      end   
   end
   
   % final output for this test block
   if errors == 0
      disp('         ... passed');
   else
      disp(sprintf('         ... ERRORS: %s', errortext));
      sumoferrors = sumoferrors + 1;
      
      % could errors have been caused by unlucky setup?
      if min(abs(record.modd)) < internalsettings.min_modd
         disp(sprintf('         Note that these errors could have been caused by low modulation depth (min recorded: %.g).',...
            min(abs(record.modd)) ));
         disp('         Check and recreate settings if this is the case (unlucky setup hitting |mod|=|unmod|).');
      end     
   end
end


return
% *******************************************************************************************************
% debugging snippets

% if errors > 0
%    pause(0.01);
% end

% clc; test_main_reader_tag(0);

% [dist_exp(1:settings{j}.reader.ranging.mfcw_calcdist.nc)', cell2mat(dist_est)']

% est = cell2mat(dist_est);
% exp = dist_exp(1:settings{j}.reader.ranging.mfcw_calcdist.nc);
% (est - exp) ./ exp

% figure; hold on;
% plot(tag_carrier(1e4:2e4), 'b')
% plot(tag_modcarrier(1e4:2e4), 'r')
% hold off; grid on;

% fi = settings{j}.tag.modulation.fc + [0, settings{j}.reader.ranging.mfcw_addseccarriers.fi];
% abs( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1), tx_temp.rho - tx_temp.drho, fi, 'linear') )
% abs( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1), tx_temp.rho + tx_temp.drho, fi, 'linear') )
% angle( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1), tx_temp.rho - tx_temp.drho, fi, 'linear') )
% angle( interp1c(linspace(0,settings{j}.rand.fs/2,settings{j}.tag.modulation.nfft/2+1), tx_temp.rho + tx_temp.drho, fi, 'linear') )
