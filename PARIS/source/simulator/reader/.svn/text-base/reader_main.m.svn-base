% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% reader - main
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
% version = reader_main()
%    Just returns the version number (string).
% settings = reader_main('initialize', settings)
%   complete and sanitize settings
% settings = reader_main('prep_query', settings)
%   prepare modulation of a query command  
% settings = reader_main('set_maxclen', settings)
%   set carrier length to maximum length
% settings = reader_main('prep_ranging', settings)
%   determine time frame for ranging operations (pilot tone)
% settings = reader_main('prep_mfcw', settings)
%   determine time frame for MFCW secondary carriers
% settings = reader_main('prep_fmcw', settings)
%   determine time frame for FMCW secondary carrier
% output = reader_main('tx_carrier', settings)
%   create and send an unmodulated carrier
% output = reader_main('tx_data', settings)
%   modulate and send prepared command
% output = reader_main('tx_mfcw', settings)
%   create and send multi-frequency continuous-wave carriers
% output = reader_main('tx_fmcw', settings)
%   create and send frequency-modulation continuous-wave carrier (UNDER CONSTRUCTION)
% output = reader_main('rx', input, settings)
%   receive and demodulate signal
% output = reader_main('mfcw_est', input, settings)
%    multi-frequency continuous-wave range estimation
% output = reader_main('fmcw_est(xcorr)', input, settings)
%    frequency-modulation continuous-wave range estimation using crosscorrelation (UNDER CONSTRUCTION)
%
%
% ***** Interface definition *****
% function output = reader_main(cmd, input, settings)
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
% ... ongoing FMCW implementation
% - unified command for ranging window (only once per simulation)
% ? more intelligent estimation window selection for ranging
%
% *******************************************************************************************************


function output = reader_main(varargin)
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
% internal settings

% maximum attenuation of components of interest (max/min) by reader input stage before a warning is issued
internalsettings.rxatt_max = 16; % dB (6 dB: power splitter!)


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
      % modulation
      settings.modulation = reader_modulation([0,0], settings.modulation); % even length for PR-ASK
      %     distribute modified settings
      settings.command.dr = settings.modulation.dr;
      % return modified settings
      output = settings;
   
   % ****************************************************************************************************
   % ACTIVE: prepare modulation of a query command
   case 'prep_query'
      %  encode command
      settings.data = reader_command('query', settings.command);
      % get needed length of carrier to encode data and sanitize settings
      settings.modulation = reader_modulation(settings.data, settings.modulation);
      settings.oscillator.length = settings.modulation.length;
      % return settings
      output = settings;

   % ****************************************************************************************************
   % ACTIVE (MISC): set carrier length to maximum length
   case 'set_maxclen'
      settings.modulation.length = settings.max_clen;
      settings.oscillator.length = settings.max_clen;
      output = settings;
      
   % ****************************************************************************************************
   % ACTIVE: determine time frame for ranging operations (pilot tone)
   case 'prep_ranging'
      % no pilot tone => take preamble
      if ~settings.command.trext
         warn('Pilote tone turned off; ranging might produce inaccurate results.');
         if settings.command.m == 1 % FM0
            settings.ranging.range_s = []; % entire signal
         else % Miller
            settings.ranging.range_s = ...
               [1, 4*settings.command.m*settings.demodulation.fs/settings.modulation.lf];
         end
      % select pilot tone as window
      else
         if settings.command.m == 1 % FM0
            settings.ranging.range_s = ...
               [1, 12*settings.demodulation.fs/settings.modulation.lf];
         else % Miller
            settings.ranging.range_s = ...
               [1, 16*settings.command.m*settings.demodulation.fs/settings.modulation.lf];
         end
      end
      % make sure the range does not exceed maximum carrier length
      % ... this will happen if no functional tags are left
      if settings.ranging.range_s(2) > settings.max_clen * settings.demodulation.fs
         warn('Truncating range to carrier length. Likely there are no replying tags left.');
         settings.ranging.range_s(2) = settings.max_clen * settings.demodulation.fs;
      end
      % return settings
      output = settings;
   
   % ****************************************************************************************************
   % ACTIVE: determine time frame for MFCW secondary carriers
   case 'prep_mfcw'
      % get position of pilot tone
      settings = reader_main('prep_ranging', settings);
      % and propagate the result to MFCW functions
      settings.ranging.mfcw_addseccarriers.range_s = settings.ranging.range_s;
      % return settings
      output = settings;
      
   % ****************************************************************************************************
   % ACTIVE: determine time frame for FMCW secondary carrier
   case 'prep_fmcw'
      % get position of pilot tone
      settings = reader_main('prep_ranging', settings);
      % and propagate the result to MFCW functions
      settings.ranging.fmcw_addfmcarrier.range_s = settings.ranging.range_s;
      % return settings
      output = settings;
      
   % ****************************************************************************************************
   % ACTIVE: create and send an unmodulated carrier
   case 'tx_carrier'
      % create carrier for modulation
      output = reader_oscillator(settings.oscillator);
      % transmitter
      output = reader_transmitter(output, settings.transmitter);

   % ****************************************************************************************************
   % ACTIVE: modulate and send prepared command
   case 'tx_data'
      % create carrier for modulation
      output = reader_oscillator(settings.oscillator);
      % modulate
      output = reader_modulation(output, settings.data, settings.modulation);
      % transmitter
      output = reader_transmitter(output, settings.transmitter);

   % ****************************************************************************************************
   % ACTIVE: create and send multi-frequency continuous-wave carriers
   case 'tx_mfcw'
      % create main carrier
      output = reader_oscillator(settings.oscillator);
      % add secondary carriers
      output = mfcw_addseccarriers(output, settings.ranging.mfcw_addseccarriers);
      % transmitter
      output = reader_transmitter(output, settings.transmitter);
      
   % ****************************************************************************************************
   % ACTIVE: create and send frequency-modulation continuous-wave carrier
   case 'tx_fmcw'
      % create main carrier
      output = reader_oscillator(settings.oscillator);
      % create and add FM carrier
      output = fmcw_addfmcarrier(output, settings.ranging.fmcw_addfmcarrier);
      % apply transmitter to sent signal
      output = reader_transmitter(output, settings.transmitter);

   % ****************************************************************************************************
   % PASSIVE: receive and demodulate signal
   case 'rx'
      % reader input stage (bandpass and power splitter)
      reader_rx = reader_receiver(input, settings.receiver);
      % IQ demodulation and sampling
      output = reader_demodulation(reader_rx, settings.demodulation, settings.oscillator); 

   % ****************************************************************************************************
   % PASSIVE: multi-frequency continuous-wave range estimation
   case 'mfcw_est'
      % select estimation window
      if ~isempty(settings.ranging.mfcw_addseccarriers)
         range_rs = round(settings.ranging.mfcw_addseccarriers.range_s * settings.demodulation.frs / settings.demodulation.fs);
      else % empty => entire signal
         range_rs = [1, size(input, 1)];
      end
      %     stationary parts only
      ind = [round(range_rs(1) + settings.ranging.window(1)*(range_rs(2)-range_rs(1))) :...
             round(range_rs(1) + settings.ranging.window(2)*(range_rs(2)-range_rs(1)))];
      %     output the window size
      msg('Estimation window: %.2f us (%.0f samples)', length(ind)/settings.demodulation.frs*1e6, length(ind));

      % estimate gains of analog input path of receiver
      apfreqres = reader_analogpathest(settings.ranging.freq, settings.analogpathest,...
         settings.transmitter, settings.receiver, settings.demodulation);
      %     split (easier to handle)
      nc = settings.ranging.mfcw_compsel.nc; % number of secondary carriers
      g_mi = apfreqres.all([1, 4+0*nc:3+1*nc]);
      g_i  = apfreqres.all([2, 4+1*nc:3+2*nc]);
      g_im = apfreqres.all([3, 4+2*nc:3+3*nc]);
      
      % warning for high attenuation
      if any(20*log10(max(abs([g_mi, g_i, g_im]))/min(abs([g_mi, g_i, g_im]))) > internalsettings.rxatt_max) 
         critwarn('Considerable attenuation of MFWC components detected (%.2g dB). Check reader input stage.',...
            20*log10(max(abs([g_mi, g_i, g_im]))/min(abs([g_mi, g_i, g_im]))));
      else
         msg('Maximum attenuation of MFWC components: %.2g dB',...
            20*log10(max(abs([g_mi, g_i, g_im]))/min(abs([g_mi, g_i, g_im]))));
      end

      % component selection 
      [c_mi, c_i, c_im] = mfcw_compsel(input, settings.ranging.mfcw_compsel);
      
      % average and remove estimates of systematic errors
      if any(any(isnan(c_mi))) || any(any(isnan(c_i))) || any(any(isnan(c_im)))
         output.avg_c_mi = c_mi;
         output.avg_c_i  = c_i;
         output.avg_c_im = c_im;
      else
         output.avg_c_mi = mean(c_mi(ind,:)) ./ g_mi;
         output.avg_c_i  = mean(c_i (ind,:)) ./ g_i;
         output.avg_c_im = mean(c_im(ind,:)) ./ g_im;
      end

      % distance estimates
      output.dist = mfcw_calcdist(output.avg_c_mi, output.avg_c_im, settings.ranging.mfcw_calcdist);
      
   % ****************************************************************************************************
   % PASSIVE: frequency-modulation continuous-wave range estimation using cross-correlation
   case 'fmcw_est(xcorr)'
      % select estimation window
      if ~isempty(settings.ranging.fmcw_addfmcarrier.range_s)
         range_rs = round(settings.ranging.fmcw_addfmcarrier.range_s * settings.demodulation.frs / settings.demodulation.fs);
      else % empty => entire signal
         range_rs = [1, size(input.rx, 1)];
      end
      %     stationary parts only
      ind = [round(range_rs(1) + settings.ranging.window(1)*(range_rs(2)-range_rs(1))) :...
             round(range_rs(1) + settings.ranging.window(2)*(range_rs(2)-range_rs(1)))];
      %     output the window size
      msg('Estimation window: %.2f us (%.0f samples)', length(ind)/settings.demodulation.frs*1e6, length(ind));
      
      % calculate maximum correlation lag
      settings.ranging.xcorr_maxlag = ceil(settings.ranging.maxdist * settings.demodulation.frs / settings.ranging.c);
      
      % interpolate for better resolution (reader sampling resolution insufficient)
      switch(lower(settings.ranging.xcfmode))
         
         % interpolate signals (slowest, high RAM usage)
         case 'sig'
            % calculate oversampling rate and output resulting resolution(s)
            settings.ranging.xcorr_osr = ceil( settings.ranging.c / (settings.ranging.res * settings.demodulation.frs) );  
            msg('FMCW xcorr resolution original / interp [cm]:   %.2f / %.2f',...
               100*settings.ranging.c./([1,settings.ranging.xcorr_osr] * settings.demodulation.frs));
            % interpolate , calculate correlation, and find maximum 
            %     interpolate signals (interp part I)
            rxi = interp(abs(input.rx(ind)), settings.ranging.xcorr_osr);
            fmi = interp(abs(input.fm(ind)), settings.ranging.xcorr_osr);
            yi = xcov(rxi, fmi, settings.ranging.xcorr_maxlag*settings.ranging.xcorr_osr, 'biased');
            %     interpolate xcorr result (interp part II) and find maximum
            [max_y, max_ind] = max(yi(settings.ranging.xcorr_maxlag*settings.ranging.xcorr_osr+1:end));
            % return modified settings and distance estimate
            output.settings = settings;
            output.dist_hat = max_ind / 2 * settings.ranging.c / (settings.ranging.xcorr_osr * settings.demodulation.frs);
            
         % interpolate signals and xcorr result (tradeoff between exactness and speed / RAM usage)
         %   ... interpolate signals by sqrt(oversampling rate)
         %   ... interpolate xcorr result by sqrt(oversampling rate) 
         %   => increased resolution in two steps by oversampling rate
         case 'sig+xcorr'
            % calculate oversampling rate and output resulting resolution(s)
            settings.ranging.xcorr_osr = ceil( sqrt(settings.ranging.c / (settings.ranging.res * settings.demodulation.frs)) );  
            msg('FMCW xcorr resolution original / interp / xcorr-interp [cm]:   %.2f / %.2f / %.2f',...
               100*settings.ranging.c./([1,settings.ranging.xcorr_osr,settings.ranging.xcorr_osr^2] *...
               settings.demodulation.frs));
            % interpolate , calculate correlation, interpolate, and find maximum 
            %     interpolate signals (interp part I)
            rxi = interp(abs(input.rx(ind)), settings.ranging.xcorr_osr);
            fmi = interp(abs(input.fm(ind)), settings.ranging.xcorr_osr);
            yi = xcov(rxi, fmi, settings.ranging.xcorr_maxlag*settings.ranging.xcorr_osr, 'biased');
            %     interpolate xcorr result (interp part II) and find maximum
            [max_y, max_ind] = max(interp(yi(settings.ranging.xcorr_maxlag*settings.ranging.xcorr_osr+1:end),...
               settings.ranging.xcorr_osr));
            % return modified settings and distance estimate
            output.settings = settings;
            output.dist_hat = max_ind / 2 * settings.ranging.c / (settings.ranging.xcorr_osr^2 * settings.demodulation.frs);
            
         % interpolate xcorr result only (fastest)
         %   ... interpolate xcorr result by oversampling rate
         case 'xcorr'
            % calculate oversampling rate and output resulting resolution(s)
            settings.ranging.xcorr_osr = ceil( settings.ranging.c / (settings.ranging.res * settings.demodulation.frs) );
            msg('FMCW xcorr resolution original / xcorr-interp [cm]:   %.2f / %.2f',...
               100*settings.ranging.c./([1,settings.ranging.xcorr_osr] * settings.demodulation.frs));
            
            % calculate correlation, interpolate for better resolution, and find maximum
            y = xcov(abs(input.rx(ind)), abs(input.fm(ind)), settings.ranging.xcorr_maxlag, 'none');
            [max_y, max_ind] = max(interp(y(settings.ranging.xcorr_maxlag+1:end), settings.ranging.xcorr_osr));
            % return modified settings and distance estimate
            output.settings = settings;
            output.dist_hat = max_ind / 2 * settings.ranging.c / (settings.ranging.xcorr_osr * settings.demodulation.frs);
            
         otherwise
            err('Unsupported mode for xcorr interpolation: "%s"', lower(settings.ranging.xcfmode));
      end
   

   % ****************************************************************************************************  
   otherwise
      err('Unsupported command: "%s".', lower(cmd));
end


return
    
% *******************************************************************************************************
% *******************************************************************************************************
% DEBUGGING SNIPPETS

% settings.ranging.mfcw_compsel.iirord = 2;
% settings.ranging.mfcw_compsel.att = 80;
% [c_mi, c_i, c_im] = mfcw_compsel(input, settings.ranging.mfcw_compsel);
% 
% close all
% 
% figure; hold on;
% plot(real(input), 'b')
% plot(imag(input), 'r')
% hold off; grid on; title('input');
% 
% figure; set(gca, 'yscale', 'log'); hold on;
% plot(abs(c_mi));
% plot(ind, abs(c_mi(ind,:)), 'r');
% hold off; grid on; title('c_{mi}'); %ylim([1e-7, 2e-5]);
% 
% figure; set(gca, 'yscale', 'log'); hold on;
% plot(abs(c_im));
% plot(ind, abs(c_im(ind,:)), 'r');
% hold off; grid on; title('c_{im}'); %ylim([1e-7, 2e-5]);


% *******************************************************************************************************
% OLD AND RUSTY
