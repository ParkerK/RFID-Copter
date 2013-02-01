% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% batch processing main script for: tag_modulation (reflection coefficient)
%
% WARNING: Settings are not checked for sanity!
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
% ***** Suitable for Batch-Processing *****
% This module is intended to be used as script (started manually and without the possibility of changing
% settings and/or the behavior). Nontheless it offers a minimalistic interface which enables automated 
% external batch processing. Please refer to the m-file itself for interface definition and behavior,
% as neither interface, nor behavior are fixed.
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

function tagchar_modulator_main(cat, rat, enr, fsr)

% *******************************************************************************************************
% initialization
if nargin == 0
   clear; close all; clc; pause(0.01);
end

% paths, etc
%     path to globalinit.m
path = fullfile(fileparts(mfilename('fullpath')), '..');
%     Matlab does not resolve symbolic links; if we're running on Unix, let the system resolve them
if isunix
   [dummy, path] = system(sprintf('readlink -f "%s"', path));
end
%     add this path and call globalinit script
restoredefaultpath; addpath(path); clear('dummy', 'path');
globalinit('sim');
% "switch to here"
cd(fileparts(mfilename('fullpath')));

if nargin == 0
   % one last chance to reconsider (only in manual mode)
   if nargin == 0
      disp(sprintf('WARNING: This may take a few hours and might overwrite existing characteristics!\n'));
      reply = input('Are you REALLY sure y/n [n]?  ', 's');
      if ~strcmpi(reply, 'y')
         return
      end
      disp(' ');
   end
end

% prevent full execution by version_system
if get_stacklevel > 0
   return
end


% *******************************************************************************************************
% standard settings (MIGHT BE OVERWRITTEN BELOW)

% TAG CONFIGURATIONS:
%   charid   zic   zmod->repl      za    picmin
%   _b-act    v5   act(v6)         v3      *
%   _b-psv    v5   psv(v7)->psv    v3      *       ... kept passive, max(Picmin) as threshold for repl.
%   _b-p2a    v5   psv(v7)->act    v3      *       ... made  active, max(Picmin) as threshold for repl.
%   _c-act    v5   act(v6)         v6      *
%   _c-psv    v5   psv(v7)->psv    v6      *       ... kept passive, max(Picmin) as threshold for repl.
%   _c-p2a    v5   psv(v7)->act    v6      *       ... made  active, max(Picmin) as threshold for repl.
%   _d-p2a    v5   psv(v7)->act    v6a     *       ... made  active, max(Picmin) as threshold for repl.

% general settings
%     basics
if nargin == 0
   settings.simchar = false; % suitability for simulator (true) or for theory only (false)
else
   settings.simchar = true; % automatic calls only to create a simulator char. (KEEP TRUE)
end
settings.writeres =  true; % also write result (characteristic) file or just workspace?
settings.charid   =  '_c-p2a'; %  ID that will be added to filename, start with '_' if not empty
settings.suffix   =    ''; % filename suffix (can be empty); generated automatically below
settings.fs       =   6e9; % Hz sampling frequency
%     system
if settings.simchar
   settings.renew_kticket =  true; % renew Kerberos ticket?
else
   settings.renew_kticket = false; % renew Kerberos ticket?
end
settings.emails        =  true; % send email if task is complete?
settings.maxfilesize   =   1.5; % GB maximum size for one file (might have problems loading otherwise)
%     folders (subfolders)
settings.folders.log =       'logs'; % log (diary)
settings.folders.wsp = 'workspaces'; % workspaceg (for plots, etc)
settings.folders.res =     'tagmod'; % results .mat-file ("tag modulator characteristic")

% data points per dimension
%     base vectors (esp. warping errors will increase dramatically for small values)
if settings.simchar
   settings.nf  = 160; % (old: 192, 384)
   settings.nm  =   6; % (old:  10,  12)
   settings.np  = 128; % (old: 128, 192)
else
   settings.nf  = 192;
   settings.nm  =   2; % light char: only linear model will be used => keep .nm=2
   settings.np  = 192; % (old: 256)
end
%     for assembly influence plot (settings.nat x settings.nat)
settings.nat = 250;

% tag characteristics
%     filenames
settings.charfiles.meas_tag_zic    = 'meas_tag_zic';      % chip input impedance
settings.charfiles.meas_tag_zmod   = 'meas_tag_zmod-psv'; % modulation impedance (modulated state)
settings.charfiles.meas_tag_za     = 'meas_tag_za-v6';    % antenna impedance
settings.charfiles.meas_tag_picmin = 'meas_tag_picmin';   % minimum operational chip input power
settings.charfiles.ads_mapping     = 'tagchar_modulator_adsm'; % assembly and detuning state mapping
%     modulation impedance
settings.zmod_max = complex(4e3, -15e3); % "unmodulated" state; modulated is taken from measurements
settings.rmod_min = 2;% Ohm minimum real part of modulation impedance ("prevent active Zmod")
%     special treatment of "passive" zmod characteristics (modulation not possible below min Pic)
%     ... Zic,mod approx Zic,unmod in this case, thus Zmod measurement close to and below min Pic subject
%         to heavy distortions => modulation impedance below picopmin replaced by zmod_max
%     ... threshold='fcn-of-f' will introduce distortions (non-differentiable reflection coefficient) 
%         due to limited resolution in combination with interpolation
%     ... Use picop_tol with care; the char. has the highest resolution close to picmin, thus a high
%         picop_tol shifts the "step" of zmod out of the highres area
settings.psv_zmod.mode      = 'act'; % {'psv': keep passive, 'act': make active}
settings.psv_zmod.threshold = 'max'; % {'max': use max(Picmin), 'fcn-of-f': use Picmin(f)}
settings.psv_zmod.picop_tol =     0; % dB tolerance area above picopmin (set to -Inf not to replace anything)
settings.psv_zmod.flen      =     1; % samples length of smoothing filter to avoid discontinuities

% power characteristics
%     limits for chip input power (extrapolation if necessary)
%     ... choose limits i.o.t. get full characteristic for mod and unmod; esp. the lower limit is
%         important here, as the power supply buffer is usually large enough to keep the tag powered at
%         picmin (unmod) during one reply => modulation below picmin(mod) is possible
%         resolution below picmin is minimal => better slightly too low than too high
settings.picmax   = 1e-1; % W
settings.picmin   = 2e-8; % W
%     limits for available power
settings.pavmax   =   10; % W (set high enough to get full char. even for mod)
settings.pavmin   =    0; % W (zero for automatic)
%     assumed "average" path-loss factor (linear indexing for log-dist model with this PLF)
if settings.simchar
   settings.plf = 3;
else
   settings.plf = 6; % light char: bias towards log scale (for ranging_detuning.m)
end
%     settings for iterations to determine optimal pav bounds
settings.npr_inop    = 0.98; % ratio of settings.np withing max(pav) through settings.pavopmin
if settings.simchar
   settings.pavit_maxit = 1000; % maximum number of iterations
   settings.pavit_bias  = 1e-3; % bias for next pav bounds 0<bias<=1 (smaller bias: more exact but slower)
else
   settings.pavit_maxit =  500; % light char: make calculation faster
   settings.pavit_bias  = 5e-2; % light char: make calculation faster
end
settings.npma        = round(settings.np/20); % length of MA filter (smoothing of power vector resolution)
   
% frequency characteristics
%     maximum/minimum frequency
settings.fmin = 1e3; % Hz
settings.fmax = settings.fs/2; % Hz
%     maximum/minimum operation frequency (regulations) => better resolution here
%     (will be adapted to match available frequency grid later)
settings.foptol = 40e6; % better resolution for operation area plus/minus settings.foptol Hz
settings.fopmax = 960e6 + settings.foptol; % Hz
settings.fopmin = 860e6 - settings.foptol; % Hz
%     settings.fopfactor times more data points within area of operation
settings.fopfactor = 2.5;
%     length of MA filter (smoothing of frequency vector resolution)
settings.nfma = round(settings.nf/5);
%     frequency resolution outside characteristic
settings.fres1 = 15; % bins/GHz for zone1 (fmin ... fchmin)
settings.fres4 =  3; % bins/GHz for zone4 (fmax ... fchmax)
%     polynomial models for extrapolation of frequency axis
settings.npoly_zic  = 2; % (unmodulated) chip impedance
settings.npoly_zmod = 2; % modulation impedance

% assembly tolerance (parallel RC)
settings.cat = 450e-15; % F
%     Q-matching (Rat-Calculation)
settings.fatm          = 915e6; % Hz frequency
settings.rat_x0        =   1e3; % Ohm initial value for numerical search for Q-matched Rat
settings.rat_min       =   1e0; % Ohm minimum parallel resistance
settings.rat_max       =   1e9; % Ohm maximum parallel resistance
settings.maxerr_qmatch =     1; % percent maximum error in matched Q (check optimization)
settings.rat_shift     =     0; % percent tolerance (shift) for Rat
%     min/max for assembly influence plot
settings.at_plot.rat_min =      1e1; % Ohm
settings.at_plot.rat_max =      1e4; % Ohm
settings.at_plot.cat_min =    1e-15; % F
settings.at_plot.cat_max = 1200e-15; % F

% detuning (antenna impedance)
%     detuning state
settings.detuning.res_en = 0; % factor for resonance boost -1..1 (-1: near water, 0: none, 1: near-metal)
settings.detuning.fshift = 0; % Hz shift in frequency (towards smaller frequencies)
%     constants for resonance boost (TEST CHANGES CAREFULLY!)
settings.detuning.weight_lin_y =   10; % nonlinearity: weight of linear part, R/X-warping
settings.detuning.weight_lin_f = 1.25; % nonlinearity: weight of linear part, f-warping
settings.detuning.weight_exp_y =    2; % nonlinearity: weight of exponential part, R/X-warping
settings.detuning.weight_exp_f =    4; % nonlinearity: weight of exponential part, f-warping
settings.detuning.weight_exp_w =    3; % nonlinearity: weight for combining trend and resonance
settings.detuning.fshift_attf  = 5e-9; % attenuation for resonance per Hz frequency shift
settings.detuning.nfilt        =    3; % samples filter length of MA filter (keep small!)

% maximum reflection coefficient (rho = 1 is not possible)
settings.rhomax     = 0.999;
settings.rhomax_tol =     1; % percent tolerance for IIR fitting

% random inverse checks rho(f, zmod, pav) -> Zic
if settings.simchar
   settings.checks = 1e4; % number of checks
else
   settings.checks = 1e2; % light char: make calculation faster
end
settings.maxerr = 0.02; % maximum relative mismatch for real/imaginary part
settings.maxavg = settings.maxerr/10; % maximum relative mismatch mean value (real/imaginary part)
%     demilitarized zone around NaN (border region) => no checks here (outliers)
settings.dmz_f = 6; % fch
settings.dmz_p = 6; % pav
settings.dmz_m = 1; % zmod

% IIR bandstop filter parameters to try
%     pronounced V-shape valley
settings.iir{1}.iirorder =      2; % filter order
settings.iir{1}.bshift   =    0.2; % 0...1/2; shift of IIR borders from char. border to center (0: none, 1/2: complete)
settings.iir{1}.cshift   =    0.6; % 0...1;   shift of IIR center towards minimum (0: none, 1: complete)
settings.iir{1}.dfmin    =   50e6; % Hz minimum initial bandwidth of IIR bandstop
%     V-U
settings.iir{2}.iirorder =      3;
settings.iir{2}.bshift   =   0.15;
settings.iir{2}.cshift   =    0.5;
settings.iir{2}.dfmin    =  100e6;
%     U-shape valley
settings.iir{3}.iirorder =      3;
settings.iir{3}.bshift   =   0.25;
settings.iir{3}.cshift   =    0.6;
settings.iir{3}.dfmin    =  150e6;
%     pronounced U-shape valley
settings.iir{4}.iirorder =      3;
settings.iir{4}.bshift   =    0.4;
settings.iir{4}.cshift   =    0.7;
settings.iir{4}.dfmin    =  150e6;

% IIR bandstop optimization (only magnitude is optimized)
%     minimum number of characteristic points to fit IIR
settings.minpts = 2; % < 2 is not allowed (slope detection)
%     options for lsqnonlin
settings.optimset = optimset('Display','off', 'LargeScale','on',...
   'TolX',1e-9, 'TolFun',5e-9, 'MaxIter',1e3, 'MaxFunEval',1e5);
%     initial point x0 [zeros, poles, gain]
settings.optimstart  = [0.98; 0.98; settings.rhomax];
%     bounds for x (upper and lower bound) [zeros, poles, gain]
settings.optimset.ub = [1-100*eps; 1-100*eps; 100]; % stable and minimum-phase
settings.optimset.lb = [0; 0; 0];
%     maximum number of optimizer retries to get |rho|<settings.rhomax
settings.maxoptretries = 2;
%     cost function weights (see function optim_costfcn for explanations)
settings.cf_w = [1, 3, 2, 0.1, 2, 3, 3, 50, 1]';
%     multiplier for weight alteration (will increase IIR residual error but improve fitting quality)
%     (e.g. in case of negative slope => more weight on |rho(f->)|=rhomax and less on slope)
settings.cf_walt = 3.5; % has to be > 1
%     multiplier for |rho|>rhomax punishment (multiplied by this factor for each retry)
settings.cf_wmult = 3; % has to be > 1 ! % 10
%     median filter to create smooth borders after IIR extrapolation
settings.bordermed_len = 5; % samples

% IIR bandstop residual error warning if maximum error above
settings.reserr_warn     = -30; % dB (bad fitting)
settings.reserr_critwarn = -10; % dB (very bad fitting)


% *******************************************************************************************************
% misc

% create directories (return values to suppress warning if directory already exists)
[dummy.stat,dummy.msg,dummy.msgid] = mkdir(settings.folders.log);
[dummy.stat,dummy.msg,dummy.msgid] = mkdir(settings.folders.wsp);
[dummy.stat,dummy.msg,dummy.msgid] = mkdir(settings.folders.res);

% configure email system
global globalsettings;
if settings.emails
   globalsettings.core.mailto      = 'daniel.arnitz@tugraz.at';
   globalsettings.core.subj_prefix = '[CHAR] '; % prefix for subject to identify email in clients
else
   globalsettings.core.mailto      = '';
   globalsettings.core.subj_prefix = '';
end

% load assembly and detuning state mappings
ads_mapping = load(settings.charfiles.ads_mapping);


% *******************************************************************************************************
% define sweep vectors
  
% sweep type (only if called as script and settings.simchar=false)
% sweepsettings.type = 'assembly';
% sweepsettings.type = 'detuning';
sweepsettings.type = 'all';
% sweepsettings.type = 'nxpanechoic_201005';
% sweepsettings.type = 'test';

% set up sweep (make sure these settings are part of the predefined detuning/assembly states)
if nargin == 0
   % simulator characteristics: a little of everything
   if settings.simchar
      cat = [450, 1250]; % fF
      rat = [-67, 0, 200]; % percent shift
      enr = [0, 0.5]; % boost of resonance (slightly weakened res. up to near-metal)
      fsr = [0, 100e6]; % frequency shift of resonance
   else
      switch(lower(sweepsettings.type))
         case 'assembly' % assembly detuning (up to massive impedances!)
            cat = [50:25:800, 850:50:1000, 1250:250:1750, 2000:1000:5000, 7500, 10000]; % fF
            rat = [-98, -95, -90, -80, -67, -50, -33, 0, 50, 100, 200, 500, 1000, 2000, 5000, 10000]; % percent shift
            enr = 0; % boost of resonance (perfectly tuned)
            fsr = 0; % frequency shift of resonance (perfectly tuned)
            sweepsettings.charid   = 'a';
         case 'detuning' % antenna detuning (up to massive detuning); assembly is good flip-chip ("perfect match")
            cat =  300; % fF
            rat =    0; % percent shift
            enr = -0.3 : 0.1 : 1.0; % boost of resonance (slightly weakened res. up to near-metal)
            fsr = 0 : 5e6 : 200e6; % frequency shift of resonance
            sweepsettings.charid   = 'd';
         case 'all' % assembly and antenna detuning (flip-chip assembly tolerance, moderate detuning)
            cat = [200:100:600]; % fF
            rat = [-50, -33, 0, 50, 100]; % percent shift
            enr = 0 : 0.1 : 0.5; % boost of resonance (none up to almost near-metal)
            fsr =  0 : 10e6 : 200e6; % frequency shift of resonance
            sweepsettings.charid   = 'ad';
         case 'nxpanechoic_201005' % SPECIAL: matching to measurements 
            cat = [450:25:600]; % fF
            rat = [-90:5:-40]; % percent shift
            enr = -0.2 : 0.1 : 0.2; % boost of resonance
            fsr =  10e6 : 2e6 : 30e6; % frequency shift of resonance
            sweepsettings.charid   = 'nxpanechoic_201005';
            % own ads_mapping (higher resolution)
            clear ads_mapping;
            ads_mapping.note = 'special mapping; for this file only';
            ads_mapping.assembly.description = 'Capacity .cat [F], shift of Ohmic part .rats [%] from optimum. Assembly state is index in .v_*.';
            ads_mapping.detuning.description = 'Antenna own-resonance boost .res_en -1..1 (0:none), frequency shift .fshift [Hz]. Detuning state is index in .v_*.';
            ads_mapping.assembly.cat  = cat*1e-15;
            ads_mapping.assembly.rats = rat;
            ads_mapping.assembly.v_cat  = sort(repmat(cat*1e-15, 1, length(rat)));
            ads_mapping.assembly.v_rats =      repmat(rat      , 1, length(cat));
            ads_mapping.detuning.res_en = enr;
            ads_mapping.detuning.fshift = fsr;
            ads_mapping.detuning.v_res_en = sort(repmat(enr, 1, length(fsr)));
            ads_mapping.detuning.v_fshift =      repmat(fsr, 1, length(enr));
            % reduce resolution of characteristic (meas. resolution is a lot smaller)
            settings.nf  =  96;
            settings.nm  =   2;
            settings.np  =  96;
            % make faster
            settings.checks = 1e1;
         case 'test'
            cat = 450; % fF
            rat =   000; % percent shift
            enr =   0; % boost of resonance (slightly weakened res. up to near-metal)
            fsr =   0; % frequency shift of resonance
            settings.charid   = '_test';
            sweepsettings.charid = '';
         otherwise
            error('Unsupported sweepsettings.type.');
      end
   end
end          

% map settings to a detuning status (quick and dirty)
sweepsettings.status_a = [];
for i = 1 : length(cat)
   for j = 1 : length(rat)
      sweepsettings.status_a = [sweepsettings.status_a,...
         find( round(ads_mapping.assembly.v_cat*1e15) == cat(i) & abs(ads_mapping.assembly.v_rats - rat(j)) < 1e-6)];
   end
end
sweepsettings.status_d = [];
for i = 1 : length(enr)
   for j = 1 : length(fsr)
      sweepsettings.status_d = [sweepsettings.status_d,...
         find( abs(ads_mapping.detuning.v_res_en - enr(i)) < 1e-6 & abs(ads_mapping.detuning.v_fshift - fsr(j)) < 1e-6)];
   end
end

% combine assembly and detuning sweeps
sweepsettings.a_full = sort(repmat(sweepsettings.status_a, 1, length(sweepsettings.status_d)));
sweepsettings.d_full =      repmat(sweepsettings.status_d, 1, length(sweepsettings.status_a));
%     check
if length(sweepsettings.a_full) ~= length(sweepsettings.d_full) || ...
      length(sweepsettings.a_full) ~= length(cat) * length(rat) * length(enr) * length(fsr)
   error('Mapping of settings to detuning/assembly states failed. Please check.');
end

% estimate file size (lightchar only)
est_filesize = length(sweepsettings.a_full) * settings.nf*settings.nm*settings.np * (16 + 8) / 1024^3; % complex double for rho_pav, double for m_pav
if ~settings.simchar  && est_filesize > settings.maxfilesize
   disp(sprintf('\n\nWARNING: The resulting characteristic will use roughly %.1f GB.\n', est_filesize));
   reply = input('Continue y/n [n]?  ', 's');
   if ~strcmpi(reply, 'y')
      return
   end
   disp(' ');
end


% *******************************************************************************************************
% call tagchar_modulator

% for all entries in sweep vectors
for i = 1 : length(sweepsettings.a_full)
   clc; disp(sprintf('%i of %i\n\n', i, length(sweepsettings.a_full))); %#ok<*DSPS>
   
   % set up states
   %     assembly
   settings.cat       = ads_mapping.assembly.v_cat ( sweepsettings.a_full(i) );
   settings.rat_shift = ads_mapping.assembly.v_rats( sweepsettings.a_full(i) );
   %     detuning
   settings.detuning.res_en = ads_mapping.detuning.v_res_en( sweepsettings.d_full(i) );
   settings.detuning.fshift = ads_mapping.detuning.v_fshift( sweepsettings.d_full(i) );
   
   % full characteristic (usable by simulator, written in individual files)
   if settings.simchar
      settings.suffix = sprintf('%s-sim_a%03i_d%03i', settings.charid, sweepsettings.a_full(i), sweepsettings.d_full(i));
      tagchar_modulator(settings);
      send_email(sprintf('Characteristic tagchar_modulator%s successfully created (%.0f of %.0f)',...
         settings.suffix, i, length(sweepsettings.a_full)), '');
   
   % light characteristic ()
   else
      settings.suffix = sprintf('%s_a%03i_d%03i', settings.charid, sweepsettings.a_full(i), sweepsettings.d_full(i));
      tagchar_mod{i}  = tagchar_modulator(settings);
      % extract and delete characteristic string from struct (no setting => identical for all) 
      characteristic = tagchar_mod{i}.characteristic;
      tagchar_mod{i} = rmfield(tagchar_mod{i}, 'characteristic');
   end
end

% write all light characteristics in one file
if ~settings.simchar
   clc; disp(sprintf('\nSaving final results...'));
   % create info header
   matfilename    = sprintf('tagchar_modulator%s_%s-light', settings.charid, sweepsettings.charid);
   characteristic = sprintf('%s\n%s',...
      'LIGHT TAG MODULATOR CHARACTERISTIC (NOT SUITABLE FOR SIMULATOR)', characteristic);
   createdon      = datestr(now, 0);
   createdby      = tagchar_mod{1}.createdby;
   usedmatfiles   = 'see individual entries of cell array for details';
   simchar        = false; % not suitable for simulator
   % save and clear workspace (save huge files in chunks; load might not work even if save -v7.3 does)
   fileinfo = whos('tagchar_mod');
   %     file too large => save in chunks
   if fileinfo.bytes/1024^3 > settings.maxfilesize % GB
      splitfile.size   = length(tagchar_mod);
      splitfile.chunks = ceil(fileinfo.bytes/1024^3/settings.maxfilesize);
      chunksize = round( splitfile.size / splitfile.chunks );
      fprintf('   File is too large (~%.1f GB). Saving in chunks: ', fileinfo.bytes/1024^3);
      for i = 1 : splitfile.chunks
         fprintf(' %i', i);
         splitfile.ind     = i;
         tagchar_mod_chunk = tagchar_mod(1+(i-1)*chunksize : min(length(tagchar_mod), i*chunksize));
         if i == 1
            matfilename_chunk = matfilename;
         else
            matfilename_chunk = [matfilename,'--',num2str(i,'%i')];
         end
         save('-v7.3', fullfile(settings.folders.res, matfilename_chunk), 'matfilename','characteristic','createdby','usedmatfiles',...
            'splitfile', 'tagchar_mod_chunk', 'ads_mapping', 'sweepsettings'); % -v7.3: large files
      end
      fprintf('\n');
   %     relatively small file => save in one file
   else
      save('-v7.3', fullfile(settings.folders.res, matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
         'tagchar_mod', 'ads_mapping', 'sweepsettings'); % -v7.3: large files
   end
   clear all;
end

% quit if this was done in batch mode
if nargin > 0
   quit;
end
