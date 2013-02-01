% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% simulator main function for MFCW gate simulation (multidimensional)
%
% 
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
% ***** Simulation Setup *****
% o) multi-reader, multi-tag
%    .) one active reader, all other readers are passive
%    .) each reader is activated once per loop pass (e.g., one tag position)
%    .) all tags reply to the query command (for interference analysis)
%    .) noiseless downlink channel (cf. signal model)
% o) tags are configured only once using a query command before ranging to speed up the process
% o) position of tag(s) is modified inside a loop (reader positions are fixed)
% o) distance reader(s) -> tag(s) is estimated using MFCW ranging during RN16 pilot tones
%    .) tags only send pilot plus header (0 bit RN16) for performance reasons
% o) the function creates a log, results and settings are written in an output mat-file
%
% o) SPECIAL MODE: tags can be used as field probes (power meters) only
%    .) tags not initialized, do not (de-)modulate either (not necessary, faster this way)
%    .) readers send a clean carrier, tag just records the available power level 
%    .) no ranging done
% 
% Warning: Tag RX path will be configured to carrier frequency of closest reader at initialization!
%
%
% ***** Interface definition *****
% function [results, settings] = sim_mfcw_ranging(settings)
%    settings   struct containing the specific simulator setup
%
%       .suffix          filename suffix for all output files (diary and results)
%       .renew_kticket   renew Kerberos ticket and AFS token?
%       .err_rep         try to repeat last loop pass in case of an error?
%       .numthreads      number of threads (will be saturated to number of CPUs available)
%       .hostname        name of host we're running on
%       .pid             process ID of this matlab session (main process)
%       .fs              sampling frequency in Hz
%       .c               speed of light in m/s
%
%       .specials        
%          .probemode        set all tags to "field probe mode" true/false, (record field strength only)
%          .probemode_clen   carrier length in "field probe mode" in s
%
%       .emails          struct containing settings for ETA emails
%          .eta_sendmail    sending of an eMail containing the ETA of a simulation
%          .eta_atloop      send ETA email at this loop index (use a few loops to get good ETA estimate) 
%          .eta_mintime     minimum overall simulation time in hrs to send ETA eMail
%
%       .folders   input/output folders (strings)
%          .root      root directory
%          .data      relative directory for data output (base: .ROOT)
%          .logs      relative directory for logging output (base: .ROOT)
%
%       .channel_global   struct containing basic global channel setups
%          .noise_on         activate noise?
%          .surf_on          reflective surfaces on?
%          .small_on         enable smallscale models?
%          .small_det        deterministic smallscale model (do not randomize impulse responses)?
%          .plf              average path-loss factor (2: free-space)
%          .n0               single-sided noise density [dBm, per ? Hz]
%          .v_dist           smallscale: distance vector in m for v_k and v_trms (!!! start at 0 m !!!)
%          .v_k              smallscale: Ricean K-factor vs. distance v_k(v_dist)
%          .v_trms           smallscale: RMS delay spread vs. distance v_trms(v_dist)
%          .pol_dim          dimension of polarization for all antennas (e.g. 2 for polarized along y in [x,y,z])
%          .surfaces         (optional) definitions of reflective surfaces, one substruct per surface
%             .ground           example reflective surface
%                .dim              dimension of normal vector (e.g. 3 for z in [x,y,z])
%                .shift            shift along normal vector (e.g. 1 means z=1 for dim=3)
%                .n2               refractive index of the reflective material (n1 assumed to be air: n1=1)
%
%       .loop    struct containing the loop setup
%          .n       number of iterations in loop
%          .pos_t   position of tags in m {[tag1 1...n]', [tag2 1...n]', ...} [x1,y1,z1; x2,y2...]
%
%       .readerpool   struct containing reader setups; cell arrays, one entry per reader
%          .n            number of readers (scalar)
%          .fc           carrier frequency in Hz
%          .ant          antenna type (characteristic filename); isotropic antenna if empty
%          .ant_rot      antenna rotation [azimuth, elevation] in degree ([0,0]: maximum in positive x-axis)
%          .ptx          EIRP transmit power in W
%          .t0           time delay in s
%          .id_ch        identical smallscale channels to this reader (monostatic antenna setup)? 
%          .virt            virtual transmitter? (0: no, 1: VTX created from TX1, ...)
%          .virt_gf         gain factor ("product of all reflection losses") of virtual transmitter
%          .virt_dim        dimension in which this VTX has been mirrored (e.g. 2 if mirror surface is xz)
%          .virt_surf       surface the VTX was last reflected in (VTX->RX must pass surface)
%          .virst_inv_surf  invert linked surface (false: VTX has to pass, true: VTX must not pass)
%          .virt_refl       number of reflections that lead to this VTX
%          .frs          reader sampling frequency in Hz
%          .quant        bits quantization
%          .pos          [x,y,z] position of readers in m
%          .mfcw         struct containing MFCW setup for each reader
%             .nc           number of secondary carriers
%             .fi           secondary carrier frequencies in Hz
%             .vari         secondary carrier variances (main carrier has variance 1)
%             .c_ord        order in which to combine the carriers (see mfcw_calcdist)
%
%       .tagpool     struct containing tag setups; cell arrays, one entry per tag
%          .n           number of tags (scalar)
%          .t0          time delay of reply in s 
%          .rn16        RN16 hex strings
%          .ant         antenna type (characteristic filename); isotropic antenna if empty
%          .ant_rot     antenna rotation [azimuth, elevation] in degree ([0,0]: maximum in positive x-axis)
%          .pwrchar     power supply characteristic {'Sofie_V2A', 'Sofie_V2B', 'Sofie_V2C', 'Sofie_AVG'}
%          .modcharid   ID for tag modulator characteristic (chip/antenna-type)
%          .adsm_file   filename of table that matches assembly (cat,rat) and detuning (enr,fsr) to a
%                       state (and thus to a tag modulator characteristic filename)
%          .cat         assembly capacity in F
%          .rat         shift for assembly resistance (tolerance) from optimal point in percent
%          .enr         boost of antenna resonance (-1..1; -1: no resonance, 0: no detuning, 1: near-metal)  
%          .fsr         Hz frequency shift of resonance (f.i.: fsr=100: shift from 940MHz to 840MHz)
%          .pos         [x,y,z] position of tags in m for initialization
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
% ! revise handling for non-probemode; do not demod/decode; do not count as receivers; channel R-R?, etc.
% ! do not initialize readerpool.n channel statistics; remove VTX from there
%
% *******************************************************************************************************


function [settings, results] = sim_mfcw_gate(settings)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   settings = version;
   return
end


% *******************************************************************************************************
% internal settings

% changes of fs to avoid low beat frequencies caused by sampling (for time-variant energy detection)
internalsettings.bnd_osr =   10; % if oversampling rate is smaller that bnd_ord ...
internalsettings.tol_osr = 1e-1; % ... it may not be closer to an integer than tol_osr


% *******************************************************************************************************
% logging

% start timer and diary
tic;
settings.diaryfilename = fullfile(settings.folders.root, settings.folders.logs, sprintf('%s.log', settings.suffix)); 
system(sprintf('rm -f %s', settings.diaryfilename));% (do not append diary)
diary(settings.diaryfilename);

% create logfile header 
disp('*******************************************************************************************************');
disp('*******************************************************************************************************');
disp('***');
disp(sprintf('*** Created by %s.m, version %s, %s on host "%s"',...
   mfilename, version, datestr(now, 'local'), settings.hostname));
disp('***');
disp('*******************************************************************************************************');
disp('*******************************************************************************************************');

% access globalsettings
global globalsettings;



% *******************************************************************************************************
% initialize simulator, load standard setup and create individuals
headline('\n\n');
headline('*******************************************************************************************************');
headline('* Initializing Simulator and Creating Individuals');

% multithreading
%     initialize Matlab pool for parallel loops
% matlabpool('local', min(settings.maxthreads, settings.numthreads));

% make sure oversampling rate ist not close to an integer
% (would result in very low beat frequency which cannot be handled by tag_demodulation)
osr = settings.fs./cell2mat(settings.readerpool.fc);
if any( and(osr < internalsettings.bnd_osr, mod(osr,1) < internalsettings.tol_osr) )
   critwarn('Low oversampling rate for carrier is close to an integer. Increasing sampling frequency.');
   % try to fix
   for i = 1 : ceil(1/internalsettings.tol_osr) - 1
      settings.fs = settings.fs + min(cell2mat(settings.readerpool.fc)) * internalsettings.tol_osr;
      osr = settings.fs./cell2mat(settings.readerpool.fc);
      if ~any( and(osr < internalsettings.bnd_osr, mod(osr,1) < internalsettings.tol_osr) )
         break;
      end
   end
   % check one last time
   if any( and(osr < internalsettings.bnd_osr, mod(osr,1) < internalsettings.tol_osr) )
      err('Unable to fix oversampling rate problem.');
   end
end

% output
disp(sprintf('\nPlease note:'));
%     renew tickets only if not in probemode (loop exec. time in probemode is in the range of seconds or below)
if settings.specials.probemode && settings.renew_kticket
   settings.renew_kticket = false;
   disp(sprintf('   .) Renewal of Kerberos/AFS tickets has been deactivated (Field Probe Mode).'));
end
%     settings content check has been switched off
if ~globalsettings.core.settings_check
   disp(sprintf('   .) Settings content check has been switched off for performance reasons.'));
end
%     sampling resolution
disp(sprintf('   .) Delays will be rounded to sampling interval, resolution is %g ns (%g cm).',...
   1e9/settings.fs, 100*settings.c/settings.fs));
%     logging "warnings"
if ~globalsettings.logging.warnings || ~globalsettings.logging.headlines || ~globalsettings.logging.messages
   disp(sprintf('\nThe following output is switched off:'));
    if ~globalsettings.logging.warnings
      disp('   .) Warning output. You will not receive any warnings!');
   end
   if ~globalsettings.logging.headlines
      disp('   .) Headlines and progress output (ETA, ...).');
   end
   if ~globalsettings.logging.messages
      disp('   .) Message output. You will not receive any detailed information.');
   end
end     
      
% load standard simulator setup
headline('\nLoading Standardsettings');
stdsettings = sim_stdsettings_mfcw(settings);

% load assembly/detuning state mapping (for all tags)
for i = 1 : settings.tagpool.n
   adsm{i} = load(settings.tagpool.adsm_file{i});
end

% global modification of standard settings
%     channel switches
stdsettings.channel.probemode      = settings.specials.probemode;
stdsettings.channel.smallscale.on  = settings.channel_global.small_on;
stdsettings.feedback.smallscale.on = settings.channel_global.small_on;
stdsettings.channel.smallscale.det  = settings.channel_global.small_det;
stdsettings.feedback.smallscale.det = settings.channel_global.small_det;
if ~settings.channel_global.noise_on; 
   stdsettings.channel.noise.type  = 'off'; % should be off anyway
   stdsettings.feedback.noise.type = 'off';
else
   stdsettings.channel.noise.type  = 'off'; % cf. signal model (should be off anyway)
   stdsettings.feedback.noise.type = 'wgn'; % cf. signal model
end
%     smallscale channel
stdsettings.channel.type   = settings.channel_global.type;
stdsettings.channel.v_dist = settings.channel_global.v_dist;
stdsettings.channel.v_k    = settings.channel_global.v_k;
stdsettings.channel.v_trms = settings.channel_global.v_trms;
%     noise
stdsettings.channel.noise.n0 = settings.channel_global.n0;
%     reflective surface(s)
if settings.channel_global.surf_on
   if ~isfield(settings.channel_global, 'surfaces')
      critwarn('Surfaces switched on in global channel settings, but none defined. Switching off.');
      settings.channel_global.surf_on = false;
   else
      stdsettings.channel.surfaces  = settings.channel_global.surfaces;
      stdsettings.feedback.surfaces = settings.channel_global.surfaces;
      %     copy polarization
      refl_names = fieldnames(stdsettings.channel.surfaces);
      for k = 1 : length(refl_names)
         stdsettings.channel.surfaces.(refl_names{k}).pol_dim  = settings.channel_global.pol_dim;
         stdsettings.feedback.surfaces.(refl_names{k}).pol_dim = settings.channel_global.pol_dim;
      end
      clear('refl_names');
   end
end

% readerpool
headline('\nReaderpool');
settings.readerpool.n_nonvirt = sum(cellfun(@(x) x==0, settings.readerpool.virt));
if find(cellfun(@(x) x==0, settings.readerpool.virt)==0, 1, 'first') - 1 ~= settings.readerpool.n_nonvirt
   err('Non-virtual transmitters have to be placed at the beginning of readerpool.');
end
% [!!!] BUG
if settings.readerpool.n_nonvirt < settings.readerpool.n && ~settings.specials.probemode
   critwarn('There is a bug in the handling of VTX in standard mode. Sorry.'); % partial fix implemented: reduced to warning
end
for i = 1 : settings.readerpool.n
   % copy from standardsettings
   settings.reader{i} = stdsettings.reader;
   % truncate fi and vari to nc
   settings.readerpool.mfcw.fi{i}   = settings.readerpool.mfcw.fi{i}(1:settings.readerpool.mfcw.nc{i});
   settings.readerpool.mfcw.vari{i} = settings.readerpool.mfcw.vari{i}(1:settings.readerpool.mfcw.nc{i});
   
   % setup modulation and complete/sanitize the settings
   settings.reader{i}.modulation.t0 = settings.readerpool.t0{i};
   settings.reader{i} = reader_main('initialize', settings.reader{i});
   
   % reader component setup
   settings.reader{i}.oscillator.fcenter = settings.readerpool.fc{i};
   settings.reader{i}.analogpathest.f0   = settings.readerpool.fc{i};
   settings.reader{i}.transmitter.ptx    = settings.readerpool.ptx{i};
   settings.reader{i}.demodulation.frs   = settings.readerpool.frs{i};
   settings.reader{i}.demodulation.fcut  =...
      3 * ( max(abs(settings.readerpool.mfcw.fi{i})) + settings.reader{i}.modulation.lf );
   settings.reader{i}.demodulation.q     = settings.readerpool.quant{i};
   
   % MFCW ranging setup (oscillator setup will be done later)
   %     general
   fm = settings.reader{i}.modulation.lf;
   fi = settings.readerpool.mfcw.fi{i};
   settings.readerpool.mfcw.fi_rs{i} = subsampling(fi, settings.reader{i}.demodulation.frs);
   if any(abs(fi - settings.readerpool.mfcw.fi_rs{i}) > eps)
      warn('Performing subsampling for reader %d.', i)
   end
   if any(abs(diff(abs(settings.readerpool.mfcw.fi_rs{i}))) <= 3*fm)
      err('Detected overlapping of carriers/sidebands. Check frequency setup (subsampling?).')
   end
   settings.reader{i}.ranging.freq = [-fm, 0, fm, fi-fm, fi, fi+fm]; % for analog input stage estimator
   %     secondary carrier generation
   settings.reader{i}.ranging.mfcw_addseccarriers.f0   = settings.readerpool.fc{i};
   settings.reader{i}.ranging.mfcw_addseccarriers.nc   = settings.readerpool.mfcw.nc{i};
   settings.reader{i}.ranging.mfcw_addseccarriers.fi   = settings.readerpool.mfcw.fi{i};
   settings.reader{i}.ranging.mfcw_addseccarriers.vari = settings.readerpool.mfcw.vari{i};
   %     component selection
   settings.reader{i}.ranging.mfcw_compsel.nc  = settings.readerpool.mfcw.nc{i};
   settings.reader{i}.ranging.mfcw_compsel.fi  = settings.readerpool.mfcw.fi_rs{i};
   settings.reader{i}.ranging.mfcw_compsel.lf  = fm;
   settings.reader{i}.ranging.mfcw_compsel.bw  = fm/2;
   settings.reader{i}.ranging.mfcw_compsel.frs = settings.reader{i}.demodulation.frs;
   %     distance estimator setup
   settings.reader{i}.ranging.mfcw_calcdist.nc    = settings.readerpool.mfcw.nc{i};
   settings.reader{i}.ranging.mfcw_calcdist.fi    = settings.readerpool.mfcw.fi{i};
   settings.reader{i}.ranging.mfcw_calcdist.c_ord = settings.readerpool.mfcw.c_ord{i};
   %     cleanup
   clear('fm', 'fi');
end

% tagpool
%     field probe mode: just parameters relevant for the channel
if settings.specials.probemode
   headline('\nTagpool (Probes)');
   critwarn('Setting probe input carrier frequency to frequency of closest reader.');
   for i = 1 : settings.tagpool.n
      % copy from standardsettings
      settings.tag{i} = stdsettings.tag;
      % set carrier frequency to fc of closest reader (warning above)
      % ... this should be considered TEMPORARY (until tag rx path is ultrawideband)
      [max_r, ind_r]  = min(get_distances(settings.readerpool.pos, settings.tagpool.pos{i}));
      settings.tag{i}.modulation.fc   = settings.readerpool.fc{ind_r};
      settings.tag{i}.demodulation.fc = settings.readerpool.fc{ind_r};
   end
   
%     normal operation
else
   headline('\nTagpool');
   critwarn('Settings tag input carrier frequency to frequency of closest reader.');
   for i = 1 : settings.tagpool.n
      % copy from standardsettings
      settings.tag{i} = stdsettings.tag;
      % individual setup
      settings.tag{i}.id.rn16             = settings.tagpool.rn16{i};
      [settings.tag{i}.adsm.as, settings.tag{i}.adsm.ds] = adstate('settings->state', adsm{i},...
         settings.tagpool.cat{i}, settings.tagpool.rat{i}, settings.tagpool.enr{i}, settings.tagpool.fsr{i});
      settings.tag{i}.modulation.charfile = ...
         sprintf('tagchar_modulator_%s-sim_a%03i_d%03i', settings.tagpool.modcharid{i}, ...
         settings.tag{i}.adsm.as, settings.tag{i}.adsm.ds); % filename: ID, assembly state, detuning state
      settings.tag{i}.modulation.t0       = settings.tagpool.t0{i};
      settings.tag{i}.modulation.pwrchar = settings.tagpool.pwrchar{i};
      
      % set carrier frequency to fc of closest reader (warning above)
      % ... this should be considered TEMPORARY (until tag rx path is ultrawideband)
      [max_r, ind_r] = min(get_distances(settings.readerpool.pos, settings.tagpool.pos{i}));
      settings.tag{i}.modulation.fc   = settings.readerpool.fc{ind_r};
      settings.tag{i}.demodulation.fc = settings.readerpool.fc{ind_r};
      
      % initialize tag (complete and sanitize settings)
      settings.tag{i} = tag_main('initialize', settings.tag{i});
   end
end

% channels reader <-> tag
headline('\nChannels Reader <-> Tag');
%     create channels
for i = 1 : settings.readerpool.n
   for j = 1 : settings.tagpool.n
         % copy from standardsettings
         %     reader -> tag (special treatment for virtual readers)
         settings.channel_rt{i,j}.probemode = stdsettings.channel.probemode;
         settings.channel_rt{i,j}.c  = stdsettings.channel.c;
         settings.channel_rt{i,j}.fs = stdsettings.channel.fs;
         if settings.readerpool.virt{i}
            settings.channel_rt{i,j}.type            = 'outdoor'; % forced long-range model for virtual readers
            settings.channel_rt{i,j}.ucg             = settings.readerpool.virt_gf{i}; % plus gain factor
            settings.channel_rt{i,j}.virtual.vtx     = true; % transmitter is virtual
            settings.channel_rt{i,j}.virtual.vrx     = false;
            settings.channel_rt{i,j}.virtual.otx     = settings.readerpool.virt{i};
            settings.channel_rt{i,j}.virtual.orx     = NaN;
            settings.channel_rt{i,j}.virtual.dim     = settings.readerpool.virt_dim{i};
            settings.channel_rt{i,j}.virtual.refl    = settings.readerpool.virt_refl{i};
            if ~isempty(settings.readerpool.virt_surf{i})
               settings.channel_rt{i,j}.virtual.surface = settings.readerpool.virt_surf{i};
               settings.channel_rt{i,j}.virtual.invsurf = settings.readerpool.virt_inv_surf{i};
            end
         else
            settings.channel_rt{i,j}.type     = stdsettings.channel.type;
            settings.channel_rt{i,j}.virtual.vtx  = false; % transmitter is real
            settings.channel_rt{i,j}.virtual.vrx  = false;
            settings.channel_rt{i,j}.virtual.otx  = NaN;
            settings.channel_rt{i,j}.virtual.orx  = NaN;
            settings.channel_rt{i,j}.virtual.dim  = NaN;
            settings.channel_rt{i,j}.virtual.refl = NaN;
         end
         settings.channel_rt{i,j}.largescale  = stdsettings.channel.largescale;
         if settings.channel_global.surf_on
            settings.channel_rt{i,j}.surfaces = stdsettings.channel.surfaces;
         end
         settings.channel_rt{i,j}.smallscale  = stdsettings.channel.smallscale;
         settings.channel_rt{i,j}.noise       = stdsettings.channel.noise;
         %     tag -> reader (special treatment for virtual readers)
         settings.channel_tr{j,i}.probemode = stdsettings.channel.probemode;
         settings.channel_tr{j,i}.c  = stdsettings.channel.c;
         settings.channel_tr{j,i}.fs = stdsettings.channel.fs;
         if settings.readerpool.virt{i}
            settings.channel_tr{j,i}.type            = 'outdoor'; % forced long-range model for virtual readers
            settings.channel_tr{j,i}.ucg             = settings.readerpool.virt_gf{i}; % plus gain factor
            settings.channel_tr{j,i}.virtual.vtx     = false;
            settings.channel_tr{j,i}.virtual.vrx     = true; % receiver is virtual
            settings.channel_tr{j,i}.virtual.otx     = NaN;
            settings.channel_tr{j,i}.virtual.orx     = settings.readerpool.virt{i};
            settings.channel_tr{j,i}.virtual.dim     = settings.readerpool.virt_dim{i};
            settings.channel_tr{j,i}.virtual.refl    = settings.readerpool.virt_refl{i};
            if ~isempty(settings.readerpool.virt_surf{i})
               settings.channel_tr{j,i}.virtual.surface = settings.readerpool.virt_surf{i};
               settings.channel_tr{j,i}.virtual.invsurf = settings.readerpool.virt_inv_surf{i};
            end
         else
            settings.channel_tr{j,i}.type     = stdsettings.channel.type;
            settings.channel_tr{j,i}.virtual.vtx  = false;
            settings.channel_tr{j,i}.virtual.vrx  = false; % receiver is real
            settings.channel_tr{j,i}.virtual.otx  = NaN;
            settings.channel_tr{j,i}.virtual.orx  = NaN;
            settings.channel_tr{j,i}.virtual.dim  = NaN;
            settings.channel_tr{j,i}.virtual.refl = NaN;
         end
         settings.channel_tr{j,i}.largescale  = stdsettings.channel.largescale;
         if settings.channel_global.surf_on
            settings.channel_tr{j,i}.surfaces = stdsettings.channel.surfaces;
         end
         settings.channel_tr{j,i}.smallscale  = stdsettings.channel.smallscale;
         settings.channel_tr{j,i}.noise       = stdsettings.channel.noise;
         % largescale
         settings.channel_rt{i,j}.largescale.f0 = settings.readerpool.fc{i};
         settings.channel_rt{i,j}.largescale.pl = settings.channel_global.plf;
         settings.channel_tr{j,i}.largescale.f0 = settings.readerpool.fc{i};
         settings.channel_tr{j,i}.largescale.pl = settings.channel_global.plf;
         % directivity
         %     reader -> tag
         settings.channel_rt{i,j}.directivity.txant = settings.readerpool.ant{i};
         settings.channel_rt{i,j}.directivity.txrot = settings.readerpool.ant_rot{i};
         settings.channel_rt{i,j}.directivity.rxant = settings.tagpool.ant{j};
         settings.channel_rt{i,j}.directivity.rxrot = settings.tagpool.ant_rot{j};
         %     tag -> reader
         settings.channel_tr{j,i}.directivity.txant = settings.tagpool.ant{j};
         settings.channel_tr{j,i}.directivity.txrot = settings.tagpool.ant_rot{j};
         settings.channel_tr{j,i}.directivity.rxant = settings.readerpool.ant{i};
         settings.channel_tr{j,i}.directivity.rxrot = settings.readerpool.ant_rot{i}; 
         % smallscale
         settings.channel_rt{i,j}.smallscale.bw   = 3 * ( 2*settings.reader{i}.modulation.lf +...
            max([0, settings.readerpool.mfcw.fi{i}]) - min([0, settings.readerpool.mfcw.fi{i}]) );
         settings.channel_tr{j,i}.smallscale.bw   = 3 * ( 2*settings.reader{i}.modulation.lf +...
            max([0, settings.readerpool.mfcw.fi{i}]) - min([0, settings.readerpool.mfcw.fi{i}]) );
         settings.channel_rt{i,j}.smallscale.fres = min(diff(sort([0, settings.readerpool.mfcw.fi{i}]))); % smallest carrier spacing
         settings.channel_tr{j,i}.smallscale.fres = min(diff(sort([0, settings.readerpool.mfcw.fi{i}])));
         % noise
         settings.channel_rt{i,j}.noise.frxs = settings.tag{j}.clock.fcenter; % receiver is tag
         settings.channel_tr{j,i}.noise.frxs = settings.reader{i}.demodulation.frs; % receiver is reader
   end
end
%     add all distance-dependent settings
settings.channel_rt = channel_newpos(settings.channel_rt, ...
   settings.readerpool.pos(1:settings.readerpool.n), settings.tagpool.pos(1:settings.tagpool.n), stdsettings.channel);
settings.channel_tr = channel_newpos(settings.channel_tr, ...
   settings.tagpool.pos(1:settings.tagpool.n), settings.readerpool.pos(1:settings.readerpool.n), stdsettings.channel);

% channels reader -> reader (not needed in probe-mode)
if ~settings.specials.probemode
   headline('\nChannels Reader -> Reader');
   %     create channels (
   for i = 1 : settings.readerpool.n
      for j = 1 : settings.readerpool.n_nonvirt
         % copy from standardsettings (special treatment for virtual readers)
         settings.channel_rr{i,j}.probemode = stdsettings.channel.probemode;
         settings.channel_rr{i,j}.c  = stdsettings.channel.c;
         settings.channel_rr{i,j}.fs = stdsettings.channel.fs;
         if settings.readerpool.virt{i}
            settings.channel_rr{i,j}.type = 'outdoor'; % forced long-range model for virtual readers
            settings.channel_rr{i,j}.ucg  = settings.readerpool.virt_gf{i}; % plus gain factor
            settings.channel_rr{i,j}.virtual.vtx  = true;
            settings.channel_rr{i,j}.virtual.otx  = settings.readerpool.virt{i};
            settings.channel_rr{i,j}.virtual.dim  = settings.readerpool.virt_dim{i};
            settings.channel_rr{i,j}.virtual.refl = settings.readerpool.virt_refl{i};
            if ~isempty(settings.readerpool.virt_surf{i})
               settings.channel_rr{i,j}.virtual.surface = settings.readerpool.virt_surf{i};
               settings.channel_rr{i,j}.virtual.invsurf = settings.readerpool.virt_inv_surf{i};
            end
         else
            settings.channel_rr{i,j}.type         = stdsettings.feedback.type;
            settings.channel_rr{i,j}.virtual.vtx  = false;
            settings.channel_rr{i,j}.virtual.otx  = NaN;
            settings.channel_rr{i,j}.virtual.dim  = NaN;
            settings.channel_rr{i,j}.virtual.refl = NaN;
         end
% %          if settings.readerpool.virt{j}
% %             settings.channel_rr{i,j}.virtual.vrx  = true;
% %             settings.channel_rr{i,j}.virtual.orx  = settings.readerpool.virt{j};
% %             settings.channel_rr{i,j}.virtual.dim  = settings.readerpool.virt_dim{j};
% %             settings.channel_rr{i,j}.virtual.refl = settings.readerpool.virt_refl{j};
% %          else
            settings.channel_rr{i,j}.virtual.vrx = false;
            settings.channel_rr{i,j}.virtual.orx = NaN;
% %          end  
         settings.channel_rr{i,j}.smallscale = stdsettings.feedback.smallscale;
         settings.channel_rr{i,j}.noise      = stdsettings.feedback.noise;
         if i == j % direct feedback to identical reader
            settings.channel_rr{i,j}.direct.gain     = 10^(-stdsettings.feedback.fb_att/20);
            settings.channel_rr{i,j}.direct.delay_s  = round(stdsettings.feedback.fb_del*settings.fs);
            settings.channel_rr{i,j}.direct.gain_dir = 1;    
         else % or normal channel to different reader
            % largescale
            settings.channel_rr{i,j}.largescale     = stdsettings.feedback.largescale;
            settings.channel_rr{i,j}.largescale.f0  = settings.readerpool.fc{i};
            settings.channel_rr{i,j}.largescale.pl  = settings.channel_global.plf;
            % reflective surfaces
            if settings.channel_global.surf_on
               settings.channel_rr{i,j}.surfaces = stdsettings.feedback.surfaces;
            end
            % directivity
            settings.channel_rr{i,j}.directivity       = stdsettings.feedback.directivity;
            settings.channel_rr{i,j}.directivity.txant = settings.readerpool.ant{i};
            settings.channel_rr{i,j}.directivity.txrot = settings.readerpool.ant_rot{i};
            settings.channel_rr{i,j}.directivity.rxant = settings.readerpool.ant{j};
            settings.channel_rr{i,j}.directivity.rxrot = settings.readerpool.ant_rot{j};
         end
         % smallscale
         settings.channel_rr{i,j}.smallscale.bw   = ...
            3 * ( max([0, settings.readerpool.mfcw.fi{i}]) - min([0, settings.readerpool.mfcw.fi{i}]) );
         settings.channel_rr{i,j}.smallscale.fres = min(diff(sort([0, settings.readerpool.mfcw.fi{i}])));
         % noise
         settings.channel_rr{i,j}.noise.frxs = settings.reader{j}.demodulation.frs;
      end
   end
   %     add all distance-dependent settings
   settings.channel_rr = channel_newpos(settings.channel_rr, ...
      settings.readerpool.pos(1:settings.readerpool.n), settings.readerpool.pos(1:settings.readerpool.n_nonvirt), stdsettings.feedback);
end

% variable cleanup
clear('max_r', 'ind_r');


% *******************************************************************************************************
% run simulator: downlink (query) ... initialization of tags (not necessary in probe mode)
headline('\n');

% field probe mode: nothing to initialize, just print a message
if settings.specials.probemode
   headline('Simulator is in FIELD PROBE MODE');
   
% normal operation: initialization of tags (query)
else
   headline('*******************************************************************************************************');
   headline('* Downlink (Query command)');
   
   % temporarily switch off multipath propagation
   for i = 1 : settings.readerpool.n
      for j = 1 : settings.tagpool.n
         tempsettings.small_on{i,j} = settings.channel_rt{i,j}.smallscale.on;
         settings.channel_rt{i,j}.smallscale.on = false;
      end
   end
   
   % prepare command
   headline('\nPreparing Query Command');
   for i = 1 : settings.readerpool.n_nonvirt
      headline('   Reader %i', i);
      settings.reader{i} = reader_main('prep_query', settings.reader{i});
   end
   % find maximum needed carrier length
   %     data only
   max_clen_s = ones(settings.tagpool.n,1) *...
      ceil(max(cellfun(@(x) x.modulation.length, settings.reader)) * settings.fs); % [samples]
   %     add the maximum overhead required by the tag implementation (group delays, ...)
   for i = 1 : settings.tagpool.n
      max_clen_s(i) = tag_main('clen_addoverhead', max_clen_s(i), settings.tag{i});
   end
   %     and distribute to all readers (add field max_clen to settings)
   settings.reader = cellfun(@(x) setfield(x, 'max_clen', max(max_clen_s)/settings.fs),...
      settings.reader, 'UniformOutput',false); %#ok<SFLD>
   
   % modulate command
   headline('\nModulate');
   for i = 1 : settings.readerpool.n_nonvirt
      if i == 1
         headline('   Reader %i', i);
         % set carrier length to maximum needed length
         settings.reader{i} = reader_main('set_maxclen', settings.reader{i});
         % create modulated carrier
         reader_carrier{i} = reader_main('tx_data', settings.reader{i});
      else
         % create dummy
         reader_carrier{i} = [0];
      end
   end
   
   % channel R -> T
   headline('\nChannel Reader -> Tag');
   tag_carrier = channel_main(reader_carrier, settings.channel_rt);
   
   % demodulate, decode, determine state
   headline('\nDemodulate + Decode, Get Link Setup');
   for i = 1 : settings.tagpool.n
      headline('   Tag %i', i);
      % re-initialize tag (set to unpowered state)
      settings.tag{i} = tag_main('re-initialize', settings.tag{i});
      % receive, check if powered
      temp = tag_main('rx', tag_carrier{i}, settings.tag{i});
      settings.tag{i} = temp.settings; % modifications to vdda
      tag_rxsignal{i} = temp.rxsignal;
      settings.tag{i}.state.powered = temp.powered;
      % set EPC state accordingly ('' or 'ready')
      settings.tag{i} = tag_main('set_epcstate', settings.tag{i});
      % if not powered => not much sense in decoding
      if ~settings.tag{i}.state.powered; continue; end
      % demodulate and decode
      temp = tag_main('rx_data', tag_rxsignal{i}, settings.tag{i});
      %    tag_decoded{i} = temp.decoded; NOT NEEDED
      settings.tag{i}.state.linkinfo = temp.linkinfo;
      % determine state of tag
      settings.tag{i} = tag_main('set_epcstate', settings.tag{i});
      % setup return link (if possible)
      settings.tag{i} = tag_main('query', settings.tag{i});
   end
   
   % return smallscale channel to its original state
   for i = 1 : settings.readerpool.n
      for j = 1 : settings.tagpool.n
         settings.channel_rt{i,j}.smallscale.on = tempsettings.small_on{i,j};
      end
   end
   
   % free memory (large vectors) and clean up variables
   clear('reader_carrier', 'tag_carrier', 'tag_rxsignal', 'temp', 'nfft', 'max_clen_s', 'tempsettings');
   
   % check for inactive tags before continuing (all tags need to be initialized here)
   if any( cellfun(@(x) ~x.state.powered, settings.tag) )
      temp = 1 : settings.tagpool.n;
      err('The following tags were not properly initialized: %s(out of tags 1..%i)',...
         sprintf('%i ', temp(cellfun(@(x) ~x.state.powered, settings.tag))), settings.tagpool.n );
      diary off;
   end
end


% *******************************************************************************************************
% run simulator: prepare uplink
headline('\n');

% field probe mode: carrier only
if settings.specials.probemode
   headline('*******************************************************************************************************');
   headline('* Prepare Uplink (field probe mode)');
   % distribute carrier length setting to all readers ("maximum carrier length")
   settings.reader = cellfun(@(x) setfield(x, 'max_clen', settings.specials.probemode_clen), settings.reader, 'UniformOutput',false); %#ok<SFLD>
   % create carrier
   for i = 1 : settings.readerpool.n_nonvirt
      % set length to maximum carrier length
      settings.reader{i} = reader_main('set_maxclen', settings.reader{i});
      % create unmodulated carrier
      reader_carrier_base{i} = reader_main('tx_carrier', settings.reader{i});
   end
   
% normal operation: RN16
else
   headline('*******************************************************************************************************');
   headline('* Prepare Uplink (RN16 reply of active tags)');
   % encode and find needed carrier length
   headline('\nEncoding RN16 and Initializing Modulation');
   for i = 1 : settings.tagpool.n
      headline('   Tag %i', i);
      settings.tag{i} = tag_main('prep_rn16', settings.tag{i});
   end
   % find maximum needed carrier length
   max_clen = max(cellfun(@(x) x.modulation.length, settings.tag));
   %     and distribute to all readers (add field max_clen to settings)
   settings.reader = cellfun(@(x) setfield(x, 'max_clen', max_clen), settings.reader, 'UniformOutput',false); %#ok<SFLD>
   
   % create carrier for modulation (including MFCW ranging signals)
   headline('\nCreating Carrier for Modulation (Including MFCW Secondary Carriers)');
   for i = 1 : settings.readerpool.n_nonvirt
      headline('   Reader %i', i);
      % setup MFCW system
      settings.reader{i} = reader_main('prep_mfcw', settings.reader{i});
      % set length to maximum carrier length
      settings.reader{i} = reader_main('set_maxclen', settings.reader{i});
      % create unmodulated MFCW carrier
      reader_carrier_base{i} = reader_main('tx_mfcw', settings.reader{i});
   end
end

% *******************************************************************************************************
% loop position and perform MFCW raging 

% allocate results
%     expected distance
results.expdist = cell(settings.loop.n,1);
results.expdist = cellfun(@(x) nan(settings.readerpool.n_nonvirt, settings.readerpool.n_nonvirt, settings.tagpool.n), results.expdist, 'UniformOutput', false);
%     channel parameters Reader -> Tag (one sender)
results.channel_rt = cell(settings.readerpool.n_nonvirt, settings.tagpool.n);
results.channel_rt = cellfun(@(x) struct(...
   'pdp_los',nan(settings.loop.n, 1), 'pdp_sum',nan(settings.loop.n, 1),...
   'k',nan(settings.loop.n, 1), 'trms',nan(settings.loop.n, 1),...
   'tmin',nan(settings.loop.n, 1), 'tmax',nan(settings.loop.n, 1),...
   'av_k',nan(settings.loop.n, 1), 'av_trms',nan(settings.loop.n, 1),...
   'fres',nan(settings.loop.n, 1), 'bw',nan(settings.loop.n, 1),...
   'n',nan(settings.loop.n, 1)), results.channel_rt, 'UniformOutput', false);
if ~settings.specials.probemode
   %     channel parameters Reader -> Tag (several receivers)
   results.channel_tr = cell(settings.readerpool.n_nonvirt, settings.tagpool.n);
   results.channel_tr = cellfun(@(x) struct(...
      'pdp_los',nan(settings.loop.n, settings.readerpool.n_nonvirt), 'pdp_sum',nan(settings.loop.n, settings.readerpool.n_nonvirt),...
      'k',nan(settings.loop.n, settings.readerpool.n_nonvirt), 'trms',nan(settings.loop.n, settings.readerpool.n_nonvirt),...
      'tmin',nan(settings.loop.n, settings.readerpool.n_nonvirt), 'tmax',nan(settings.loop.n, settings.readerpool.n_nonvirt),...
      'av_k',nan(settings.loop.n, settings.readerpool.n_nonvirt), 'av_trms',nan(settings.loop.n, settings.readerpool.n_nonvirt),...
      'fres',nan(settings.loop.n, settings.readerpool.n_nonvirt), 'bw',nan(settings.loop.n, settings.readerpool.n_nonvirt),...
      'n',nan(settings.loop.n, settings.readerpool.n_nonvirt)), results.channel_tr, 'UniformOutput', false);
   %     tag RX power
   results.rx_pwr = nan(settings.loop.n, settings.readerpool.n_nonvirt, settings.tagpool.n);
   %     power levels recorded by tag
   results.tag_pwr = cell(settings.tagpool.n, 1);
   results.tag_pwr = cellfun(@(x) struct('pav',nan(settings.loop.n,settings.readerpool.n_nonvirt),...
      'pin',nan(settings.loop.n, settings.readerpool.n_nonvirt), 'pic',nan(settings.loop.n, settings.readerpool.n_nonvirt),...
      'fnc',nan(settings.loop.n, settings.readerpool.n_nonvirt), 'vdda',nan(settings.loop.n, settings.readerpool.n_nonvirt)), results.tag_pwr, 'UniformOutput', false);
   %     MFCW ranging
   results.mfcw = cell(settings.readerpool.n_nonvirt, 1);
   for i = 1 : settings.readerpool.n_nonvirt
      numcomb = max(floor((settings.readerpool.mfcw.nc{i}+2)/2)*2, size(settings.reader{i}.ranging.mfcw_calcdist.c_ord,1));
      results.mfcw{i} = struct(...
         'dist_hat',nan(settings.loop.n, settings.readerpool.n_nonvirt, numcomb),...
         'avg_c_mi',nan(settings.loop.n, settings.readerpool.n_nonvirt, numcomb),...
         'avg_c_i', nan(settings.loop.n, settings.readerpool.n_nonvirt, numcomb),...
         'avg_c_im',nan(settings.loop.n, settings.readerpool.n_nonvirt, numcomb));
   end
   clear('numcomb');
end

% start calculation
LoopInd = 0; eta = get_eta();
while LoopInd < settings.loop.n
   LoopInd = LoopInd + 1;
   
   % make sure the loop can be restarted if anything happens
   try
      % time?
      eta = get_eta(eta, LoopInd, settings.loop.n);
      %     output header if not in probemode (would slow down simulation considerably)
      headline('\n');
      headline('*******************************************************************************************************');
      headline('* Looping Positions and Performing Ranging (MFCW)');
      headline('* %s, run %i of %i (on %s, pid %.0f)', settings.suffix, LoopInd, settings.loop.n, settings.hostname, settings.pid );
      headline('* Time: %2i days %2i hrs %2i min %2i sec', eta.time_past);
      if LoopInd > 1
         headline('* ETA:  %2i days %2i hrs %2i min %2i sec (End: %s)', eta.time_left, datestr(eta.datenum_stop, 0));
      end
      %     send ETA email (if requested, this is a "long" simulation and the proper loopindex has been reached)
      if settings.emails.eta_sendmail && (LoopInd == settings.emails.eta_atloop) && ...
            ( eta.sec_left / 3600 >= settings.emails.eta_mintime )
         send_email(sprintf('ETA: %s for %s', datestr(eta.datenum_stop, 0), settings.suffix), '');
      end
      % renew Kerberos ticket and AFS authentication token if requested
      if settings.renew_kticket
         [status, system_msg] = renew_ticket();
         if status == 0
            if isempty(system_msg)
               headline('\nSuccessfully renewed Kerberos ticket and AFS authentication token.');
            else
               headline('\nSuccessfully renewed Kerberos ticket and AFS authentication token [%s].', system_msg);
            end
         else
            err('Kerberos ticket and/or AFS authentication token renewal failed:\n****************************************\n%s', system_msg);
         end
         clear('status', 'system_msg');
      end
      
      % modify position of tags
      headline('\nModifying Positions');
      for i = 1 : settings.tagpool.n
         settings.tagpool.pos{i} = settings.loop.pos_t{i}(LoopInd, :);
      end
      %     set channels to new distances
      settings.channel_rt = channel_newpos(settings.channel_rt, ...
         settings.readerpool.pos(1:settings.readerpool.n), settings.tagpool.pos(1:settings.tagpool.n), stdsettings.channel);
      settings.channel_tr = channel_newpos(settings.channel_tr, ...
         settings.tagpool.pos(1:settings.tagpool.n), settings.readerpool.pos(1:settings.readerpool.n), stdsettings.channel);
      %     calculate expected distances (R1->T1->R1, R1->T1->R2, and so on); round to sampling interval
      for r1 = 1 : settings.readerpool.n_nonvirt
         for r2 = 1 : settings.readerpool.n_nonvirt
            for t1 = 1 : settings.tagpool.n
               results.expdist{LoopInd}(r1, r2, t1) = round( 1/2 *... % 1/2: we estimate distance r->t, not r->t->r below
                  (settings.channel_rt{r1,t1}.largescale.dist + settings.channel_tr{t1,r2}.largescale.dist) / ... 
                  settings.c*settings.fs ) * settings.c/settings.fs;
            end
         end
      end
      %     clean up variables
      clear('r1', 'r2', 't1');
      
      % select active reader (all other readers are passive)
      for ActRdr = 1 : settings.readerpool.n_nonvirt
         % ... only one active reader
         headline('\n\n****************************************');
         headline('* Active Reader: %i', ActRdr);
         reader_carrier = cell(1, settings.readerpool.n_nonvirt);
         reader_carrier = cellfun(@(x) [], reader_carrier, 'UniformOutput',false);
         reader_carrier{ActRdr} =  reader_carrier_base{ActRdr};
         
         % channel R -> T
         headline('\nChannel Reader -> Tag');
         %     create identical up/downlink channels to/from active reader
         %     ... [???] more elegant implementation without loops
         %     ... clock * rand should be ok even if rand sequence is accidentally nonrandom (identical seeds)
         if settings.readerpool.id_ch{ActRdr}
            for j = 1 : settings.tagpool.n
               settings.channel_rt{ActRdr,j}.smallscale.seed = round(mod(sum(1e7*clock), 2^32-1) * rand);
               settings.channel_tr{j,ActRdr}.smallscale.seed = settings.channel_rt{ActRdr,j}.smallscale.seed;
            end
         end
         %     apply channels
         [tag_carrier, channel_stat] = channel_main(reader_carrier, settings.channel_rt);
         %     record some PDP and smallscale parameters
         for j = 1 : settings.tagpool.n
            results.channel_rt{ActRdr,j}.pdp_los(LoopInd) = channel_stat{ActRdr,j}.pdp_los; % estimated PDP LOS component
            results.channel_rt{ActRdr,j}.pdp_sum(LoopInd) = channel_stat{ActRdr,j}.pdp_sum; % estimated PDP sum
            results.channel_rt{ActRdr,j}.k      (LoopInd) = channel_stat{ActRdr,j}.k; % estimated Ricean K-Factor
            results.channel_rt{ActRdr,j}.trms   (LoopInd) = channel_stat{ActRdr,j}.trms; % estimated RMS delay spread
            results.channel_rt{ActRdr,j}.tmin   (LoopInd) = channel_stat{ActRdr,j}.tmin; % minimum delay (LOS)
            results.channel_rt{ActRdr,j}.tmax   (LoopInd) = channel_stat{ActRdr,j}.tmax; % minimum delay (last significant NLOS)
            results.channel_rt{ActRdr,j}.cir_d  {LoopInd} = channel_stat{ActRdr,j}.cir_d; % channel impulse response: delay vector in s
            results.channel_rt{ActRdr,j}.cir_g  {LoopInd} = channel_stat{ActRdr,j}.cir_g; % channel impulse response: gain vector
            if isfield(channel_stat{ActRdr,j}, 'av_k')
               results.channel_rt{ActRdr,j}.av_k   (LoopInd) = channel_stat{ActRdr,j}.av_k; % (theoretical) average Ricean K-Factor
               results.channel_rt{ActRdr,j}.av_trms(LoopInd) = channel_stat{ActRdr,j}.av_trms; % (theoretical) average RMS delay spread
               results.channel_rt{ActRdr,j}.fres   (LoopInd) = channel_stat{ActRdr,j}.fres; % frequency resolution
               results.channel_rt{ActRdr,j}.bw     (LoopInd) = channel_stat{ActRdr,j}.bw; % bandwidth
               results.channel_rt{ActRdr,j}.n      (LoopInd) = channel_stat{ActRdr,j}.n; % number of taps
            end
         end
         
         % field probe mode: nothing else to do in this loop (channels are recorded)
         if settings.specials.probemode
            continue;
         end
         
         % record RX power levels (remove intial transient); tag power levels are recorded below
         delay_s_max = round( max(cellfun(@(x) x.tmax(LoopInd)/2, results.channel_rt)) * settings.fs );
         results.rx_pwr(LoopInd, ActRdr, :) = cellfun(@(x) var(x(delay_s_max+1:end)), tag_carrier);        
         
         % modulate RN16 or just reflect
         headline('\nModulation of RN16 (only active tags)');
         for i = 1 : settings.tagpool.n
            headline('   Tag %i', i);
            tx_temp = tag_main('tx_data_active', tag_carrier{i}, settings.tag{i});
            tag_modcarrier{i} = tx_temp.txsignal; % signal
            settings.tag{i}   = tx_temp.settings; % modifications to vdda
            % record power levels
            results.tag_pwr{i}.pav (LoopInd, ActRdr) = tx_temp.pav_hat; % available power
            results.tag_pwr{i}.pin (LoopInd, ActRdr) = tx_temp.pin_hat; % input power (not reflected)
            results.tag_pwr{i}.pic (LoopInd, ActRdr) = tx_temp.pic_hat; % chip input power
            results.tag_pwr{i}.fnc (LoopInd, ActRdr) = tx_temp.powered; % functional (powered)?
            results.tag_pwr{i}.vdda(LoopInd, ActRdr) = tx_temp.vdda;    % power supply voltage
         end
         
         % channel T -> R
         headline('\nChannel Tag -> Reader');
         [reader_modcarrier, channel_stat] = channel_main(tag_modcarrier, settings.channel_tr);
         %     record some PDP and smallscale parameters
         for i = 1 : settings.tagpool.n
            for j = 1 : settings.readerpool.n_nonvirt
               results.channel_tr{i,j}.pdp_los(LoopInd, ActRdr) = channel_stat{i,j}.pdp_los; % estimated PDP LOS component
               results.channel_tr{i,j}.pdp_sum(LoopInd, ActRdr) = channel_stat{i,j}.pdp_sum; % estimated PDP sum
               results.channel_tr{i,j}.k      (LoopInd, ActRdr) = channel_stat{i,j}.k; % estimated Ricean K-Factor
               results.channel_tr{i,j}.trms   (LoopInd, ActRdr) = channel_stat{i,j}.trms; % estimated RMS delay spread
               results.channel_tr{i,j}.tmin   (LoopInd, ActRdr) = channel_stat{i,j}.tmin; % minimum delay (LOS)
               results.channel_tr{i,j}.cir_d  {LoopInd, ActRdr} = channel_stat{i,j}.cir_d; % channel impulse response: delay vector in s
               results.channel_tr{i,j}.cir_g  {LoopInd, ActRdr} = channel_stat{i,j}.cir_d; % channel impulse response: gain vector
               if isfield(channel_stat{i,j}, 'av_k')
                  results.channel_tr{i,j}.av_k   (LoopInd, ActRdr) = channel_stat{i,j}.av_k; % (theoretical) average Ricean K-Factor
                  results.channel_tr{i,j}.av_trms(LoopInd, ActRdr) = channel_stat{i,j}.av_trms; % (theoretical) average RMS delay spread
                  results.channel_tr{i,j}.fres   (LoopInd, ActRdr) = channel_stat{i,j}.fres; % frequency resolution
                  results.channel_tr{i,j}.bw     (LoopInd, ActRdr) = channel_stat{i,j}.bw; % bandwidth
                  results.channel_tr{i,j}.n      (LoopInd, ActRdr) = channel_stat{i,j}.n; % number of taps
               end
            end
         end
         
         % channel R -> R ("feedback")
         headline('\nChannel Reader -> Reader');
         reader_modcarrier = cellfun(@(a,b) a+b, reader_modcarrier, cellfun(@(a,b) a(1:length(b)),...
            channel_main(reader_carrier, settings.channel_rr), reader_modcarrier, 'UniformOutput',false),...
            'UniformOutput',false);
         
         % demodulation and decoding
         headline('\nDemodulation');
         for i = 1 : settings.readerpool.n_nonvirt
            headline('   Reader %i', i);
            reader_baseband{i} = reader_main('rx', reader_modcarrier{i}, settings.reader{i});
         end
         
         % "memory management" and variable cleanup
         %     free memory (large vectors)
         clear('tag_carrier', 'tag_modcarrier', 'reader_modcarrier', 'tx_temp');
         %     clean up variables
         clear('nfft', 'max_clen_s');
         
         % ranging (multifrequency continuous-wave radar)
         headline('\nMFCW Ranging');
         for i = 1 : settings.readerpool.n_nonvirt
            headline('   Reader %i', i);
            % MFCW
            mfcw_est = reader_main('mfcw_est', reader_baseband{i}, settings.reader{i});
            % record
            results.mfcw{i}.dist_hat(LoopInd, ActRdr, 2:end) = cell2mat(mfcw_est.dist);
            results.mfcw{i}.avg_c_mi(LoopInd, ActRdr, :)     = mfcw_est.avg_c_mi;
            results.mfcw{i}.avg_c_i (LoopInd, ActRdr, :)     = mfcw_est.avg_c_i ;
            results.mfcw{i}.avg_c_im(LoopInd, ActRdr, :)     = mfcw_est.avg_c_im;
            % output expected/estimated
            msg('Estimated distances:');
            for j = 2 : settings.reader{i}.ranging.mfcw_calcdist.nc + 1
               msg('   %.0fFCW estimate (no model): %8.3f m', j, mfcw_est.dist{j});
            end
            msg('Perfect estimates (distance (R->T->R)/2):');
            for j = 1 : settings.tagpool.n
               msg('   Reader %.0f -> Tag %.0f -> Reader %.0f: %8.3f m', ActRdr, j, i, results.expdist{LoopInd}(ActRdr, i, j));
            end
         end
      end
      
   % ... this requires user interaction
   catch ME
      % automatic mode: option to repeat the last loop (if error can be removed)
      if settings.err_rep
         disp(sprintf('\n\n********* Simulation halted due to an error. Press a button to restart previous loop, [Strg]+[C] to abort. *********\n'));
         % try to send an email
         try
            disp('   Trying to send an eMail...');
            send_email('User interaction required (caught an error).',...
               sprintf('The following error was caught by %s on %s; the simulation has been halted.\n\n******************************\n\n%s',...
               mfilename, datestr(now, 0), getReport(ME, 'extended', 'hyperlinks','off')));
         catch %#ok<CTCH>
            disp('      ... nope, didn''t work. This seems to be serious.');
         end
         disp(sprintf('\n*********\n%s', getReport(ME, 'extended', 'hyperlinks','off')));
         disp(sprintf('********* Simulation halted due to an error. Press a button to restart previous loop, [Strg]+[C] to abort. *********'));
         % pause simulation and restart last loop if requested
         pause;
         LoopInd = LoopInd - 1;
      % "debug" mode: throw the error, enter debug mode if activated
      else
         rethrow(ME);
      end
   end
      
end

% save results and stop diary
results.misc.simtime = toc;
system(sprintf('rm -f %s', fullfile(settings.folders.root, settings.folders.data, settings.suffix)));% delete old file
save('-v7.3', fullfile(settings.folders.root, settings.folders.data, settings.suffix), 'settings', 'results');
disp(sprintf('\n\n\n%s, time: %i days %i hrs %i min %i sec', settings.suffix, sec2dhms(results.misc.simtime)));
diary off
