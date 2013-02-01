% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - protocol ("large scale" link-) timings as defined in EPCglobal Class1 Gen2 v1.0.9
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
% version = epc_protocol_timings()
%    Just returns the version number (string).
% timings = epc_protocol_timings(reader_modulation, reader_command)
%    Returns basic EPCglobal Class1 Gen2 protocol timings for the given setup in READER_MODULATION and
%    READER_COMMAND.
%    Note that only query/RN16-timings are implemented so far.
%
%
% ***** Function definition *****
% timings = epc_protocol_timings(reader_modulation, reader_command)
%    reader_modulation   struct containing reader modulator settings
%    reader_command      struct containing reader command settings
%
%    timings   struct with substructs containing basic protocol timings
%       .ft             frequency tolerance of the tag subcarrier (nominal temp, not including variations)
%       .rtcal          length of data-0 plus data-1 in s
%       .tpri           period of a single subcarrier cycle
%       .idle           idle times between commands/replies (min==max if only one time was specified)
%          .t1             time reader command -> tag reply [min, max] in us
%          .t2             time tag reply -> reader command [min, max] in us
%          .t3             time reader command -> reader command in case there was no tag reply [min, max] in us
%          .t4             time reader command -> reader command [min, max] in us
%       .reader         timing details concerning the reader command
%          .len            length of reader command in us
%       .tag            timing details concerning the tag reply
%          .len            length of tag reply [min, max] in us
%          .pilot          position of pilot tone within tag reply [start, stop] times .tag.len
%       .comm           timing overview for the entire communication
%          .reader_start   beginning of reader command [min, max] in us
%          .reader_stop    end of reader command  [min, max] in us
%          .tag_start      beginning of tag reply  [min, max] in us
%          .tag_stop       end of tag reply  [min, max] in us
%          .pilot_start    beginning of pilot tone (unmodulated subcarrier) within tag reply  [min, max] in us
%          .pilot_stop     end of pilot tone (unmodulated subcarrier) within tag reply  [min, max] in us
%          .len            overall length of one round command-reply  [min, max] in us
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

function timings = epc_protocol_timings(reader_modulation, reader_command)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   timings = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% frequency tolerances (lookup-table)

switch(reader_modulation.dr)
   case 64/3
      if reader_modulation.lf == 640e3
         timings.ft = 0.15;
      elseif 320e3 < reader_modulation.lf && reader_modulation.lf < 640e3
         timings.ft = 0.22;
      elseif reader_modulation.lf == 320e3
         timings.ft = 0.1;
      elseif 256e3 < reader_modulation.lf && reader_modulation.lf < 320e3
         timings.ft = 0.12;
      elseif reader_modulation.lf == 256e3
         timings.ft = 0.1;
      elseif 160e3 <= reader_modulation.lf && reader_modulation.lf < 256e3
         timings.ft = 0.1;
      elseif 107e3 <= reader_modulation.lf && reader_modulation.lf < 160e3
         timings.ft = 0.07;
      elseif  95e3 <= reader_modulation.lf && reader_modulation.lf < 107e3
         timings.ft = 0.05;
      else
         critwarn('Unsupported combination DR / LF. Setting frequency tolerance to FT=0');
         timings.ft = 0;
      end
   case 8
      if 320e3 < reader_modulation.lf && reader_modulation.lf <= 465e3
         timings.ft = 0.19;
      elseif reader_modulation.lf == 320e3
         timings.ft = 0.1;
      elseif 256e3 < reader_modulation.lf && reader_modulation.lf < 320e3
         timings.ft = 0.12;
      elseif reader_modulation.lf == 256e3
         timings.ft = 0.1;
      elseif 160e3 < reader_modulation.lf && reader_modulation.lf < 256e3
         timings.ft = 0.1;
      elseif reader_modulation.lf == 160e3
         timings.ft = 0.07;
      elseif 107e3 < reader_modulation.lf && reader_modulation.lf < 160e3
         timings.ft = 0.07;
      elseif  40e3 < reader_modulation.lf && reader_modulation.lf < 107e3
         timings.ft = 0.04;
      else
         critwarn('Unsupported combination DR / LF. Setting frequency tolerance to FT=0');
         timings.ft = 0;
      end
   otherwise
      critwarn('Unsupported Divide Ratio. Setting frequency tolerance to FT=0.');
      timings.ft = 0;
end


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% generic
timings.rtcal = reader_modulation.rtcal;
timings.tpri  = 1 / reader_modulation.lf;


% idle periods [min, max] assuming minimum timings (max==min in case max==Inf is allowed)
%     T1
timings.idle.t1 = ones(1,2) * max(timings.rtcal, 10*timings.tpri);
timings.idle.t1 = timings.idle.t1 .* [(1-timings.ft), (1+timings.ft)] + [-2e-6, 2e-6];
%     T2
timings.idle.t2 = [3, 20] * timings.tpri;
%     T4
timings.idle.t4 = [1, 1] * 2*timings.rtcal;
%     T3
timings.idle.t3 = [1, 1] * max(0, timings.idle.t4(1)-timings.idle.t1(1));


% reader command
timings.reader.len = reader_modulation.length;


% tag reply 
%     get tag data vector (bits at subcarrier level)
tag.m     = reader_command.m;
tag.trext = reader_command.trext;
switch lower(reader_command.cmd)
   case 'query'
      tag.nbits = 16;
   otherwise
      err('Reader command "%s" not supported.', reader_command.cmd);
end
tag.data = tag_encoding(zeros(1,tag.nbits), tag, false);
%     length of tag reply
timings.tag.len = length(tag.data)/2 * timings.tpri * [1-timings.ft, 1+timings.ft];
%     position of pilot tone (or initial subcarrier signal) [start,end] times timings.tag.len
if reader_command.trext % pilot tone active
   if reader_command.m == 1 % FM0
      timings.tag.pilot = [0, 12] * reader_command.m * 2 / length(tag.data) ;
   else % Miller
      timings.tag.pilot = [0, 16] * reader_command.m * 2 / length(tag.data);
   end
else % no (additional) pilot tone
   if reader_command.m == 1 % FM0
      timings.tag.pilot = [0, 0]; % [start, end] times timings.tag.len
   else % Miller
      timings.tag.pilot = [0, 4] * reader_command.m  * 2 / length(tag.data);
   end
end


% [min, max] of reader command, tag reply and pilot tone
timings.comm.reader_start = timings.idle.t4;
timings.comm.reader_stop  = timings.comm.reader_start + timings.reader.len;
timings.comm.tag_start    = timings.comm.reader_stop  + timings.idle.t1;
timings.comm.tag_stop     = timings.comm.tag_start + timings.tag.len;
timings.comm.pilot_start  = timings.comm.tag_start + timings.tag.len * timings.tag.pilot(1);
timings.comm.pilot_stop   = timings.comm.tag_start + timings.tag.len * timings.tag.pilot(2);
%     overall length
timings.comm.len = timings.comm.tag_stop + timings.idle.t2;

