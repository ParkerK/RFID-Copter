% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% load .mat-file and print version information
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
% data = loadmat(filename)
%    Loads the .mat-file addressed by FILENAME and prints version information.
% data = loadmat(filename, leadin)
%    Loads the .mat-file addressed by FILENAME and prints version information. Does not use MESSAGE to
%    display the version information, but SPRINTF with LEADIN printed before the version information.
%    Leadin will be ignored if it is not a string/character.
%
%
% ***** Global Variables *****
% globalsettings
%    .core.mat_en_cache   enable matfile cache; load will only be called once per file if enabled
%    .core.mat_stubborn   ignore missing files; return empty struct if file is missing
% loadmat_cache   cache for matfile data (has to be enabled); field names are file names without dots
%
%
% ***** Interface definition *****
% function data = loadmat(filename)
%    filename    filename of .mat-file to load
%    leadin      (optional) string to print before the version information
%
%    data        loaded data including version info header
% 
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
% ? hash for generation of fieldname
%
% *******************************************************************************************************

function data = loadmat(filename, leadin)
% global variables
global globalsettings
global loadmat_cache;

% default for fromcache_msg
fromcache_msg = '';

% CACHE ENABLED
if globalsettings.core.mat_en_cache
   % does the cache exist yet?
   if ~isstruct(loadmat_cache)
      loadmat_cache = struct();
   end

   % create field name for cache (struct does not like mathematical signs for an obvious reason)
   fieldname = filename;
   fieldname( ((fieldname=='.') + (fieldname=='-') + (fieldname=='+') +...
               (fieldname=='*') + (fieldname=='/') + (fieldname=='\'))>0 ) = [];
   
   % load from cache
   if isfield(loadmat_cache, fieldname)
      data = loadmat_cache.(fieldname);
      fromcache_msg = ' from cache'; % for message
   % load file from disk and write to cache
   else 
      try 
         data = load(filename);
      catch ME
         data = struct();
      end
      loadmat_cache.(fieldname) = data;
   end

% CACHE DISABLED
else
   try data = load(filename); catch ME; data = struct(); end
end

% output and exception handling
%     unable to load this file and ignoring of errors turned off
if exist('ME', 'var') && ~globalsettings.core.mat_stubborn
   rethrow(ME);
   
%     if we were unable to load this file
elseif isempty(fieldnames(data))
   if nargin==2 && ischar(leadin) % custom
      disp(sprintf('%sWARNING: Unable to load %s.mat; returning a 1-by-1 struct with no fields%s',...
         leadin, filename, fromcache_msg));
   else % use message system
      critwarn('Unable to load %s.mat; returning a 1-by-1 struct with no fields%s',...
         filename, fromcache_msg);
   end
   
%     version information if the loading was successfull
else
   if nargin==2 && ischar(leadin) % custom
      disp(sprintf('%sUsing %s.mat%s (created by %s)', leadin, data.matfilename, fromcache_msg, data.createdby));
   else % use message system
      msg('Using %s.mat%s (created by %s)', data.matfilename, fromcache_msg, data.createdby);
   end
end

