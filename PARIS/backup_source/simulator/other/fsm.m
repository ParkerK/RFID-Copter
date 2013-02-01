% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - general finite state machine implementation
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
% version = fsm()
%    Just returns the version number (string).
% carrier = fsm(input, transitions, basisfcns, initialstate, symbolset)
%    Returns the output vector of the state machine with input vector INPUT.
%
%
% ***** Function definition *****
% output = fsm(input, transitions, basisfcns, initialstate, symbolset)
%    input         input signal
%    transitions   transition matrix of size length(symbolset) x states
%                  [state1, state2, ...]
%                  in case the symbol is -1, transitions(2) for 0 and so on
%    basisfcns     basis function matrix of size length(basisfunctions) x states
%                  [state1, state2, ...]
%    initialstate  state at first input value input(1) (states are named
%                  1,2,3,... i.e. column number of transitions or basisfcns)
%    symbolset     optional list of allowed input symbols (default: [0;1] for binary)
%
%    output        output signal (encoded by FSM)
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

function output = fsm(input, transitions, basisfcns, initialstate, symbolset)
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

% number of input variables
if nargin < 4
    symbolset = [0; 1]; % default: binary
end
if nargin < 3
    criterr('Not enough input arguments.')
end

% length of input vector
if isempty(input)
    criterr('Length of input signal is zero.')
end

% size of transitions and basisfcns
if size(transitions, 1) == 0 || size(transitions, 2) == 0
    err('Matrix transitions has zero size');
end
if size(basisfcns, 1) == 0 || size(basisfcns, 2) == 0
    err('Matrix basisfcns has zero size');
end
if size(transitions, 2) ~= size(basisfcns, 2)
    err('Matrices transitions and basisfcns must have same amount of columns (states)');
end

% length of symbolset
if length(symbolset) < 2
    warn('Less than two valid symbols defined: using binary symbolset [0,1] instead');
end

% make all input vectors column vectors
input = input(:);
symbolset = symbolset(:);


% *******************************************************************************************************
% preparations

% length of basis functions
N = size(basisfcns,1);

% prepeare output vector
output = zeros(N*length(input), 1);

% add dummy to input vector (just for statemachine termination)
input = [input; 0]; % add dummy (just for statemachine)

% initial state
state = initialstate;


% *******************************************************************************************************
% state machine

for i = 1 : length(input)-1
    
    % output (when entering state)
    output((i-1)*N+1:i*N) = basisfcns(:,state);
    
    % find out index of in(i) in symbolset
    [tf, index] = ismember(input(i+1), symbolset);
        % if in(i) is not a valid symbol
    if tf == 0
        err('Found invalid symbol in input.')
    end
        
    % transition to next state
    state = transitions(index, state);
end
