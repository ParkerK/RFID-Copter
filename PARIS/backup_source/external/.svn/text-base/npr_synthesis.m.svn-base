%NPR_SYNTHESIS Near perfect reconstruction synthesis
%  Y = NPR_SYNTHESIS(C,X) synthesizes a wideband signal Y from a
%  number of subbands stored in X. Each subband is a row in X.
%  C is a two dimensional array, containing the filter coefficients.
%  The number of rows in X must be twice the number of rows in C. 
%
% (c) 2007 Wessel Lubberhuizen
%  All rights reserved.
%  
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%  
%      * Redistributions of source code must retain the above copyright 
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright 
%        notice, this list of conditions and the following disclaimer in 
%        the documentation and/or other materials provided with the distribution
%        
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
%
% modified by Daniel Arnitz 06/2009: fft replaced by ifft

function [y] = npr_synthesis(coeff,x)

% number of channels
N=size(coeff,1);

% number of slices
M=size(x,2);

% split into even and odd channels
x = reshape(x,2,N,M);

y1 = squeeze(x(1,:,:));
y2 = squeeze(x(2,:,:));

% apply dft
y1 = ifft(y1,[],1)*N;
y2 = ifft(y2,[],1)*N;

% apply channel filters
for i=1:N
    y1(i,:) = filter(coeff(i,:),1,y1(i,:));
    y2(i,:) = filter(coeff(i,:),1,y2(i,:));
end

% apply frequency shift
for i=1:N;
    y2(i,:) = y2(i,:) * exp(-sqrt(-1)*pi*(i-1)/N);
    y2(i,2:2:M) = -y2(i,2:2:M);
end

% combine filter results
y = reshape(y1-y2,1,M*N);
