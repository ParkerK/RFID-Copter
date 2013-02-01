%NPR_ANALYSIS Near perfect reconstruction analysis
%  Y = NPR_ANALYSIS(COEFF,X) separates the input signal X into
%  subbands. Each subband is a row in X. COEFF is a two dimensional array, 
%  containing the filter coefficients. The number of rows in X will be
%  be twice the number of rows in COEFF. 
%
%  See also npr_synthesis, npr_coeff
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
% modified by Daniel Arnitz 06/2009: ifft replaced by fft


function y = npr_analysis(coeff,x)

% number of channels
N=size(coeff,1);

% number of slices
M=ceil(length(x)/N);

% create polyphase input signals
x1=reshape(x,N,M);
x2=x1;
for i=1:N;
    x2(i,:) = x2(i,:) * exp(sqrt(-1)*pi*(i-1)/N);
    x2(i,2:2:M) = -x2(i,2:2:M);
end
% apply channel filters
coeff = fliplr(coeff);
for i=1:N
    x1(i,:) = filter(coeff(i,:),1,x1(i,:));
    x2(i,:) = filter(coeff(i,:),1,x2(i,:));    
end

% apply dft
x1 = fft(x1,[],1)*N;
x2 = fft(x2,[],1)*N;

% assemble even and odd channels
y = [reshape(x1,1,N,M); reshape(x2,1,N,M)];
y = reshape(y,2*N,M);
