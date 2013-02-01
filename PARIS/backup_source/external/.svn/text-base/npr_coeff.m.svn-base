%NPR_COEFF generates near NPR filter bank coefficients
%  COEFF = NPR_COEFF(N,L) generates the filter coefficients
%  for a near perfect reconstruction filter bank.
%  The number of subbands will be N, and L is the number of 
%  filter taps used per subband. The output COEFF will have 
%  size (N/2,L).
%
%  The prototype is constructed starting with an equiripple
%  approximation to a 'root raised error function' shaped filter.
%
%  NPR_COEFF(N,L) with no output arguments plots the magnitude
%  and the prototype filter in the current figure window.
%
%  See also npr_analysis, npr_synthesis, npr_coeff
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

function coeff = npr_coeff(N,L,K)

% K determines the overlap (higher is less).
if ~exist('K')
    K = 7.5;
end
N = N / 2;

M = N;

if M>2
  % start with lower N,  because the
  % equiripple design will take far too much time,
  % or does not even converge for large N.
  M = 2;
end

F=linspace(0,1,512*L*M);
% The prototype is based on a root raised error function
A=sqrt(1-0.5*(1+erf(K*(M*F-0.5))));

% d=fdesign.arbmag('N,F,A',L*M-1,F,A);
% Hd = design(d,'equiripple');
% b = Hd.numerator;


% REPLACED
b = firls(L*M-1,F,A);



if N>M 
  Q = N/M;
  % Interpolate, using fourier interpolation method
  f1=fft(b,length(b));
  f2=zeros(1,L*N);
  % copy the lower frequency half. 
  n=0:L*M/2;
  f2(1+n)=f1(1+n);
  % if we want to make the impulse response symmetric in time,
  % it needs to be time-shifted a little.
  % we do it using a phase correction in the frequency domain
  f2(1+n)=f2(1+n).*exp(-sqrt(-1)*(Q-1)*pi*n/L/N);
  % recreate the upper part of the spectrum from the lower part
  f2(L*N-n)=conj(f2(2+n));
  % back to time domain
  b=real(ifft(f2));
end

% force the impulse response to be symmetric
% and hence linear phase
b = b + fliplr(b);

% normalize
b=b/sum(b);


if nargout==0
figure;
    subplot(4,1,1);
    t=0:length(b)-1;
    plot(t,b/max(b));
    axis([0 max(t) min(b)/max(b)*2 1.1]);
    title('prototype filter - impulse response');
    xlabel('sample');
    ylabel('value');
    grid on;
    subplot(4,1,2);
    f=2*(0:1023)/1024/N;
    h=freqz(b,1,pi*f);
    plot(f*N,20*log10(abs(h)));
    title('prototype filter - frequency response');
    xlabel('f / f_{channel}')
    ylabel('power (dB)');
    grid on;
    subplot(4,1,3);
    plot(f*N,180/pi*unwrap(angle(h)));
    title('prototype filter - frequency response');
    xlabel('f / f_{channel}')
    ylabel('phase (deg)');
    grid on;
    subplot(4,1,4);
    f=(0:1024)/1024/N;
    h1=freqz(b,1,pi*f);
    h2=fliplr(h1);
    plot(f*N,10*log10(abs(h1).^2.'+abs(h2).^2.'));
    title('power complementary check');
    xlabel('f / f_{channel}')
    ylabel('error (dB)');
    grid on;
else
    % reshape
    coeff = reshape(b,N,L);
end

