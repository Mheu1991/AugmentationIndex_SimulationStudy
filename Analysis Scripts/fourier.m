function[coef]=fourier(signal,sf,ncof,f);
%FOURIER FOURIER computes the first NCOF Fourier coefficients of a signal 
%	SIGNAL, with basic frequency F and samplefrequency SF.
%	The signal should meet the following requirements:
%
%	1:	The signal must be periodic.
%	
%	2:	The given sample should be an integer number of periods.
%
%	call:
%
%	[coef] = fourier(signal,sf,ncof,f);
%
%	in:	signal:	guess what !
%		sf:	sampling frequency			[Hz]
%		ncof:	number of required fourier coefficients [-]
%		f:	basic frequency of the signal		[Hz]
%
%	out:	coef:	set of complex Fourier coefficients

coef	= zeros(1,ncof+1);
coef(1)	= mean(signal);
%
% Eliminate DC-component:
%
signal	= signal - mean(signal);
signal	= signal(:);
%
% Compute the number of periods in this sample:
%
T 	= 1/f;
dt 	= 1/sf;
tmax	= length(signal) * dt;
np	= round(tmax / T);
%
% Compute the Fourier coefficients:
%
t 	= linspace(1/sf,length(signal)/sf,length(signal));
t	= t(:);
%
for cnt1 = 2 : ncof + 1,
	coef(cnt1) = 2*f/np * ...
		sum(dt * signal .* exp(-2 * pi * j * f * (cnt1-1) .* t));
end;
