function[signal] = ifourier(coef, ns, f);
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
%	[signal] = ifourier(coef, ns, f);
%
%	in:	coef:	Fourier coefficients odf the signal (complex)
%		ns:	number of required samples		[-]
%		f:	basic frequency of the signal		[Hz]
%
%	out:	signal:	you already know ..

% build time axis:

t = linspace(1/(ns*f),1/f,ns);

% DC-component:

signal = ones(1,ns)*coef(1);

% harmonics:

for cnt = 2:length(coef),
  signal = signal + coef(cnt) * exp(2*pi*j*(cnt-1)*f.*t);
end;

signal = real(signal);