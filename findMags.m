function [mag, phi, freq] = findMags(D, Freq, measure_window)
fprintf('Finding magnitudes for file %s\n', D.file);

% check vs. csv data
if(D.nsweeps == length(Freq))
	% calculate range of bins (samples) for measurement (skip ramp on/off)
	measure_bins = ms2bin(measure_window(1), D.cal.Fs)	...
							: ms2bin(measure_window(2), D.cal.Fs);
	% allocate mag, phi arrays
	mag = zeros(D.nsweeps, 1);
	phi = zeros(D.nsweeps, 1);
	freq = Freq;
	for n = 1:D.nsweeps
		% compute pll magnitude and phase at this frequency
		[mag(n), phi(n)] = fitsinvec(D.data{n}(measure_bins), ...
										1, D.cal.Fs, Freq(n));
	end

else
	% check vs. csv data
	fprintf('D data %s has %d sweeps\n', D.file, D.nsweeps);
	fprintf('CSV file has %d stimuli\n', length(Freq));
	% auto detect tones from raw unattenuated data
	% tolerance (in Hz) for finding peak frequency in autodetect mode
	fprintf('automatically determining test frequency\n');

	measure_bins = ms2bin(measure_window(1), D.cal.Fs)	...
								: ms2bin(measure_window(2), D.cal.Fs);
	% allocate mag, phi arrays
	mag = zeros(D.nsweeps, 1);
	phi = zeros(D.nsweeps, 1);
	freq = zeros(D.nsweeps, 1);
	for n = 1:D.nsweeps
		% get spectrum of data
		[~, ~, fmax, ~] = daqdbfft(D.data{n}', D.cal.Fs, length(D.data{n}));
		% store automagically determined frequency
		freq(n) = fmax;
		% compute pll magnitude and phase at this frequency
		[mag(n), phi(n)] = fitsinvec(D.data{n}(measure_bins), ...
											1, D.cal.Fs, freq(n));

	end

end