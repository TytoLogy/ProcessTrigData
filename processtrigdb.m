%------------------------------------------------------------------------
% script to process triggered capture data
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% Files:
%------------------------------------------------------------------------
% 
% Direct output from datawave (output from PA5 attenuator):
% 	DW_Atten_FullList_20190930_1.bin         
% 	DW_Atten_FullList_20190930_1.mat         
% 
% Direct output from datawave (output from datawave D/A, no attenuation):
% 	DW_Raw_FullList_20190930_1.bin           
% 	DW_Raw_FullList_20190930_1.mat
% 
% NICal Calibration files:
% 	FilterOnly_Fc100Khz_1V_5k-115k.cal       
% 	LCY_5k-80k_1V_20dBatten.cal              
% 	LCY_5k-80k_1V_20dBatten.fig              
% 	LCY_5k-80k_1V_PA5_20dBatten.cal          
% 
% Recorded output data using calibration mic:
% 	LCY_Calibration_FullList_20190930_1.bin  
% 	LCY_Calibration_FullList_20190930_1.mat  
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Revisions:
%	2 Oct 2019 (SJS): modified for use with I. Kristaponyte data
%	3 Oct 2019 (SJS): processing data
%-----------------------------------3-------------------------------------
%------------------------------------------------------------------------



%% Settings
% file with list of test information
infofile = 'calibrate_full_list_atten_30Higher.csv';

% path and name of triggered data file (.bin)
datapath = '/Users/sshanbhag/Work/Data/Audio/Calibration/Inga/30Sep2019';
% data files
LCYfile = 'LCY_Calibration_FullList_20190930_1.bin';
DWrawfile = 'DW_Raw_FullList_20190930_1.bin';
DWattfile = 'DW_Atten_FullList_20190930_1.bin';
% filtering high pass cutoff (Hz);
HP.Fc = 3500;
% filtering low pass cutoff (Hz);
LP.Fc = 80000;
% filter order
HP.forder = 5;
LP.forder = 5;

% windowing to eliminate transients at onset, offset
ramp_ms = 0.5;

% measurement window (ms)
measure_window = [5 95];

%% read in frequency, attenuation, desired level and scaling from file
% columns are:
% trial,freq,mV,atten,dB_SPL
% read in using csvread, skip first row (header)
listdata = csvread(infofile, 1, 0);
% 
Trial = listdata(:, 1);
Freq = listdata(:, 2);
mV = listdata(:, 3);
Atten = listdata(:, 4);
Level = listdata(:, 5);


%% read in raw data for LCY
LCY = readBinData(fullfile(datapath, LCYfile));

%% read in atten DW data
% DWa = readBinData(fullfile(datapath, DWattfile));

%% process LCY data

% get # of sweeps
LCY.nsweeps = length(LCY.data);
% check vs. csv data
fprintf('LCY data %s has %d sweeps\n', LCYfile, LCY.nsweeps);
fprintf('CSV file has %d stimuli\n', length(listdata));
% SPL conversion
VtoPa = 1 ./ ...
		(LCY.cal.Gain(1) * invdb(LCY.cal.MicGain(1)) * LCY.cal.MicSensitivity);

% first, filter the data

fprintf('Filtering %s data\n', LCYfile);
% Nyquist frequency
fnyq = LCY.cal.Fs/2;

% build a highpass filter for processing the data
[HP.fcoeffb, HP.fcoeffa] = butter(HP.forder, HP.Fc/fnyq, 'high');
% build a lowpass filter for processing the data
[LP.fcoeffb, LP.fcoeffa] = butter(LP.forder, LP.Fc/fnyq, 'low');

% loop through sweeps, apply ramp and filter
figure(1)
dt = 1/LCY.cal.Fs;
for n = 1:LCY.nsweeps
	% plot raw data
	subplot(2, 1, 1)
	tvec = 1000 * dt * (0:(length(LCY.data{n}) - 1));
	plot(tvec, LCY.data{n});
	ylabel('Raw (V)');
	title(sprintf('Sweep %d', n));
	% filter data
	% apply short ramp and highpass filter
	tmp = filtfilt(HP.fcoeffb, HP.fcoeffa, ...
					sin2array(LCY.data{n}', ramp_ms, LCY.cal.Fs));
	% apply lowpass filter
	tmp = filtfilt(LP.fcoeffb, LP.fcoeffa, tmp);
	% plot filtered data
	subplot(212)
	plot(tvec, tmp);
	xlabel('Time (ms)');
	ylabel('Filtered (V)');
	drawnow
end

%% read in raw DW data
DWrawfile = 'DW_RawfromDataWave_FullList_20191003_1.bin'
DWr = readBinData(fullfile(datapath, DWrawfile));

%% process raw, unattenuated data

% get # of sweeps
DWr.nsweeps = length(DWr.data);
% check vs. csv data
fprintf('DWr data %s has %d sweeps\n', DWrawfile, DWr.nsweeps);
fprintf('CSV file has %d stimuli\n', length(listdata));

% first, filter the data
fprintf('Filtering DWr data\n');
% Nyquist frequency
fnyq = DWr.cal.Fs/2;
% build a highpass filter for processing the data
[HP.fcoeffb, HP.fcoeffa] = butter(HP.forder, HP.Fc/fnyq, 'high');
% build a lowpass filter for processing the data
[LP.fcoeffb, LP.fcoeffa] = butter(LP.forder, LP.Fc/fnyq, 'low');

% loop through sweeps, apply ramp and filter
figure(1)
dt = 1/DWr.cal.Fs;
for n = 1:DWr.nsweeps
	% plot raw data
	subplot(2, 1, 1)
	tvec = 1000 * dt * (0:(length(DWr.data{n}) - 1));
	plot(tvec, DWr.data{n});
	ylabel('Raw (V)');
	title(sprintf('DWr Sweep %d', n));
	% filter data
	% apply short ramp and highpass filter
	tmp = filtfilt(HP.fcoeffb, HP.fcoeffa, ...
					sin2array(DWr.data{n}', ramp_ms, DWr.cal.Fs));
	% apply lowpass filter
	tmp = filtfilt(LP.fcoeffb, LP.fcoeffa, tmp);
	% plot filtered data
	subplot(212)
	plot(tvec, tmp);
	xlabel('Time (ms)');
	ylabel('Filtered (V)');
	drawnow
	DWr.data{n} = tmp;
end

%% auto detect tones from raw unattenuated data
% tolerance (in Hz) for finding peak frequency in autodetect mode
FreqDetectWidth = 21;
calfreq = 0;
measure_bins = ms2bin(measure_window(1), DWr.cal.Fs)	...
							: ms2bin(measure_window(2), DWr.cal.Fs);
for n = 1:DWr.nsweeps
	% get spectrum of data
	[tmpfreqs, tmpmags, fmax, magmax] = daqdbfft(DWr.data{n}', DWr.cal.Fs, ...
																	length(DWr.data{n}));
	
	freq(n) = fmax;
	% compute pll magnitude and phase at this frequency
	[mag(n), phi(n)] = fitsinvec(DWr.data{n}(1:ms2bin(100, DWr.cal.Fs)), ...
										1, DWr.cal.Fs, freq(n));

end


csvwrite('DWraw_vals.csv', [freq', mag']);

%%
figure

dbvals.freqs = freqs;
if length(dbvals.freqs) ~= length(dbvals.dbvals)
	error('mismatch in # freqs and dbvals');
end
plot(freqs*0.001, dbvals.dbvals, '.-')
xlabel('freqs (kHz)')
ylabel('dB SPL')
grid on
grid minor
%%

fp = fopen('LSY-A_10V_9Jul2019.txt', 'wt');

for n = 1:length(freqs)
	fprintf(fp, '%f\t%f\n', freqs(n), dbvals.dbvals(n));
end

fclose(fp);
