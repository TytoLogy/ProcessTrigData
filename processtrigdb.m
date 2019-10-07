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
%	7 Oct 2019 (SJS): modifications for newer 3 Oct data set
%------------------------------------------------------------------------
%------------------------------------------------------------------------



%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Settings
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%---------------------------------------------------------
% file with list of test information
%---------------------------------------------------------
infofile = 'calibrate_full_list_atten_30Higher.csv';
%---------------------------------------------------------
% path to data files triggered data file (.bin)
%---------------------------------------------------------
datapath = '/Users/sshanbhag/Work/Data/Audio/Calibration/Inga/03Oct2019';
%---------------------------------------------------------
% data files
%---------------------------------------------------------
% recorded with microphone
MICfile = 'DW_WithAtten_fromMicrophone_FullList_20191003_1.bin';
%---------------------------------------------------------
% raw output from datawave (no attenuation)
%---------------------------------------------------------
DWrawfile = 'DW_RawfromDataWave_FullList_20191003_2.bin';
%---------------------------------------------------------
% attenuated
%---------------------------------------------------------
DWattfile = 'DW_Atten_FullList_20190930_1.bin';

%---------------------------------------------------------
% Settings for processing data
%---------------------------------------------------------
% filtering high pass cutoff (Hz);
HP.Fc = 3500;
% filtering low pass cutoff (Hz);
LP.Fc = 80000;
% filter order
HP.forder = 5;
LP.forder = 5;
% ramp on/off duration (ms) to eliminate transients at onset, offset
ramp_ms = 0.5;
% measurement window (ms)
% assuming a 100 ms duration stimulus with 5 ms onset/offset ramps
measure_window = [5 95];

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% read in frequency, attenuation, desired level and scaling from file
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% read in using csvread, skip first row (header)
listdata = csvread(infofile, 1, 0);
% Assign data to individual arrays
%	columns in the infofile are:
%		trial,freq,mV,atten,dB_SPL
Trial = listdata(:, 1);
Freq = listdata(:, 2);
mV = listdata(:, 3);
Atten = listdata(:, 4);
Level = listdata(:, 5);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% MIC data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%---------------------------------------------------------
%% read in acquired data for MIC (microphone)
%---------------------------------------------------------
% read in and filter data
MIC = readAndFilterTrigData(	'file', fullfile(datapath, MICfile), ...
										'filterband', [HP.Fc LP.Fc], ...
										'filterorder', HP.forder, ...
										'showData', 'y');
% check vs. csv data
fprintf('MIC data %s has %d sweeps\n', MICfile, MIC.nsweeps);
fprintf('CSV file has %d stimuli\n', length(listdata));
if MIC.nsweeps ~= length(listdata)
	warning('Mismatch between # of triggered data and # in csv file')
end

%---------------------------------------------------------
%% process MIC data
%---------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Raw (unattenuated, straight from DataWave) data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%---------------------------------------------------------
%% read in raw DW data
%---------------------------------------------------------
% read in and filter data
DWr = readAndFilterTrigData(	'file', fullfile(datapath, DWrawfile), ...
										'filterband', [HP.Fc LP.Fc], ...
										'filterorder', HP.forder, ...
										'showData', 'y');
% check vs. csv data
fprintf('DWr data %s has %d sweeps\n', DWrawfile, DWr.nsweeps);
fprintf('CSV file has %d stimuli\n', length(listdata));
if DWr.nsweeps ~= length(listdata)
	warning('Mismatch between # of triggered data and # in csv file')
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Attenuated DataWave signal data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%---------------------------------------------------------
%% read in atten DW data
%---------------------------------------------------------
% read in and filter data
DWa = readAndFilterTrigData(	'file', fullfile(datapath, DWattfile), ...
										'filterband', [HP.Fc LP.Fc], ...
										'filterorder', HP.forder, ...
										'showData', 'y');
% check vs. csv data
fprintf('DWa data %s has %d sweeps\n', DWattfile, DWr.nsweeps);
fprintf('CSV file has %d stimuli\n', length(listdata));
if DWa.nsweeps ~= length(listdata)
	warning('Mismatch between # of triggered data and # in csv file')
end



%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% computations
%------------------------------------------------------------------------
%------------------------------------------------------------------------




%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% DW raw data: find magnitudes
%------------------------------------------------------------------------
%------------------------------------------------------------------------

[DWr.mag, DWr.phi, DWr.freq] = findMags(DWr, Freq, measure_window);


%% check vs. csv data
if(DWr.nsweeps == length(listdata))
	% calculate range of bins (samples) for measurement (skip ramp on/off)
	measure_bins = ms2bin(measure_window(1), DWr.cal.Fs)	...
							: ms2bin(measure_window(2), DWr.cal.Fs);
	% allocate mag, phi arrays
	DWr.mag = zeros(DWr.nsweeps, 1);
	DWr.phi = zeros(DWr.nsweeps, 1);
	DWr.freq = Freq;
	for n = 1:DWr.nsweeps
		% compute pll magnitude and phase at this frequency
		[DWr.mag(n), DWr.phi(n)] = fitsinvec(DWr.data{n}(measure_bins), ...
										1, DWr.cal.Fs, Freq(n));
	end

	[~, fbase] = fileparts(DWrawfile);
	csvwrite([fbase '_vals.csv'], [DWr.freq, DWr.mag]);	
else
	% check vs. csv data
	fprintf('DWr data %s has %d sweeps\n', DWrawfile, DWr.nsweeps);
	fprintf('CSV file has %d stimuli\n', length(listdata));
	% auto detect tones from raw unattenuated data
	% tolerance (in Hz) for finding peak frequency in autodetect mode
	fprintf('automatically determining test frequency\n');
	measure_bins = ms2bin(measure_window(1), DWr.cal.Fs)	...
								: ms2bin(measure_window(2), DWr.cal.Fs);
	% allocate mag, phi arrays
	DWr.mag = zeros(DWr.nsweeps, 1);
	DWr.phi = zeros(DWr.nsweeps, 1);
	DWr.freq = zeros(DWr.nsweeps, 1);
	for n = 1:DWr.nsweeps
		% get spectrum of data
		[tmpfreqs, tmpmags, fmax, magmax] = daqdbfft(DWr.data{n}', DWr.cal.Fs, ...
																		length(DWr.data{n}));

		DWr.freq(n) = fmax;
		% compute pll magnitude and phase at this frequency
		[DWr.mag(n), DWr.phi(n)] = fitsinvec(DWr.data{n}(measure_bins), ...
											1, DWr.cal.Fs, DWr.freq(n));

	end

% write to csv file
[~, fbase] = fileparts(DWrawfile);
csvwrite([fbase '_vals.csv'], [DWr.freq, DWr.mag]);
%plot values
figure
plot(DWr.freq*0.001, DWr.mag, '.')
xlabel('freqs (kHz)')
ylabel('Mag (V)')
grid on
grid minor


%{
% ??????

fp = fopen('LSY-A_10V_9Jul2019.txt', 'wt');

for n = 1:length(freqs)
	fprintf(fp, '%f\t%f\n', freqs(n), dbvals.dbvals(n));
end

fclose(fp);

%}
