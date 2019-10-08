%------------------------------------------------------------------------
%------------------------------------------------------------------------
% processIKrig
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% script to process triggered capture data from NICal for Inga's data
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% Input Files (specific data files will change...):
%------------------------------------------------------------------------
% 
% List of freqs, atten, dblevels etc. from datawave:
%	calibrate_full_list_atten_30Higher.csv 
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
datapath = '~/Work/Data/Audio/Calibration/Inga/03Oct2019';
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
DWattfile = 'DW_WithAtten_fromDataWave_FullList_20191003_1.bin';

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
%% MIC data: find magnitudes
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[MIC.mag, MIC.phi, MIC.freq] = findMags(MIC, Freq, measure_window);

% determine unique frequency values
freq = unique(MIC.freq);
nfreq = length(freq);
att = unique(Atten);
natt = length(att); %#ok<*NASGU>
lev = unique(Level);
nlev = length(lev);

fprintf('%s\n', MIC.file);
fprintf('Frequencies (Hz) tested:\n');
fprintf('\t%d\n', freq);
fprintf('\n');
fprintf('Levels (dB SPL) tested:\n');
fprintf('\t%d\n', lev);
fprintf('\n');
fprintf('Attenuations (dB) tested:\n');
fprintf('\t%d\n', att);
fprintf('\n');

%% compute dB SPL: 
%	(1) convert MIC.mag (peak sinusoid value, in Volts) to RMS volts by
%		 multiplying by sqrt(2)/2 (RMS value of 1 cycle of a sinusoid)
%	(2) convert to Pa using scaling factor (based on cal mic sensitivity)
%	(3) use 20log10(x/2e-5 Pa) to convert to dB SPL
MIC.db = dbspl(MIC.cal.VtoPa * (sqrt(2)/2) * MIC.mag);

% write to csv file
[~, fbase] = fileparts(MIC.file);
% csvwrite([fbase '_vals.csv'], [MIC.freq, MIC.mag, MIC.db]);
% use array2table and writetable to write csv file with column labels
tmp = array2table( [MIC.freq, MIC.mag, MIC.db], 'VariableNames', ...
				{'freq', 'peak_mag', 'dbspl'});
writetable(tmp, [fbase '_vals.csv'], 'WriteVariableNames', true);

% Plot through levels
figure
lstr = cell(nlev, 1);
for l = 1:nlev
	% compute index into arrays
	startx = (1 + (l - 1)*nfreq);
	endx = (l*nfreq);
	indx = startx:endx;
	if l == 1
		plot(0.001*Freq(startx:endx), MIC.db(startx:endx), '-.');
	else
		hold on
		plot(0.001*Freq(startx:endx), MIC.db(startx:endx), '-.');
		hold off
	end
	lstr{l} = sprintf('%d dB SPL', lev(l));
end
legend(lstr)
xlabel('Frequency (kHz)')
ylabel('dB SPL');
grid on
grid minor
set(gcf, 'Name', MIC.file);
title(MIC.file, 'Interpreter', 'none');
% save as .fig and .pdf
saveas(gcf, [fbase '.fig'])
print(gcf, [fbase '.pdf'], '-dpdf')
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% DW raw data: find magnitudes
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[DWr.mag, DWr.phi, DWr.freq] = findMags(DWr, Freq, measure_window);
%write to csv file
[~, fbase] = fileparts(DWr.file);
% csvwrite([fbase '_vals.csv'], [DWr.freq, DWr.mag]);
% use array2table and writetable to write csv file with column labels
tmp = array2table( [DWr.freq, DWr.mag], 'VariableNames', ...
				{'freq', 'peak_mag'});
writetable(tmp, [fbase '_vals.csv'], 'WriteVariableNames', true);

% determine unique frequency values
freq = unique(DWr.freq);
nfreq = length(freq);
att = unique(Atten);
natt = length(att);
lev = unique(Level);
nlev = length(lev);

fprintf('%s\n', DWr.file);
fprintf('Frequencies (Hz) tested:\n');
fprintf('\t%d\n', freq);
fprintf('\n');
fprintf('Levels (dB SPL) tested:\n');
fprintf('\t%d\n', lev);
fprintf('\n');
fprintf('Attenuations (dB) tested:\n');
fprintf('\t%d\n', att);
fprintf('\n');

% Plot through levels
figure
lstr = cell(nlev, 1);
for l = 1:nlev
	% compute index into arrays
	startx = (1 + (l - 1)*nfreq);
	endx = (l*nfreq);
	indx = startx:endx;
	if l == 1
		plot(0.001*Freq(startx:endx), DWr.mag(startx:endx), '-.');
	else
		hold on
		plot(0.001*Freq(startx:endx), DWr.mag(startx:endx), '-.');
		hold off
	end
	lstr{l} = sprintf('%d dB SPL', lev(l));
end
legend(lstr)
xlabel('Frequency (kHz)')
ylabel('Peak Volts');
grid on
grid minor
set(gcf, 'Name', DWr.file);
title(DWr.file, 'Interpreter', 'none');
% save as .fig and .pdf
saveas(gcf, [fbase '.fig'])
print(gcf, [fbase '.pdf'], '-dpdf')

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% DW atten data: find magnitudes
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[DWa.mag, DWa.phi, DWa.freq] = findMags(DWa, Freq, measure_window);

% compute rms dB re: 1V
DWa.db = db( (sqrt(2)/2) * DWa.mag);

% write to csv file
[~, fbase] = fileparts(DWa.file);
% csvwrite([fbase '_vals.csv'], [DWa.freq, DWa.mag DWa.db]);
% use array2table and writetable to write csv file with column labels
tmp = array2table([DWa.freq, DWa.mag DWa.db], 'VariableNames', ...
				{'freq', 'peak_mag', 'dBV'});
writetable(tmp, [fbase '_vals.csv'], 'WriteVariableNames', true);

% determine unique frequency values
freq = unique(DWa.freq);
nfreq = length(freq);
att = unique(Atten);
natt = length(att);
lev = unique(Level);
nlev = length(lev);

fprintf('%s\n', DWa.file);
fprintf('Frequencies (Hz) tested:\n');
fprintf('\t%d\n', freq);
fprintf('\n');
fprintf('Levels (dB SPL) tested:\n');
fprintf('\t%d\n', lev);
fprintf('\n');
fprintf('Attenuations (dB) tested:\n');
fprintf('\t%d\n', att);
fprintf('\n');

% Plot through levels
figure
lstr = cell(nlev, 1);
for l = 1:nlev
	% compute index into arrays
	startx = (1 + (l - 1)*nfreq);
	endx = (l*nfreq);
	indx = startx:endx;
	if l == 1
		plot(0.001*Freq(startx:endx), DWa.mag(startx:endx), '-.');
	else
		hold on
		plot(0.001*Freq(startx:endx), DWa.mag(startx:endx), '-.');
		hold off
	end
	lstr{l} = sprintf('%d db SPL', lev(l));
end
legend(lstr)
xlabel('Frequency (kHz)')
ylabel('Peak Volts');
grid on
grid minor
set(gcf, 'Name', DWa.file);

figure
lstr = cell(nlev, 1);
for l = 1:nlev
	% compute index into arrays
	startx = (1 + (l - 1)*nfreq);
	endx = (l*nfreq);
	indx = startx:endx;
	if l == 1
		plot(0.001*Freq(startx:endx), DWa.db(startx:endx), '-.');
	else
		hold on
		plot(0.001*Freq(startx:endx), DWa.db(startx:endx), '-.');
		hold off
	end
	lstr{l} = sprintf('%d db SPL', lev(l));
end
legend(lstr)
xlabel('Frequency (kHz)')
ylabel('db (Volts RMS)');
grid on
grid minor
set(gcf, 'Name', DWa.file);
title(DWa.file, 'Interpreter', 'none');
% save as .fig and .pdf
saveas(gcf, [fbase '.fig'])
print(gcf, [fbase '.pdf'], '-dpdf')





%{
% ??????

fp = fopen('LSY-A_10V_9Jul2019.txt', 'wt');

for n = 1:length(freqs)
	fprintf(fp, '%f\t%f\n', freqs(n), dbvals.dbvals(n));
end

fclose(fp);

%}
