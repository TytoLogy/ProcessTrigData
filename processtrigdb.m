% script to process triggered calibration data 

% dbvals = processTriggeredDataSPL(	...
% 										'base', <base file name for .daq files>,
% 										'path',	   'E:\Galazyuk\1Jul2019', ...
% 										'sense', 0.316, 
% 										'gain', 0, 
% 										'window', 50,
% 										'lpfc', 105000,
% 										'lporder', 3,
% 										'hpfc', 500,
% 										'hporder', 3 )

%% write freqs to text file
fstep = 500;
freqs = 1000:fstep:100000;
fp = fopen('freqs.txt', 'wt');
for n = 1:length(freqs)
	fprintf(fp, '%f\n', freqs(n));
end
fclose(fp);


%%
fstep = 500;
freqs = 1000:fstep:100000;

dbvals = processTriggeredBinDataAG('mode', 'rms');

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
