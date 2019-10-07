function varargout = readAndFilterTrigData(varargin)

%------------------------------------------------------------------------
% define default values for D (data) struct
%------------------------------------------------------------------------
% data file and path to data file (empty will result in querying user for
% information)
D.file = '';
D.path = '';
% filtering settings
% highpass cutoff freq
D.HP.Fc = 3000;
% lowpass cutoff
D.LP.Fc = 100000;
% filter order
D.HP.order = 5;
D.LP.order = 5;
% ramp on/off to eliminate onset/offset transients, ms
D.filter_ramp_ms = 0.5;

% show data traces while processing?
showData = true;

%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
% loop through arguments
index = 1;
while index <= nargin
	% act according to tag
	switch lower(varargin{index})
		case 'file'
			% set inputfile and path
			[fpath, fname, fext] = fileparts(varargin{index+1});
			D.path = fpath;
			% need to combine file and extension
			D.file = [fname fext];
			% increment index by 2 places
			index = index + 2;

		case 'filterband'
			tmpf = varargin{index+1};
			if ~isnumeric(tmpf)
				error('%s: filter freqs must be a number', mfilename);
			elseif any(tmpf) < 0
				error('%s: cutoff freqs must greater than 0', mfilename);
			else
				D.HP.Fc  = tmpf(1);
				D.LP.Fc  = tmpf(2);
			end
			index = index + 2;

		case 'filterorder'
			tmpf = varargin{index+1};
			if ~isnumeric(tmpf)
				error('%s: filter order must be a number', mfilename);
			elseif tmpf < 0
				error('%s: filter order must greater than 0', mfilename);
			else
				D.HP.order  = tmpf;
				D.LP.order  = tmpf;
			end
			index = index + 2;
			
		case 'hpfc'
			tmpf = varargin{index+1};
			if ~isnumeric(tmpf)
				error('%s: HP cutoff freq must be a number', mfilename);
			elseif tmpf < 0
				error('%s: HP cutoff freq must greater than 0', mfilename);
			else
				D.HP.Fc  = tmpf;
			end
			index = index + 2;
			
		case 'hporder'
			tmpf = varargin{index+1};
			if ~isnumeric(tmpf)
				error('%s: HP filter order must be a number', mfilename);
			elseif tmpf < 0
				error('%s: HP filter order must greater than 0', mfilename);
			else
				D.HP.order  = tmpf;
			end
			index = index + 2;

		case 'lpfc'
			tmpf = varargin{index+1};
			if ~isnumeric(tmpf)
				error('%s: LP cutoff freq must be a number', mfilename);
			elseif tmpf < 0
				error('%s: LP cutoff freq must greater than 0', mfilename);
			else
				D.LP.Fc  = tmpf;
			end
			index = index + 2;
			
		case 'lporder'
			tmpf = varargin{index+1};
			if ~isnumeric(tmpf)
				error('%s: LP filter order must be a number', mfilename);
			elseif tmpf < 0
				error('%s: LP filter order must greater than 0', mfilename);
			else
				D.LP.order  = tmpf;
			end
			index = index + 2;
			
		case 'ramptime'
			tmpf = varargin{index+1};
			if ~isnumeric(tmpf)
				error('%s: ramp on/off time must be a number', mfilename);
			elseif tmpf < 0
				error('%s: ramp on/off time must greater than 0', mfilename);
			else
				D.filter_ramp_ms  = tmpf;
			end
			index = index + 2;
			
		case 'showdata'
			showData = checkLogicVal(varargin{index+1});
			index = index + 2;

		otherwise
			error('%s: invalid option %s', mfilename, varargin{index});
	end	
end


%------------------------------------------------------------------------
% get datafile if not provided as input
%------------------------------------------------------------------------
if nargin == 0
	[D.file, D.path] = uigetfile('*.bin', ...
											'Select triggered acquisition data file');
	if isempty(D.file)
		fprintf('Cancelled\n');
		varargout{1} = D;
		return
	end
end

%------------------------------------------------------------------------
% read in binary data
%------------------------------------------------------------------------
tmp = readBinData(fullfile(D.path, D.file));
% store fields from tmp in D...
D.cal = tmp.cal;
D.data = tmp.data;
clear tmp

%------------------------------------------------------------------------
% process data
%------------------------------------------------------------------------

% get # of sweeps
D.nsweeps = length(D.data);
fprintf('data file %s has %d sweeps\n', D.file, D.nsweeps);

% calculate voltage to Pascal conversion factor - used to calculate dB SPL
D.cal.VtoPa = 1 ./ ...
		(D.cal.Gain(1) * invdb(D.cal.MicGain(1)) * D.cal.MicSensitivity);

% filter the data
fprintf('Filtering %s data\n', D.file);
% Nyquist frequency
fnyq = D.cal.Fs/2;
% build a highpass filter for processing the data
[D.HP.fb, D.HP.fa] = butter(D.HP.order, D.HP.Fc/fnyq, 'high');
% build a lowpass filter for processing the data
[D.LP.fb, D.LP.fa] = butter(D.LP.order, D.LP.Fc/fnyq, 'low');
% build a highpass filter for processing the data
[D.HP.fb, D.HP.fa] = butter(D.HP.order, D.HP.Fc/fnyq, 'high');
% build a lowpass filter for processing the data
[D.LP.fb, D.LP.fa] = butter(D.LP.order, D.LP.Fc/fnyq, 'low');

% loop through sweeps, apply ramp and filter
if showData
	figure
	dt = 1/D.cal.Fs;
	set(gcf, 'Name', D.file);
end


for n = 1:D.nsweeps
	% filter data
	% apply short ramp and highpass filter
	tmp = filtfilt(D.HP.fb, D.HP.fa, ...
					sin2array(D.data{n}', D.filter_ramp_ms, D.cal.Fs));
	% apply lowpass filter
	tmp = filtfilt(D.LP.fb, D.LP.fa, tmp);
	% plot data
	if showData
		% plot raw data
		subplot(2, 1, 1)
		tvec = 1000 * dt * (0:(length(D.data{n}) - 1));
		plot(tvec, D.data{n});
		ylabel('Raw (V)');
		title(sprintf('Sweep %d', n));	
		% plot filtered data	
		subplot(212)
		plot(tvec, tmp);
		xlabel('Time (ms)');
		ylabel('Filtered (V)');
		drawnow
	end
	% store filtered data
	D.data{n} = tmp;
end

% assign output
varargout{1} = D;
