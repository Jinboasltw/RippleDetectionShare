function y = f_smooth(signal, windowlength, overlap)
if mod(windowlength,2)==0
    windowlength = windowlength + 1;
    overlap = overlap + 1;
end
delta = windowlength - overlap;
%% AVG FILTER SETUP
coeffFilterWin = ones(1, windowlength)/windowlength;
fDelay = round((length(coeffFilterWin)-1)/2);
%% APPLY AVG FILTER
indices = 1:delta:length(signal);
% Zeropad signal
signal(end+1:indices(end)+fDelay) = 0;
y = zeros(1, length(signal));
signalWinFiltered = filter(coeffFilterWin, 1, signal);
y = signalWinFiltered(fDelay+1:end);