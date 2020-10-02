% Compute array subtraction to estimate the similarity between vectors x
% and y. Use cross correlation to time align x and y with the shorter
% vector being zero padded to match the longer. Function computes the
% sum of the absolute difference between the two vectors after attempting time
% alignment.

function outputDifferences = xCorrSimilarity(x,y)



if length(y) == length(x)
    % no resizing necessary if arrays match in size
elseif abs(length(y)-length(x)) > 0.01*max([length(x) length(y)])
    % if size difference is too great throw an error
    error('Array size difference exceeds tolerance')
elseif length(y) > length(x)
    x = [x zeros(1,length(y)-length(x))];    
else
    y = [y zeros(1,length(x)-length(y))];
end

[r,lags] = xcorr(y,x);
peakInd = lags(abs(r)==max(abs(r)));

if peakInd == 0
elseif peakInd > 0
    y = y(peakInd+1:end);
    x = x(1:length(y));   
else   
    x = x(abs(peakInd)+1:end);
    y = y(1:length(x));
end

outputDifferences = sum(abs(x-y));