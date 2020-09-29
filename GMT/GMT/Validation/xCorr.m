
load('AzBio_3sent_validation.mat')

% defaultData = elData(1,:);
% temp = [ones(1,15) elData(1,1:end-15)];

% defaultData = elData(1,:);
% temp = [ones(1,15) elData(1,:)];

temp = elData(1,:);
defaultData = [ones(1,15) elData(1,:)];


if length(defaultData) == length(temp)    
    a = defaultData;
    b = temp;
elseif length(defaultData) > length(temp)
    a = defaultData;
    b = [temp zeros(1,length(defaultData)-length(temp))];    
else
    a = [defaultData zeros(1,length(temp)-length(defaultData))];
    b = temp;
end

[r,lags] = xcorr(a,b);
peakInd = lags(abs(r)==max(abs(r)));



if peakInd == 0
    c = defaultData;
    d = temp;
elseif peakInd > 0
    c = defaultData(peakInd+1:end);
    d = temp(1:length(c));
   
else   
    d = temp(abs(peakInd)+1:end);
    c = defaultData(1:length(d));
end

plot([0:length(a)-1],a,'LineWidth',2,'Color','c')
hold on
plot([0:length(b)-1],b,'LineWidth',2,'Color','m')

plot([0:length(c)-1],c,'--','LineWidth',2,'Color','y')
plot([0:length(d)-1],d,':','LineWidth',2,'Color','k')
hold off

legend({'validation data','submission','validation shift','submission shift'})

xlim([6800 7800])


outputDifferences = sum(abs(c-d));