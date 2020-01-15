function MneurToBin=NeurToBinMatrix(neuralLocsOct,nFFT,Fs)


fGrid=(0:floor(nFFT/2))*(Fs/nFFT);
fBinsOct=log2(fGrid(2:end));

% plot(fBinsOct);
% plot(1./diff(fBinsOct));%density: num of fft bins per octave
binCountPerOct=1./diff(fBinsOct);
%binCountPerOct=[binCountPerOct 2*binCountPerOct(end)-binCountPerOct(end-1)];
% stem(1:(floor(nFFT/2)-1),binCountPerOct);
%linear regression
coef=[ones(floor(nFFT/2)-1,1) (1:(floor(nFFT/2)-1))']\binCountPerOct';
% hold on
% plot(1:(floor(nFFT/2)-1),coef(1)+coef(2)*(1:(floor(nFFT/2)-1))')

scl=coef(1)+coef(2).* db2mag(neuralLocsOct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tri=@(t) heaviside(-2*abs(t)+1);
% t=-2:.01:2;
% plot(t,tri(t));
% figure;
% stem(2.^neuralLocsOct,ones(1,length(neuralLocsOct)),'Marker','none');
% hold on
% stem(2.^fBinsOct,1.2*ones(1,length(fBinsOct)),'Marker','none','color','red');
% set(gca,'XScale','log')


% tmp=[neuralLocsOct(1)-(neuralLocsOct(2)-neuralLocsOct(1)),...
%      neuralLocsOct,...
%      neuralLocsOct(end)+(neuralLocsOct(end)-neuralLocsOct(end-1))];
% kernWidth=min(...
%     tmp(2:(end-1))-tmp(1:(end-2)),...
%     tmp(3:(end  ))-tmp(2:(end-1)));    
    
nNeuralLocs=length(neuralLocsOct);
MneurToBin=zeros(floor(nFFT/2),nNeuralLocs);

I=zeros(1,floor(nFFT/2));
for k=1:nNeuralLocs
    tmp=abs(fBinsOct-neuralLocsOct(k));
    [~,I(k)]=min(tmp);
    %MneurToBin(:,k)=tmp;%.*(abs(tmp)<(10*kernWidth(k)));%tri((fBinsOct-neuralLocsOct(k))/(2*kernWidth(k)));
    MneurToBin(I(k),k)=1;%./scl(k);%/binCountPerOct(k);
end

load('preemph.mat');
[~,I]=max(emphDb);
emphDb((I+1):end)=emphDb(I);
emphDb=-emphDb;emphDb=emphDb-emphDb(1);%||||||||||||||||||||||||||||||||
scl=interp1([0 ;emphF],[0 ;emphDb],(1:floor(nFFT/2))*Fs/nFFT);

MneurToBin=bsxfun(@times,MneurToBin,db2mag(scl)'.^1);
MneurToBin(isnan(MneurToBin))=0;

%  figure;
%  subplot(1,2,1);
%  imagesc(MneurToBin);
%  subplot(1,2,2);
%  plot(I)

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------
% fBins=fBins-min(fBins);
% fBins=fBins/max(fBins);
% fBins=fBins*.8*(max(neuralLocs)-min(neuralLocs));
% fBins=fBins+min(neuralLocs);
%-------------------------------

% nNeuralLocs=length(neuralLocsOct);
% MneurToBin=zeros(floor(nFFT/2),nNeuralLocs);
% 
% % for iBin=1:floor(nFFT/2)
% %     [~,iLoc] = min(abs(fBinsOct(iBin)-neuralLocsOct));
% %     MneurToBin(iBin,iLoc) = 1;
% % end
% 
% for iLoc = 1:nNeuralLocs
%     [~, iBin] = min(abs(fBinsOct-neuralLocsOct(iLoc)));
%     MneurToBin(iBin,iLoc) = 1;
% end
% 
% binEn = sum(MneurToBin,2)';
% iNoEnergy = find(binEn == 0);
% for iBin = iNoEnergy
%     [~,iLoc] = min(abs(fBinsOct(iBin)-neuralLocsOct));
%     MneurToBin(iBin,iLoc) = 1;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % binEn = sum(MneurToBin,2)';
% % MneurToBin = MneurToBin ./ (binEn(:)*ones(1,nNeuralLocs));
% % neurEn=sum(MneurToBin);
% % MneurToBin=bsxfun(@times,MneurToBin,1./neurEn);
% 
% figure;
% imagesc(MneurToBin);
% 
% return;

function x=heaviside(x)
x(x<0)=0;