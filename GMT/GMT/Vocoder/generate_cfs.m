function bands= generate_cfs(lo,hi,nBands)
%GENERATE_CFS Summary of this function goes here
%   Detailed explanation goes here

density = nBands/(hz2erb(hi)-hz2erb(lo));
bands = [];
for i = 1:nBands
    bands(i) = erb2hz(hz2erb(lo)+(i-.5)/density);
end
end

function hz = erb2hz(erb)
tmp = exp((erb-43)/11.17);
hz = (0.312-14.675*tmp)/(tmp-1)*1000;
end

function erb = hz2erb(hz)
erb = 11.17*log((hz+312)/(hz+14675))+43;
end