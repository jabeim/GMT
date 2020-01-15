function win=Leos_window(nFFT)

w=nFFT;
ftsize=nFFT;


if rem(w, 2) == 0 % force window to be odd-len
w = w + 1;
end
win = zeros(1, ftsize);
halff = ftsize/2; % midpoint of win
halflen = (w-1)/2;
acthalflen = min(halff, halflen);
halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);

