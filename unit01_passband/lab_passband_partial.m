%% Lab:  Upconversion and downconversion

% In this lab, we will demonstrate how to:

% * Generate a random baseband signal from symbols

% * Simulate upconversion and downconversion

% * Simulate a simple passband filter

% * Measure and plot the PSD of signals 

% 

% NYU students:  Complete all the sections labeled |%TODO|.  Publish the 

% file, print to PDF and submit the PDF in NYU Classes. Do not submit the 

% MATLAB files.



%% Create a sequence of symbolss

% To illustrate the process of up- and down-conversion, we will create a

% simple sequence of complex symbols.  This type of sequence is used

% commonly digital transmissions where each symbol encodes bits to be

% transmitted.  We will discuss this later.  

%

% Right now, generate a vector |xsym| of |nsym=1024| 

% random symbols with values $(\pm 1 \pm 1)/\sqrt{2}$ complex symbol values (i.e.

% QPSK).   You may use the command |randsample|. 

nsym = 1024;



% TODO

%    sym = ...
x = zeros(nsym);
x = randsample(2,nsym,1);
x = x - 1;
x = x';
I = sqrt(-1);
sym = (x+I*x)/sqrt(2);





%% Create the baseband signal

% Let $x(t)$ be the tranmitted baseband signal $x(t)$ generated from 

% the sequence of complex symbols |xsym|.  Assume the symbol rate is

% |fysm=20| MHz.  Create a vector |x| from $x(t)$ sampled at |nov=16| 

% times the symbol rate.  Also, create a corresponding time vector |t|.

fsym = 20;  % Frequency in MHz

nov = 16;   % Over-sampling rate 

fsamp = nov*fsym;  % Sampling rate f or x



% TODO

%     x = ...

%     t = ...

x = zeros([1,1024*16]);
t = zeros([1,1024*16]);
for i=1:1024
    for j=1:16
        x(:,16*(i-1)+j)=sym(:,i);
        t(:,16*(i-1)+j)=(16*(i-1)+j-1)/fsamp;
    end
end




% TODO.

% Plot |x(t)| vs. |t| for the first 8 symbols.  Use separate plots

% for the real and imaginary components.  Label the x-axis with the 

% correct units.

x_real = real(x);
x_imag = imag(x);
plot(t,x_real);
xlabel('time/us')
xlim([0,t(1,8*16)]);
plot(t,x_imag);
xlabel('time/us')
xlim([0,t(:,8*16)]);

%% Plot the PSD of the TX signal

% Properly measurig the PSD is somewhat tricky.  Fortunately,

% Matlab has an excellent routines for computing the PSD:

%

%   [Px, fx] = pwelch(x,hamming(512),[],[],fsamp*1e6,'centered');

%

% This routine uses the Welch algorithm.  Use this routine to compute and

% plot the PSD.  Label your axes correctly.  Plot the PSD in dBm/Hz.



% TODO

[Px, fx] = pwelch(x,hamming(512),[],[],fsamp*1e6,'centered');

plot(fx,10*log10(Px))

%% Upconvert

% Create a real passband signal |xp| by upconverting |x| with a carrier

% frequency of |fc=80| MHz.  Plot the PSD of |x| and |xp| on the same plot,

% so you can compare the two.

fc = 80;


% TODO

%     xp = ...

%     [Pp, fp] = ...

xp = x_real.*cos(2*pi*fc.*t)-x_imag.*sin(2*pi*fc.*t);

[Pp, fp] = pwelch(xp,hamming(512),[],[],fsamp*1e6,'centered');

hold on

plot(fp,10*log10(Pp),'r');

legend('PSD of x','PSD of xp')

hold off

Comment
% The baseband signal is transformed to passband signal



%% Passband channel

% We now simulate a passband channel.  Suppose the passband channel is

% described by the linear system:

%

%    dyp/dt = f0*(xp-yp)

%

% for |f0=25| MHz.  We can approximate this in discrete-time by:

% 

%    yp(t+1) = yp(t) + f0*tsamp*(xp(t)-yp(t)).

% 

% Simulate the discrete-time filter with the |filter| command to 

% create a sampled version of the output |yp|.  



% TODO 

%     yp = filter(...)

f0 = 25;

b = 25/320;
a = [1,-1+25/320];
yp = filter(b,a,xp);



% Plot the PSD of |yp|. Also plot the expected PSD of |yp| using

%

%    S_{y_p}(f) = S_{x_p}(f)|G_p(f)|^2,

%

% where $G_p(f)$ is the passband frequency response.  Note that you
% can 

% compute the frequency response using the |freqs| command.  You may notice

% a small discrepancy between the measured and expected frequency response.

% This is due to the digital implementation of the filter.



% TODO:

%     Gp = ...

[Py, fy] = pwelch(yp,hamming(512),[],[],fsamp*1e6,'centered');

Gp = freqs(b,a,fp);

Spy = Pp .* (Gp .* conj(Gp));

plot(fy,10*log(Py));
hold on
plot(fp,10*log(Spy),'r');
legend('PSD of yp','PSD of Syp(f)','location','east')
hold off


%% Design a digital downconversion filter

% We will design a simple digital filter to filter the downconverted

% signal.  Use the <matlab:doc('cheby1') publish> function to create a

% digital fourth order filter with  |f|<= 12.5 MHz.  Use a passband ripple

% of 0.5 dB.  



% TODO

[blpf,alpf] = cheby1(4,0.5,12.5/320,'low');


%% Plot the filter frequency response

% Use the <matlab:doc('freqz') publish> to plot the frequency response of

% the digital filter.  Label your axes 
% TODO

[Hcheb,f] = freqz(blpf,alpf);
freqz(blpf,alpf);

%% Downconvert the signal digital domain

% Downconvert the signal |yp| by mixing and filtering it with the filter.

% Plot the PSD of the received signal.



% TODO

yi = 2*yp.*cos(2*pi*fc*t);
yq = -2*yp.*sin(2*pi*fc*t);
y = yi + I*yq;
y = filter(blpf,alpf,y);
ymm=abs(y);
[Pyy,fyy] = pwelch(y,hamming(512),[],[],fsamp*1e6,'centered');
plot(fyy,10*log(Pyy));
hold on
plot(fx,10*log(Px),'r');
legend('PSD of y','PSD of x','location','northeast')
hold off

Comment
% After the downconversion, the PSD of y is mostly in the baseband. 
