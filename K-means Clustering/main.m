clear 
close all
clc

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

load ps5_data.mat

%% Filtering out low frequencies
x = RealWaveform;
f0 = 30000; % sampling rate of waveform (Hz)
fstop = 250; % stop frequency (Hz)
fNyquist = f0/2; % the Nyquist limit
n = length(x);
fall = linspace(-fNyquist,fNyquist,n);
desiredresponse = ones(n,1);
desiredresponse(abs(fall)<=fstop) = 0;
xfiltered = real(ifft(fft(x).*fftshift(desiredresponse)));

%% Plotting the waveforms
t = 1/f0*(1:1:n);
figure();
subplot(2,1,1);
plot(t,x);
title('Waveform of Electrode data');
xlabel('Time (s)');
ylabel('Amplitude (uV)');
subplot(2,1,2);
plot(t,xfiltered);
hold on
plot(t,250*ones(length(xfiltered),1));
title('Waveform of Filtered Electrode data');
xlabel('Time (s)');
ylabel('Amplitude (uV)');

%% Detecting threshold crossings
data_threshold = sign(xfiltered-250);
data_cross = zeros(n,1);
for i = 1:1:n-1
    data_cross(i+1) = data_threshold(i+1)-data_threshold(i);
end
data_cross = sign(data_cross - realmin);
spikecurve = [];
[r, c] = find(data_cross == 1);

figure(2);
Fs = 30;
t = 1/Fs*(1:1:31);

for i = 1:1:length(r)
    temp = xfiltered(r(i)-10:r(i)+20);
    spikecurve = [spikecurve temp];
    plot(t,spikecurve(:,i));
    hold on
end
plot(t,250*ones(size(spikecurve,1)));
title('Windowed Neural Spikes');
xlabel('Spike Interval (ms)');
ylabel('Amplitude (uV)');

%% K-Means and Plots
[clust1,cent1,obj1] = my_kmeans(spikecurve',size(InitTwoClusters_1,2),InitTwoClusters_1);
plots(spikecurve,clust1,cent1,obj1);

[clust2,cent2,obj2] = my_kmeans(spikecurve',size(InitTwoClusters_2,2),InitTwoClusters_2);
plots(spikecurve,clust2,cent2,obj2);

[clust3,cent3,obj3] = my_kmeans(spikecurve',size(InitThreeClusters_1,2),InitThreeClusters_1);
plots(spikecurve,clust3,cent3,obj3);

[clust4,cent4,obj4] = my_kmeans(spikecurve',size(InitThreeClusters_1,2),InitThreeClusters_2);
plots(spikecurve,clust4,cent4,obj4);


