%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Digital Processing of Speech and Audio Signals Project %%%%%%%%%%%%%%%
% "Study and analysis of CIS stimulating strategy for Kannada sentences"%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code by : %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% S. Ashwin Hebbar (16EC133) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Sampath Koti(16EC138) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
wavDir = 'E:\NITK\Courseware\5th Semester\DPSAS\wav\'; % Path to folder where .wav files are stored.
p1=5; %% Audio File Number
[signal,Fs]= audioread([wavDir,sprintf('%d', p1),'.wav']);

% % TO PLOT INPUT SIGNAL
n = linspace(1,length(signal),length(signal));
figure;
plot(n/Fs, signal);

%%
num = 10; %Number of channels


order = 3;  %Bandpass Filter Order

fmin = 50;    %minimum frequency
fmax = 2500;  %mximum frequency
Fc = linspace(log10(fmin), log10(fmax), num); % Fc center frequency varying linearly in log scale
Fc = power(10, Fc);
%Fc = linspace(log10(1+fmin/1000), log10(1+fmax/1000), num); % Fc center frequency varying linearly in log scale
%Fc = 1000*(power(10, Fc)-1);
Q = 4; % Quality factor

% Custom bnd pass filter - Gave better results
% fc1=[11,78,112,181.5,276.8,422.6,984.5,1502.7,2293.7,2893];
% fc2=[60,100,155,250,400,955,1485,2200,2800,3500];
% Fc = (fc2 + fc1)/2;
b = zeros(num, 2*order + 1);
a = zeros(num, 2*order + 1);
y = zeros(num, length(signal));
yabs = zeros(num, length(signal));
ylp = zeros(num, length(signal));
ylp_rec = zeros(num, length(signal));
ylp_hil = zeros(num, length(signal));

[blp, alp]=butter(8, 200/(Fs/2));   %low pass filter 


figure;
for i = 1:num
    %%%%%%%%%%% Band Pass Filter %%%%%%%%%%
    [b(i,:),a(i,:)] = butter(order, [Fc(i)*(1 - (1/(2*Q))), Fc(i)*(1+(1/(2*Q)))]/(Fs/2));
    [h,w]= freqz(b(i,:),a(i,:));
    
    NormFactor = max(abs(h));   %normalization fctor for band pass filters
     w1 = find(w*Fs/(2*pi) > fmax*1.3);   %ploting band pass filters
   plot(w(1:w1(1))*Fs/(2*pi), (abs(h(1:w1(1)))/NormFactor));
     hold on;
    y(i,:) = filter(b(i,:), a(i,:), signal)/NormFactor;
   
    %%%%%%%%%% Rectification %%%%%%%%%
    yabs(i,:) = abs(y(i,:));
    yabs(i,:) = yabs(i,:) .* yabs(i,:);
    %%%%%%%%%% LPF %%%%%%%%%%
    ylp_rec(i,:) = filter(blp,alp,yabs(i,:));
    
    
    %%%%%%%%%%%%%%Hilbert Transform%%%%
    ylp_hil(i,:)=abs(hilbert(y(i,:)));
    
    
    
    %%%%%%%%%%%% Non Linear Compression %%%%%%%%%    
    ylp(i,:) = log10(abs(1+255*ylp_hil(i,:)))/log10(256);
    

    %
    
end

hold off;
%%
figure;
subplot(311);
plot(y(4,:));
subplot(312);
plot(ylp_rec(4,:));
subplot(313);
plot(ylp_hil(4,:));;

figure;
plot(ylp(4,:));
%%
%%%%%%%% Biphase pulse modulation %%%%%
tpulse = 0.004; % Time period of biphase pulse(in seconds)
t_pos = (tpulse/num)/2;
n_pos = fix(Fs * t_pos);
n_pulse = n_pos*2*num;

biphase = zeros(num, length(signal));
signal_final = zeros(length(signal),1);
signal_bi_final = zeros(length(signal),1);
channel = zeros(num, length(signal));
channel_bi = zeros(num, length(signal));

%%%% Biphase pulse %%%%%
for j = 1:length(signal)
    index = 1+rem(fix(j/(2*n_pos)), num);
    value = (-2)*rem(fix(j/n_pos),2)+1;
    biphase(index, j) = value;
end
figure;
subplot(411);
plot(biphase(1,1:1000));
subplot(412);
plot(biphase(2,1:1000));
subplot(413);
plot(biphase(num-1,1:1000));
subplot(414);
plot(biphase(num,1:1000));
%%
% [blp2,alp2]=butter(8, 4000/(Fs/2));

n = 1:length(signal);
mod = zeros(num, length(signal));

for j = 1:num
    channel1(j,:) = (ylp(j,:));
    channel_bi1(j,:) = ylp(j,:) .* biphase(j,:);
    n1 = n * (2*pi*Fc(j)/Fs);
    mod(j,:) = cos(n1);
    channel(j,:) = channel1(j,:) .* mod(j,:);
    channel_bi(j,:) = channel_bi1(j,:) .* mod(j,:);
    signal_final = signal_final + channel(j,:)';
    signal_bi_final = signal_bi_final + channel_bi(j,:)';
end

figure;
subplot(211);
plot(channel(4,:));
subplot(212);
plot(channel_bi(4,1:10000));
%%
signal_final = signal_final/max(signal_final);
figure;
plot(signal_final);
%%
%%% Reconstructed signal without biphase %%%
audiowrite([wavDir,'output\',sprintf('%d', p1),'output','.wav'], signal_final/max(signal_final), Fs);

signal_bi_final = signal_bi_final/max(signal_bi_final);
[blp1,alp1]=butter(10, 1000/(Fs/2));
signal_bi_final = filter(blp1, alp1, signal_bi_final);
signal_bi_final = signal_bi_final/max(signal_bi_final);
figure;
plot(signal_bi_final);
%%
%%% Reconstructed signal with biphase %%%
%%audiowrite([wavDir,'output\',sprintf('%d', p1),'outputbi','.wav'], signal_bi_final/max(signal_bi_final), Fs);