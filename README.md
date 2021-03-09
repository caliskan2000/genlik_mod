# genlik_mod

% Amaç 
% AM modulasyonu ile yayın yapan iki radyo istasyonu, kanal ve radyolarını 
% bu kanallara ayarlamış iki dinleyicinin olduğu bir senaryo için 
% kanal ortamı sadece radyo istasyonlarından gönderilen işaretlerin 
% toplandığı ideal bir durum söz konusudur. 
% Radyo A’nın modülatörü 10kHz ve Radyo B modülatörü 30kHz frekanslarına sahiptirler. 
% Örnekleme frekansı 110 kHz tercih edilmiştir. Bu senaryoyu gerçekleştiren bir MATLAB programı yazınız. 

close all
clc
clear all

% Initial values
f1=10000;  %10KHz
f2=30000;  %30KHz
fs=110000; %110KHz
ts=1/fs;
ns=[0:ts:(1-ts)];  % zaman ekseni
eks=[-54999:1:55000];  % to obtain 1x110000

%%
% Producing sound signals
num_samp=22000;   %number of samples for a sound-signal
ts1=1/num_samp;   
ns1=[0:ts1:(1-ts1)];    %total duration is 1 second

% a=cos(2*pi*ns1*(1000*randn(1,1)))+sin(2*pi*ns1*(1000*randn(1,1))); 
a=cos(2*pi*ns1*(1000*.9))+sin(2*pi*ns1*(1000*.9)); 
% randomly generated sound

% c=cos(2*pi*ns1*(1000*randn(1,1)))+sin(2*pi*ns1*(1000*randn(1,1)));
c=cos(2*pi*ns1*(1000*4))+sin(2*pi*ns1*(1000*4));
% randomly generated sound

% audiowrite('sesa.wav',a,22000);   
% audiowrite('sesc.wav',c,22000);
% 
% x=audioread('sesa.wav');
% y=audioread('sesc.wav');

a=interp(a,5);
c=interp(c,5);
sound(a,fs)
pause(2)
sound(c,fs)

% plot of x & y in time domain

figure
subplot(2,1,1)
plot(ns,a,'r');
xlabel('Hertz,Hz');
ylabel('Genlik (Volt)');
title('x(t) zaman izgesinde');
axis([0 .001 -3 3])

subplot(2,1,2)
plot(ns,c,'r');
xlabel('Hertz,Hz');
ylabel('Genlik (Volt)');
title('y(t) zaman izgesinde');
axis([0 .003 -3 3])
x=a;
y=c;
% plot of x and y in frequency domain
xf=fft(x)/length(x);
xfm=abs(xf);
xfm=fftshift(xfm);

yf=fft(y)/length(y);
yfm=abs(yf);
yfm=fftshift(yfm);

figure
subplot(2,1,1)
plot(eks,xfm,'r');
xlabel('Hertz,Hz');
ylabel('|X(F)|');
title('x(t) frekans cevabi');
axis([-6000 6000 0 0.5])


subplot(2,1,2)
plot(eks,yfm,'r');
xlabel('Hertz,Hz');
ylabel('|Y(F)|');
title('y(t) frekans cevabi');
axis([-6000 6000 0 0.5])

% Channel
m=x.*cos(2*pi*ns*f1);
n=y.*cos(2*pi*ns*f2);

vr=1;
nr=sqrt(vr)*randn(size(m));
z=n+m+nr;
% sound(z,fs)
% plot of m(t), n(t) and z(t) in frequency domain
mf=fft(m)/length(m);
mfm=abs(mf);
mfm=fftshift(mfm);
figure
subplot(3,1,1)
plot(eks,mfm,'r');
xlabel('Hertz,Hz');
ylabel('|M(F)|');
title('m(t) frekans cevabi');
axis([-55000 55000 0 0.4])

nf=fft(n)/length(n);
nfm=abs(nf);
nfm=fftshift(nfm);
subplot(3,1,2)
plot(eks,nfm,'r');
xlabel('Hertz,Hz');
ylabel('|N(F)|');
title('n(t) frekans cevabi');
axis([-55000 55000 0 0.4])

zf=fft(z)/length(z);
zfm=abs(zf);
zfm=fftshift(zfm);
subplot(3,1,3)
plot(eks,zfm,'r');
xlabel('Hertz,Hz');
ylabel('|Z(F)|');
title('z(t) frekans cevabi');
axis([-55000 55000 0 0.4])

% DEMODULATION PROCESS
k=z.*cos(2*pi*ns*f1);
L=z.*cos(2*pi*ns*f2);

kf=fft(k)/length(k);
kfm=abs(kf);
kfm=fftshift(kfm);
figure
subplot(2,1,1)
plot(eks,kfm,'r');
xlabel('Hertz,Hz');
ylabel('|K(F)|');
title('k(t) frekans cevabi');
axis([-110000 110000 0 0.5])

Lf=fft(L)/length(L);
Lfm=abs(Lf);
Lfm=fftshift(Lfm);
subplot(2,1,2)
plot(eks,Lfm,'r');
xlabel('Hertz,Hz');
ylabel('|L(F)|');
title('L(t) frekans cevabi');
axis([-110000 110000 0 0.5])

% H1 Filter Design for signal A
wg=[1200/(fs/2)];             
wd=[2000/(fs/2)];          
gddb=1;
sddb=30;
[N,Wn]=buttord(wg,wd,gddb,sddb);
[B,A]=butter(N,Wn);
[H,W] = freqz(B,A,fs/2+1);
eH=flipud(H);
H=[eH(1:55000);H];

figure
subplot(2,1,1)
plot([-55000:1:55000],abs(H))
xlabel('Hertz,Hz');
ylabel('|H1(F)|');
axis([-10000 10000 0 1.5])
title('h1(t) filtresinin frekans cevabi');


% Filtering for signal x
suz_x=filter(B,A,k);
suz_xf=fft(suz_x)/length(suz_x);
suz_xfm=abs(suz_xf);
suz_xfm=fftshift(suz_xfm);

% H2 Filter Design for signal y
wg=[7000/(fs/2)];             
wd=[9000/(fs/2)];          
gddb=1;
sddb=30;
[N,Wn]=buttord(wg,wd,gddb,sddb);
[B,A]=butter(N,Wn);

[H,W] = freqz(B,A,fs/2+1);
eH=flipud(H);
H=[eH(1:55000);H];

subplot(2,1,2)
plot([-55000:1:55000],abs(H))
xlabel('Hertz,Hz');
ylabel('|H2(F)|');
axis([-10000 10000 0 1.5])
title('h2(t) filtresinin frekans cevabi');


% % % Filtering for signal y
suz_y=filter(B,A,L);
suz_yf=fft(suz_y)/length(suz_y);
suz_yfm=abs(suz_yf);
suz_yfm=fftshift(suz_yfm);

sound(suz_x,fs)
pause(2)
sound(suz_y,fs)

% plot of x and y
figure
subplot(2,1,1)
plot(eks,suz_xfm,'r');
xlabel('Hertz,Hz');
ylabel('|X(F)|');
axis([-6000 6000 0 0.6]);
title('Hoparlördeki x(t) frekans cevabi');

subplot(2,1,2)
plot(eks,suz_yfm,'r');
xlabel('Hertz,Hz');
ylabel('|Y(F)|');
axis([-6000 6000 0 0.6]);
title('Hoparlördeki y(t) frekans cevabi');
