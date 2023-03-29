%%%%%%%% WCL paper  %%%%%%%%
%%%%%%%% By.cdc     %%%%%%%%
%%%%%%%% 2022.02.14 %%%%%%%%
% 1. Consider a CSGC ambient RF source, a cluster of cooperative MIMO ambient
% tags, a multi-antenna reader.
% 2. Using a linear channel model.
clear; clc
close all

%******** main program ********
%% System model parameters.
M = 2;                      % number of ambient Tags
N = 2;                      % number of receive antennas

K = 1000;                   % number of legacy symbols covered by one backscatter symbol
Tlegacy = 4e-6;             % 4us, a symbol time of legacy transmission
Tsym = K*Tlegacy;           % 4ms, a symbol time of backscatter transmission

% note*:
% 1. consider WiFi 11b/g/n/ac transmission time is 
% T = Tsym + Tgi = 3.2us + 0.8us = 4us, and switch time of ADG902 is about 10ns.
% 2. in many papers, K could be 1 to 1000 and more, up to transmission rate
% and type of ambient source.c

numSym = 1000;              % number of total symbol
Nsamp = 1000;               % 1000 samples per backscatter symbol (>= 30 and more)

t = 0:Tsym/Nsamp:Tsym*(numSym-1/Nsamp);         % transmission time sequence, length = Nsamp

Ps_dbm = 10;                % dBm, transmission power of ambient source
fc = 2.4e9;                 % Hz, carrier frequency

% Modulation type
mod = 'ook';

MonteCarlo = 1e4;
data = [];
for BER = -5:2.5:25

    ber2 = 0;    ber1 = 0;    
    ber3 = 0;    ber4 = 0;
parfor ii = 1:MonteCarlo

%% Source node
% The ambient RF source obey CSGC distribution
Ps = 1e-3*(10^(Ps_dbm/10)); % Watt
% source signal at transmit antenna, s ~ CN(0,1)
s = sqrt(1/2)*(randn(1,length(t))+1j*randn(1,length(t))); 

%% Channel model
% S --> T, assuming the distance is close for RF-EH.
% consider only path loss, Friis formula
Gs = 1;                     % antenna gain at source
Gt = 1;                     % antenna gain at tag
lambda = 3e8 / fc;          % wave length of carrier
d = 10;
d_st1 = d;                  % m, distance between S and T
d_st2 = d;                  %+(-dd+2*dd*rand);  % the distance of #2 tag
d_st3 = d;                  % the distance of #3 tag
d_st4 = d;                  % the distance of #4 tag

v1 = 2;                     % path loss exponent between source and tag
Pst1 = Ps*(Gs*Gt*lambda^2)/((4*pi)^2*d_st1^v1);       % power at #1 tag
Pst2 = Ps*(Gs*Gt*lambda^2)/((4*pi)^2*d_st2^v1);       % power at #2 tag
Pst3 = Ps*(Gs*Gt*lambda^2)/((4*pi)^2*d_st3^v1);       % power at #3 tag
Pst4 = Ps*(Gs*Gt*lambda^2)/((4*pi)^2*d_st4^v1);       % power at #4 tag

% according to "Ambient Backscatter: Wireless Communication Out of Thin
% Air", the traditional T2T communication range is 1.5-2.5ft, i.e.,
% 0.4572-0.762m. donating dd = 0.2m
% so, tag2 place between 30cm-70cm
% small-scall fading, CN(0,1); length is ?
hst1 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));   % complete i.i.d. channels
hst2 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));   
hst3 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));  
hst4 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));  


% S --> R, assuming block Raleigh fading
Gr = 1;                     % antenna gain at receiver
d_sr = d;                   % m, distance between S and R
v2 = 2;                     % path loss exponent between source and receiver
Psr = Ps*(Gs*Gr*lambda^2)/((4*pi)^2*d_sr^v2);       % power at receiver 
% small-scall fading, CN(0,1)
hsr1 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
hsr2 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
% The add to 4x4.
hsr3 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
hsr4 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));

% d=0.1
d_tra1 = 0.10;  d_tra2 = 0.13;         
d_trb1 = 0.12;  d_trb2 = 0.10;          
d_trc1 = 0.16;  d_trc2 = 0.12; 
d_trd1 = 0.21;  d_trd2 = 0.16;
% The add to 4x4.
d_tra3 = 0.16;  d_tra4 = 0.21;
d_trb3 = 0.12;  d_trb4 = 0.16;
d_trc3 = 0.10;  d_trc4 = 0.12;
d_trd3 = 0.12;  d_trd4 = 0.10;


v3 = 2;                     % path loss exponent between tags and receiver
Ptra1 = Pst1*(Gt*Gr*lambda^2)/((4*pi)^2*d_tra1^v3); % power at receiver by tags
Ptra2 = Pst1*(Gt*Gr*lambda^2)/((4*pi)^2*d_tra2^v3); % power at receiver by tags
Ptrb1 = Pst2*(Gt*Gr*lambda^2)/((4*pi)^2*d_trb1^v3); % power at receiver by tags
Ptrb2 = Pst2*(Gt*Gr*lambda^2)/((4*pi)^2*d_trb2^v3); % power at receiver by tags

Ptrc1 = Pst3*(Gt*Gr*lambda^2)/((4*pi)^2*d_trc1^v3);
Ptrc2 = Pst3*(Gt*Gr*lambda^2)/((4*pi)^2*d_trc2^v3);
Ptrd1 = Pst4*(Gt*Gr*lambda^2)/((4*pi)^2*d_trd1^v3);
Ptrd2 = Pst4*(Gt*Gr*lambda^2)/((4*pi)^2*d_trd2^v3);

% The add to 4x4.
Ptra3 = Pst1*(Gt*Gr*lambda^2)/((4*pi)^2*d_tra3^v3);
Ptra4 = Pst1*(Gt*Gr*lambda^2)/((4*pi)^2*d_tra4^v3);
Ptrb3 = Pst2*(Gt*Gr*lambda^2)/((4*pi)^2*d_trb3^v3);
Ptrb4 = Pst2*(Gt*Gr*lambda^2)/((4*pi)^2*d_trb4^v3);

Ptrc3 = Pst3*(Gt*Gr*lambda^2)/((4*pi)^2*d_trc3^v3);
Ptrc4 = Pst3*(Gt*Gr*lambda^2)/((4*pi)^2*d_trc4^v3);
Ptrd3 = Pst4*(Gt*Gr*lambda^2)/((4*pi)^2*d_trd3^v3);
Ptrd4 = Pst4*(Gt*Gr*lambda^2)/((4*pi)^2*d_trd4^v3);

% small-scall fading, CN(0,1); length equals num of backscatter symbols
htra1 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htra2 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrb1 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrb2 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));

htrc1 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrc2 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrd1 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrd2 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));

% The add to 4x4
htra3 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htra4 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrb3 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrb4 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));

htrc3 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrc4 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrd3 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));
htrd4 = sqrt(0.5)*(randn(1,numSym)+1j*randn(1,numSym));

%% BD operation
% create symbol sequence at tags
bckSym = randi([0 1], 1, numSym);       % random backscatter symbols
if strcmp(mod, 'bpsk')
    bckSym = 2*bckSym - 1;
end

alpha1 = 1.00;                            % backscatter coefficient 1, ranging from [0-1]
alpha2 = 1.00;                            % backscatter coefficient 2, ranging from [0-1]
alpha3 = 1.00;
alpha4 = 1.00;

snr_dB = BER;                            % dB, ranging from [0 20]
sigma = sqrt(Psr/(10^((snr_dB)/10)));  % ???
w1 = (sigma/sqrt(2))*(randn(1,length(t))+1j*randn(1,length(t)));
w2 = (sigma/sqrt(2))*(randn(1,length(t))+1j*randn(1,length(t)));
w3 = (sigma/sqrt(2))*(randn(1,length(t))+1j*randn(1,length(t)));
w4 = (sigma/sqrt(2))*(randn(1,length(t))+1j*randn(1,length(t)));

% Psr = 0;  % when DLI is cancelled, the ber1 will be 0, but ber2 > 0

x_2x2 = zeros([1,1200]);
x_4x2 = zeros([1,1200]);
x_4x4 = zeros([1,1200]);
x_1x1 = zeros([1,1200]);


% optimal to matrix
S = reshape(s,Nsamp,numSym).';
W1 = reshape(w1,Nsamp,numSym).';
W2 = reshape(w2,Nsamp,numSym).';
W3 = reshape(w3,Nsamp,numSym).';
W4 = reshape(w4,Nsamp,numSym).';

% channel estimation 
    p0 = 0; p1 = 1;
    if strcmp(mod, 'bpsk')
        p0 = -1;
    end
    
    N = .2*Nsamp;
    ce0 = ones([numSym, N]);
    ce1 = ones([numSym, N]);
    symL = ones([numSym, Nsamp-2*N]);

X = [p0*ce0, p1*ce1, bckSym.'.*symL];

% 1x1
g_1 = sqrt(Psr)*hsr1;
H_1x1 = sqrt(Ptra1)*alpha1*hst1.*htra1;

% signal model of 1x1
Y_1x1 = g_1.'.*S + H_1x1.'.*X.*S + W1;
% measurement
u_1x1_0 = 1/N*sum(abs(Y_1x1(:,1:N)).^2,2);    u_1x1_1 = 1/N*sum(abs(Y_1x1(:,N+1:2*N)).^2,2);
z_1x1_ = 1/(Nsamp-2*N)*sum(abs(Y_1x1(:,2*N+1:end)).^2,2);
% Maximum Likelihood Function
re_1x1_0 = sqrt(Nsamp)./u_1x1_0.*exp(-(z_1x1_-u_1x1_0).^2*(Nsamp)./u_1x1_0.^2);
re_1x1_1 = sqrt(Nsamp)./u_1x1_1.*exp(-(z_1x1_-u_1x1_1).^2*(Nsamp)./u_1x1_1.^2);
% Threhold decision
x_1x1_ = (re_1x1_0<re_1x1_1);


% 2x2
g_2 = [sqrt(Psr)*hsr1; sqrt(Psr)*hsr2].';
H_2x2 = [   sqrt(Ptra1)*alpha1*hst1.*htra1;
            sqrt(Ptrb1)*alpha2*hst2.*htrb1;
            sqrt(Ptra2)*alpha1*hst1.*htra2;
            sqrt(Ptra2)*alpha2*hst2.*htrb2  ].';

% Signal model of 2x2
Y1_2x2 = g_2(:,1).*S+(H_2x2(:,1)+H_2x2(:,2)).*X.*S+W1;
Y2_2x2 = g_2(:,2).*S+(H_2x2(:,3)+H_2x2(:,4)).*X.*S+W2;
% Measurement 
% For Y1_2x2
u1_2x2_0 = 1/N*sum(abs(Y1_2x2(:,1:N)).^2,2);    u1_2x2_1 = 1/N*sum(abs(Y1_2x2(:,N+1:2*N)).^2,2);
z1_2x2 = 1/(Nsamp-2*N)*sum(abs(Y1_2x2(:,2*N+1:end)).^2,2);
% For Y2_2x2
u2_2x2_0 = 1/N*sum(abs(Y2_2x2(:,1:N)).^2,2);    u2_2x2_1 = 1/N*sum(abs(Y2_2x2(:,N+1:2*N)).^2,2);
z2_2x2 = 1/(Nsamp-2*N)*sum(abs(Y2_2x2(:,2*N+1:end)).^2,2);
% Maximum Likelihood Function
re_2x2_0 = sqrt(Nsamp)./u1_2x2_0.*exp(-(z1_2x2-u1_2x2_0).^2*(Nsamp)./u1_2x2_0.^2)* ...
    sqrt(Nsamp)./u2_2x2_0.*exp(-(z2_2x2-u2_2x2_0).^2*(Nsamp)./u2_2x2_0.^2);
re_2x2_1 = sqrt(Nsamp)./u1_2x2_1.*exp(-(z1_2x2-u1_2x2_1).^2*(Nsamp)./u1_2x2_1.^2)* ...
    sqrt(Nsamp)./u2_2x2_1.*exp(-(z2_2x2-u2_2x2_1).^2*(Nsamp)./u2_2x2_1.^2);
% Threhold decision
x_2x2_ = (re_2x2_0<re_2x2_1);

% 4x2
g_2;
H_4x2 = [   sqrt(Ptra2)*alpha1*hst1.*htra2;
            sqrt(Ptrb2)*alpha2*hst2.*htrb2;
            sqrt(Ptrc2)*alpha3*hst3.*htrc2;
            sqrt(Ptrd2)*alpha4*hst4.*htrd2;
            sqrt(Ptra3)*alpha1*hst1.*htra3;
            sqrt(Ptrb3)*alpha2*hst2.*htrb3;
            sqrt(Ptrc3)*alpha3*hst3.*htrc3;
            sqrt(Ptrd3)*alpha4*hst4.*htrd3  ].';

% Signal model of 4x2
Y1_4x2 = g_2(:,1).*S+(H_4x2(:,1)+H_4x2(:,2)+H_4x2(:,3)+H_4x2(:,4)).*X.*S+W2;
Y2_4x2 = g_2(:,2).*S+(H_4x2(:,5)+H_4x2(:,6)+H_4x2(:,7)+H_4x2(:,8)).*X.*S+W3;
% Measurement 
% For Y1_4x2
u1_4x2_0 = 1/N*sum(abs(Y1_4x2(:,1:N)).^2,2);    u1_4x2_1 = 1/N*sum(abs(Y1_4x2(:,N+1:2*N)).^2,2);
z1_4x2 = 1/(Nsamp-2*N)*sum(abs(Y1_4x2(:,2*N+1:end)).^2,2);
% For Y2_4x2
u2_4x2_0 = 1/N*sum(abs(Y2_4x2(:,1:N)).^2,2);    u2_4x2_1 = 1/N*sum(abs(Y2_4x2(:,N+1:2*N)).^2,2);
z2_4x2 = 1/(Nsamp-2*N)*sum(abs(Y2_4x2(:,2*N+1:end)).^2,2);
% Maximum Likelihood Function
re_4x2_0 = sqrt(N*10)./u1_4x2_0.*exp(-(z1_4x2-u1_4x2_0).^2*(N*10)./u1_4x2_0.^2)* ...
    sqrt(Nsamp)./u2_4x2_0.*exp(-(z2_4x2-u2_4x2_0).^2*(Nsamp)./u2_4x2_0.^2);
re_4x2_1 = sqrt(N*10)./u1_4x2_1.*exp(-(z1_4x2-u1_4x2_1).^2*(N*10)./u1_4x2_1.^2)* ...
    sqrt(Nsamp)./u2_4x2_1.*exp(-(z2_4x2-u2_4x2_1).^2*(Nsamp)./u2_4x2_1.^2);
% Threhold decision
x_4x2_ = (re_4x2_0<re_4x2_1);

% 4x4
g_4 = [sqrt(Psr)*hsr1; sqrt(Psr)*hsr2; sqrt(Psr)*hsr3; sqrt(Psr)*hsr4].';
H_4x4 = [   sqrt(Ptra1)*alpha1*hst1.*htra1;
            sqrt(Ptrb1)*alpha2*hst2.*htrb1;
            sqrt(Ptrc1)*alpha3*hst3.*htrc1;
            sqrt(Ptrd1)*alpha4*hst4.*htrd1;

            sqrt(Ptra2)*alpha1*hst1.*htra2;
            sqrt(Ptrb2)*alpha2*hst2.*htrb2;
            sqrt(Ptrc2)*alpha3*hst3.*htrc2;
            sqrt(Ptrd2)*alpha4*hst4.*htrd2;

            sqrt(Ptra3)*alpha1*hst1.*htra3;
            sqrt(Ptrb3)*alpha2*hst2.*htrb3;
            sqrt(Ptrc3)*alpha3*hst3.*htrc3;
            sqrt(Ptrd3)*alpha4*hst4.*htrd3;

            sqrt(Ptra4)*alpha1*hst1.*htra4;
            sqrt(Ptrb4)*alpha2*hst2.*htrb4;
            sqrt(Ptrc4)*alpha3*hst3.*htrc4;
            sqrt(Ptrd4)*alpha4*hst4.*htrd4; ].';
% Signal model of 4x4
Y1_4x4 = g_4(:,1).*S + sum(H_4x4(:,1:4),2).*X.*S   + W1;
Y2_4x4 = g_4(:,2).*S + sum(H_4x4(:,5:8),2).*X.*S   + W2;
Y3_4x4 = g_4(:,3).*S + sum(H_4x4(:,9:12),2).*X.*S  + W3;
Y4_4x4 = g_4(:,4).*S + sum(H_4x4(:,13:16),2).*X.*S + W4;
% Measurement 
% For Y1_4x4
u1_4x4_0 = 1/N*sum(abs(Y1_4x4(:,1:N)).^2,2);    u1_4x4_1 = 1/N*sum(abs(Y1_4x4(:,N+1:2*N)).^2,2);
z1_4x4 = 1/(Nsamp-2*N)*sum(abs(Y1_4x4(:,2*N+1:end)).^2,2);
% For Y2_4x2
u2_4x4_0 = 1/N*sum(abs(Y2_4x4(:,1:N)).^2,2);    u2_4x4_1 = 1/N*sum(abs(Y2_4x4(:,N+1:2*N)).^2,2);
z2_4x4 = 1/(Nsamp-2*N)*sum(abs(Y2_4x4(:,2*N+1:end)).^2,2);
% For Y3_4x4
u3_4x4_0 = 1/N*sum(abs(Y3_4x4(:,1:N)).^2,2);    u3_4x4_1 = 1/N*sum(abs(Y3_4x4(:,N+1:2*N)).^2,2);
z3_4x4 = 1/(Nsamp-2*N)*sum(abs(Y3_4x4(:,2*N+1:end)).^2,2);
% For Y2_4x2
u4_4x4_0 = 1/N*sum(abs(Y4_4x4(:,1:N)).^2,2);    u4_4x4_1 = 1/N*sum(abs(Y4_4x4(:,N+1:2*N)).^2,2);
z4_4x4 = 1/(Nsamp-2*N)*sum(abs(Y4_4x4(:,2*N+1:end)).^2,2);
% Maximum Likelihood Function
re_4x4_0 = sqrt(Nsamp)./u1_4x4_0.*exp(-(z1_4x4-u1_4x4_0).^2*(Nsamp)./u1_4x4_0.^2)* ...
    sqrt(Nsamp)./u2_4x4_0.*exp(-(z2_4x4-u2_4x4_0).^2*(Nsamp)./u2_4x4_0.^2)* ...
    sqrt(Nsamp)./u3_4x4_0.*exp(-(z3_4x4-u3_4x4_0).^2*(Nsamp)./u3_4x4_0.^2)* ...
    sqrt(Nsamp)./u4_4x4_0.*exp(-(z4_4x4-u4_4x4_0).^2*(Nsamp)./u4_4x4_0.^2);
re_4x4_1 = sqrt(Nsamp)./u1_4x4_1.*exp(-(z1_4x4-u1_4x4_1).^2*(Nsamp)./u1_4x4_1.^2)* ...
    sqrt(Nsamp)./u2_4x4_1.*exp(-(z2_4x4-u2_4x4_1).^2*(Nsamp)./u2_4x4_1.^2)* ...
    sqrt(Nsamp)./u3_4x4_1.*exp(-(z3_4x4-u3_4x4_1).^2*(Nsamp)./u3_4x4_1.^2)* ...
    sqrt(Nsamp)./u4_4x4_1.*exp(-(z4_4x4-u4_4x4_1).^2*(Nsamp)./u4_4x4_1.^2);
% Threhold decision
x_4x4_ = (re_4x4_0<re_4x4_1);


if strcmp(mod, 'bpsk')
    bckSym = (bckSym + 1)/2;
end
ber1 = ber1+biterr(x_1x1_.',bckSym)/numSym;
ber2 = ber2+biterr(x_2x2_.',bckSym)/numSym;
ber3 = ber3+biterr(x_4x2_.',bckSym)/numSym;
ber4 = ber4+biterr(x_4x4_.',bckSym)/numSym;


end
    ber_2x2 = ber2 / MonteCarlo;
    ber_4x2 = ber3 / MonteCarlo;
    ber_1x1 = ber1 / MonteCarlo;
    ber_4x4 = ber4 / MonteCarlo;

[BER, ber_1x1, ber_2x2, ber_4x2, ber_4x4]
data = [data; [BER, ber_1x1, ber_2x2, ber_4x2, ber_4x4]];
end
dlmwrite("tmp.txt",data);
