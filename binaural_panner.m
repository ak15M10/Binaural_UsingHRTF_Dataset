
% Date: 17/04/2023
% By: Akash kumar (D20ME101)
% Spatial Audio assignment:02

% Modified date : 18/04/2023
% Make it interactive
% ask user for input.

% azimuthal = 90(Left); 270(Right).

clc
clear
close all
clear sound;
% Load the HRTF dataset
Hrtf01=load('ReferenceHRTF.mat');
[music01, fs01]=audioread('nightingale_near_street-23902.mp3');
[music02,fs02]=audioread('warm-memories-135965.mp3');
%E
% Sampling frequency of HRIR dataset()
h_fs=Hrtf01.sampleRate;% Sampling rate
X_s=Hrtf01.sourcePosition; % Source position
Hrir_2=Hrtf01.hrtfData;% Impulse response (left +right)


% convert stereo to mono audio.
music11=music01(:,2);
music22=music02(:,2);
% Resample the audio signal (upsample to Hrtf Fs)
music22=resample(music22,h_fs,fs02);
music11=resample(music11,h_fs,fs01);

% music22=0.5*(music02(:,2)+music02(:,1));
% sound(music11,fs01)
% sound(music22,fs02)
% clear sound;


%% Visualize thhe location of Source position
Azimuth_src=X_s(:,1);
elevation_src=X_s(:,2);
xs=cosd(elevation_src).*cosd(Azimuth_src);
ys=cosd(elevation_src).*sind(Azimuth_src);
zs=sind(elevation_src);
X_available= [xs,ys,zs] ; % position vector of all available data points


figure(1)
%
p1=plot3(xs,ys,zs,'g.');
hold on
plot3(0,0,0,'*','LineWidth',12);plot3(0,0,0,'ko','LineWidth',5);
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
grid on
title('Location of source taken for HRIR')


%% TAke input from user
% Make it user understandable
% Angular position of the source
% source type and its angular position
% example: music 01, 02,03

% Print a message " elevation range:(-30 to 80 degree) "
% azimuthal = " 0 to 358 degree "

% fprintf('Welcome to 3D audio panner');
disp('Welcome to the 3D audio panner');
disp(' ');
disp('Please use Headphone for better experience');
disp(' ');
disp('front:Azimuth =0 deg');
disp('Left :Azimuth =90 deg');
disp('Right:Azimuth =270 deg');

n_src=1; % number of sources
% Default location
theta_nth=0;
phi_nth=0;
%
Source_inform=zeros([],4);
X_picked=zeros([],3);

X_desired=zeros([]);

available_music=1:2;

while true
    src_No=input ('Press 1 for "nightingale nearstreet" or Press 2 for "warm-memories": ');
    theta_nth = input('Enter the azimuthal angle (0-360 deg): ');
    phi_nth = input('Enter the elevation angle (-30 deg to 80 deg): ');
    gain=input('Enter the source strength (0 to 1): ');


    if phi_nth >=-30 && phi_nth<=80  && gain >=0 && gain<=1 && ismember(src_No,available_music) ==true
        Source_inform(n_src,:)=[theta_nth,phi_nth,gain,src_No];
        X_desired(n_src,:)= [cosd(phi_nth).*cosd(theta_nth),...
            cosd(phi_nth).*sind(theta_nth),sind(phi_nth)];
        n_src=n_src+1;% number of sources
    else
        disp('!!!!! Error !!!!!!! Error !!!!!! Error !!!!!');
        disp('!!!!! Error !!!!!!! Error !!!!!! Error !!!!!');
        disp('!!!!! Error !!!!!!! Error !!!!!! Error !!!!!');

        disp('Please enter a valid input');
    end

    in1=input('If you want to place another source press any key else press E for exit: ','s');
    if in1=='E'
        break;
    end
end

%
figure(1)
for ns=1:n_src-1
    %     plot3(X_picked(ns,1),X_picked(ns,2),X_picked(ns,3),'msq');
    hold on
   p2= plot3(X_desired(ns,1),X_desired(ns,2),X_desired(ns,3), 'r*','DisplayName','Desired point');
    %     legend('data points','center','','desired','closest picked point');
end
%%




% % visualize the error plot
% figure(2)
% plot(abs_error,'.');
% title('Error b/w desired point and all available data points')
% % find the index position of the closest point



%% Rendering
% start wait bar
% find the nearest point on dataset to the desired location
% Step 1: find the error (simple/less computation)(absolute error)
% Step 2: find the index position of min error. peak it
Binaural1=cell(1,n_src-1); % 3D vector (3rd Dim. for different source)
for ns=1:n_src-1
    Error=(X_available-X_desired(ns,:));
    abs_error=zeros(length(Azimuth_src),1);
    for ii=1:length(Azimuth_src)
        abs_error(ii)=norm(Error(ii,:));
    end
    [min_er,indx_i]=min(abs_error);

    phi_p=X_s(indx_i,2);theta_p=X_s(indx_i,1); % p_ picked point on dataset
    X_picked(ns,:)= [cosd(phi_p).*cosd(theta_p),cosd(phi_p).*sind(theta_p),sind(phi_p)];

    %  Plot  % Plot %
    figure(1)
    hold on
    p3=plot3(X_picked(ns,1),X_picked(ns,2),X_picked(ns,3),'msq','DisplayName','Actual point');
    %         legend('data points','center','','desired','closest picked
    %         point'); % Time consuming operation

    if Source_inform(ns,4) ==1
        Signal1=Source_inform(ns,3)*music11;
    elseif Source_inform(ns,4) ==2
        Signal1=Source_inform(ns,3)*music22; % gain * signal
    end

    % Extract the HRir of the picked point
    %%need to verify for left right
    Hrir_xx_left =Hrir_2(:,indx_i,1);
    Hrir_xx_right=Hrir_2(:,indx_i,2);

    M=length(Signal1);
    N=length(Hrir_xx_left);
    % L=M+N-1;
    sig1=[Signal1;zeros(N-1,1)];
    hir_Left =[ Hrir_xx_left;zeros(M-1,1)];
    hir_Right=[Hrir_xx_right;zeros(M-1,1)];

    FfT_left =fft(hir_Left);
    Fft_right=fft(hir_Right);

    fft_signal=fft(sig1);
    Binaural_left= ifft(FfT_left.*fft_signal);
    Binaural_right=ifft(Fft_right.*fft_signal);

    Binaural_left=Binaural_left(1:M);
    Binaural_right=Binaural_right(1:M);

    Binaural=[Binaural_left,Binaural_right];
    Binaural1(:,ns)={Binaural};
end

figure(1)
legend([p1 p2 p3],{'Data points', 'Desired points','Actual points'})
% Add all the binaural signals
% zero padding to the smallest signal
% Step01: find the maximum length among all
% initialize with the zeros, then insert the respective signals
sz=zeros(1,2);
for  ns=1:n_src-1
    sz(ns,:)=size(Binaural1{ns});
end
max_len=max(max(sz));

Audio_sig=zeros(max_len,2);
ith_Sig=zeros(max_len,2);
for  ns=1:n_src-1
    ith_Sig(1:length(Binaural1{ns}),:)=cell2mat(Binaural1(:,ns));
    Audio_sig=ith_Sig+Audio_sig;
end

% % Visualize the impulse response
% figure(3)
% plot(Hrir_xx_left,'.-b'); hold on;
% plot(Hrir_xx_right,'r.-');
% title('HRIR');
% legend('Left hrir','Right hrir');
% grid on;

% do convolution without resampling the input audio


audiowrite('Binaual rendered04.wav',Audio_sig,h_fs)
sound(Audio_sig,h_fs)
% clear sound;

%% Enter a message to reduce the sound level for ear protection
