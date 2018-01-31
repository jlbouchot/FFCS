% A littile script to test the haar-wavelet-FF reconstruction

% Make sure we don't have any nnoisy stuff in the memory
clear variables;
close all;
clc;

% We'll work with the simple doppler example
load('noisdopp.mat');

% Keep the intersting things somewhere
N = length(noisdopp);
H = haarmtx(N);

noisdopp = noisdopp(:);

n = ceil(log2(N));

% Define our measurement matrix
smax = 25;
m = ceil(smax*log(N));
A = 1/sqrt(m)*randn(m,N);

noise_level = 0.01;

lowest_lvl = 3; % Everything below will be considered the absolute low frequency component of the signal
nb_stations=n-lowest_lvl+2; % don't forget the low frequency! 
y = zeros(m,nb_stations);

% % % %% L1 synthesis "classical compressed sensing" approach
% % % % Define our projection matrices and measurement vectors
% % % for one_proj=lowest_lvl:n
% % %     firstIdx = 2^(one_proj-1)+1; lastIdx = min(2^one_proj,N);
% % %     Pi{one_proj-lowest_lvl+1} = H(firstIdx:lastIdx, :)'*H(firstIdx:lastIdx,:);
% % %     Di{one_proj-lowest_lvl+1} = eye(N);
% % %     y(:,one_proj-lowest_lvl+1) = A*Pi{one_proj-lowest_lvl+1}*noisdopp + noise_level*randn(m,1);
% % % end
% % % 
% % % % Keep all the low frequencies in a single measurement vector
% % % firstIdx = 1; lastIdx = 2^(lowest_lvl-1);
% % % Pi{nb_stations} = H(firstIdx:lastIdx, :)'*H(firstIdx:lastIdx,:);
% % % Di{nb_stations} = H'; Di{nb_stations}(:,[1:firstIdx-1,lastIdx+1:end]) = 0;
% % % y(:,nb_stations) = A*Pi{nb_stations}*noisdopp + noise_level*randn(m,1);
% % % 
% % % 
% % % xhat = noisdopp + l1anaFFrecovery(A, y, noise_level*m*ones(n,1), Pi, Di, eye(N));
% % % % xhat = l1anaFFrecovery(A, y, noise_level*sqrt(m)*ones(n,1), Pi, Di, eye(N));
% % % 
% % % % Now let's check with a CS reconstruction! 
y_CS = A*noisdopp + noise_level*randn(m,1);
% % % x_CS = l1minl2data(A, y_CS, noise_level*m, eye(N)); % (A, y, noise,D)
% % % 
% % % figure
% % % subplot(3,1,1), plot(noisdopp), title('Original signal')
% % % subplot(3,1,2), plot(xhat), title('Reconstructed signal')
% % % subplot(3,1,3), plot(x_CS), title('Single sensor compressed sensing')


%% L1 analysis approach
y = zeros(m,nb_stations);

% Define our projection matrices and measurement vectors
for one_proj=lowest_lvl:n
    firstIdx = 2^(one_proj-1)+1; lastIdx = min(2^one_proj,N);
    Pi{one_proj-lowest_lvl+1} = H(firstIdx:lastIdx, :)'*H(firstIdx:lastIdx,:);
    Di{one_proj-lowest_lvl+1} = H'; Di{one_proj-lowest_lvl+1}(:,[1:firstIdx-1,lastIdx+1:end]) = 0;
    y(:,one_proj-lowest_lvl+1) = A*Pi{one_proj-lowest_lvl+1}*noisdopp + noise_level*randn(m,1);
end

% Keep all the low frequencies in a single measurement vector
firstIdx = 1; lastIdx = 2^(lowest_lvl-1);
Pi{nb_stations} = H(firstIdx:lastIdx, :)'*H(firstIdx:lastIdx,:);
Di{nb_stations} = H'; Di{nb_stations}(:,[1:firstIdx-1,lastIdx+1:end]) = 0;
y(:,nb_stations) = A*Pi{nb_stations}*noisdopp + noise_level*randn(m,1);


xhat = l1anaFFrecovery(A, y, noise_level*m*ones(n,1), Pi, Di, eye(N));
% xhat = l1anaFFrecovery(A, y, noise_level*sqrt(m)*ones(n,1), Pi, Di, eye(N));

% Now let's check with a CS reconstruction! 
% y_CS = A*noisdopp + noise_level*randn(m,1);
x_CS = l1minl2data(A, y_CS, noise_level*m, H'); % (A, y, noise,D)

figure
subplot(3,1,1), plot(noisdopp), title('Original signal')
subplot(3,1,2), plot(xhat), title('Reconstructed signal')
subplot(3,1,3), plot(x_CS), title('Single sensor compressed sensing')


%% Graphs and saves
% Some units to help automatic savings 
x0=100;
y0=100;
width=550;
height=300;
xmin = 0;
xmax = length(noisdopp);
ymin = -10;
ymax = 10;

figure
plot(noisdopp), xlabel('Timestamps', 'fontsize',18), ylabel('Intensity', 'fontsize',18)
set(gcf,'units','points','position',[x0,y0,width,height])
xlim([xmin xmax])
ylim([ymin, ymax])
savefig('originalSignal')
saveas(gcf, 'originalSignal.eps')
saveas(gcf, 'originalSignal.png')
figure
plot(xhat), xlabel('Timestamps', 'fontsize',18), ylabel('Intensity', 'fontsize',18)
set(gcf,'units','points','position',[x0,y0,width,height])
xlim([xmin xmax])
ylim([ymin, ymax])
savefig('fusionFrameRecovery')
saveas(gcf, 'fusionFrameRecovery.eps')
saveas(gcf, 'fusionFrameRecovery.png')
figure
plot(x_CS), xlabel('Timestamps', 'fontsize',18), ylabel('Intensity', 'fontsize',18)
xlim([xmin xmax])
ylim([ymin, ymax])
set(gcf,'units','points','position',[x0,y0,width,height])
savefig('l1analysisRecovery')
saveas(gcf, 'l1analysisRecovery.eps')
saveas(gcf, 'l1analysisRecovery.png')