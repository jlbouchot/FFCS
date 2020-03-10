% Script to display the noise robustness of the approach when dealing with
% erasure. 

clear all
close all
clc


preMesurementNoises = 0:0.001:0.026;
nbTests = length(preMesurementNoises);

l2_fl1a = zeros(nbTests,1);
l2_l1a = zeros(nbTests,1);
l2_fcs = zeros(nbTests,1);
l2_cs = zeros(nbTests,1);

ssim_fl1a = zeros(nbTests,1);
ssim_l1a = zeros(nbTests,1);
ssim_fcs = zeros(nbTests,1);
ssim_cs = zeros(nbTests,1);

psnr_fl1a = zeros(nbTests,1);
psnr_l1a = zeros(nbTests,1);
psnr_fcs = zeros(nbTests,1);
psnr_cs = zeros(nbTests,1);

for oneTest=1:nbTests
    curNoise = preMesurementNoises(oneTest);
    
    curPath = fullfile('.','results', 'NoiseBehaviour', 'db4', ['CurrNoise', num2str(curNoise)], 'allResults.mat');
    load(curPath);
    
    l2_fl1a(oneTest) = cameraman.l2error.fl1a;
    l2_l1a(oneTest) = cameraman.l2error.l1a;
    l2_cs(oneTest) = cameraman.l2error.cs;
    l2_fcs(oneTest) = cameraman.l2error.fcs;
    
    ssim_fl1a(oneTest) = cameraman.ssim.fl1a;
    ssim_l1a(oneTest) = cameraman.ssim.l1a;
    ssim_cs(oneTest) = cameraman.ssim.cs;
    ssim_fcs(oneTest) = cameraman.ssim.fcs;
    
    psnr_fl1a(oneTest) = cameraman.psnr.fl1a;
    psnr_l1a(oneTest) = cameraman.psnr.l1a;
    psnr_cs(oneTest) = cameraman.psnr.cs;
    psnr_fcs(oneTest) = cameraman.psnr.fcs;
    
end

figure
hold on
plot(preMesurementNoises,l2_fl1a, 'r.')
plot(preMesurementNoises,l2_l1a, 'b.')
plot(preMesurementNoises,l2_fcs, 'k.')
plot(preMesurementNoises,l2_cs, 'c.')
legend('Fused L1 Analysis', 'L1 Analysis', 'Fused Compressed sensing', 'Compressed sensing')
xlabel('Percentage of dead pixels')
ylabel('L2 error of the recovery')


figure
hold on
plot(preMesurementNoises,ssim_fl1a, 'r.')
plot(preMesurementNoises,ssim_l1a, 'b.')
plot(preMesurementNoises,ssim_fcs, 'k.')
plot(preMesurementNoises,ssim_cs, 'c.')
legend('Fused L1 Analysis', 'L1 Analysis', 'Fused Compressed sensing', 'Compressed sensing')
xlabel('Percentage of dead pixels')
ylabel('SSIM score of the recovery')

figure
hold on
plot(preMesurementNoises,psnr_fl1a, 'r.')
plot(preMesurementNoises,psnr_l1a, 'b.')
plot(preMesurementNoises,psnr_fcs, 'k.')
plot(preMesurementNoises,psnr_cs, 'c.')
legend('Fused L1 Analysis', 'L1 Analysis', 'Fused Compressed sensing', 'Compressed sensing')
xlabel('Percentage of dead pixels')
ylabel('PSNR of the recovery')