close all
load("AMP_Iter2BER.mat")
load("AMP_Iter2SNRs.mat")
semilogy(shitSNR,shitSER);hold on

load("AMP_Iter4BER.mat")
load("AMP_Iter4SNRs.mat")
semilogy(shitSNR,shitSER);hold on

load("AMP_Iter6BER.mat")
load("AMP_Iter6SNRs.mat")
semilogy(shitSNR,shitSER);hold on

load("AMP_Iter8BER.mat")
load("AMP_Iter8SNRs.mat")
semilogy(shitSNR,shitSER);hold on

load("AMP_Iter10BER.mat")
load("AMP_Iter10SNRs.mat")
semilogy(shitSNR,shitSER);hold on

title("AMP的SNR-BER曲线")
xlabel("SNR/dB")
ylabel("BER")

legend('4次迭代','2次迭代','6次迭代','8次迭代','10次迭代');    