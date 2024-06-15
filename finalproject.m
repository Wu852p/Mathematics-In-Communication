clear all;
close all;

SNRindB1 = 0:2:15;
SNRindB2 = 0:0.1:15;
M = 64;  % 64-QAM
k = log2(M);


N = 100000; 
Es = 1;  % 每symbol能量
Eb = Es / k;

for i = 1:length(SNRindB1)
    snr_in_dB = SNRindB1(i);
    snr = 10^(snr_in_dB / 10);  % 每bit的SNR
    N0 = Eb / snr;  
    sgma = sqrt(N0 / 2);  

    % 生成隨機數據源
    dsource = randi([0, M-1], N, 1);  % 生成0到63之间的均匀随机数

    % 使用QAM調制
    qam_sig = qammod(dsource, M, 'UnitAveragePower', true);

    % 添加噪聲
    noise = sgma * (randn(N, 1) + 1i * randn(N, 1));
    r = qam_sig + noise;

    % 使用QAM解調
    demod_sig = qamdemod(r, M, 'UnitAveragePower', true);

    % 計算錯誤率
    numoferr = sum(demod_sig ~= dsource);
    smld_err_prb(i) = numoferr / N;
end

% 理論錯誤率計算
for i = 1:length(SNRindB2)
    SNR = exp(SNRindB2(i) * log(10) / 10);  % signal-to-noise ratio
    % theoretical symbol error rate
    theo_err_prb(i) = 4 * (1 - 1/sqrt(M)) * qfunc(sqrt((3 * log2(M) * SNR) / (M - 1)));
end

% 繪製結果
semilogy(SNRindB1, smld_err_prb, '*');
hold on;
semilogy(SNRindB2, theo_err_prb);
hold off;
legend('Simulated symbol error rate', 'Theoretical symbol error rate');
xlabel('E_b/N_0 in dB', 'fontsize', 16, 'fontname', 'Helvetica');
ylabel('Error Probability', 'fontsize', 16, 'fontname', 'Helvetica');
title('Performance of a 64-QAM system from Monte Carlo simulation', 'fontsize', 16, 'fontname', 'Helvetica');
fname = '64QAM.png';
print(fname, '-dpng');
