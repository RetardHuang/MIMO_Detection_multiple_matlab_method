clear
seed = 1;
rng(seed);
M = 4; % Alphabet size
len = 1e4;
Ntx = 4;

SNR_dB = 12;
SNR = [];
SER = [];
while(1)
    block = 0;
    Num = 0;
    while(Num<300&&block<1e6)
        x = randi([0,M-1],len,1);
        % Use 16-QAM modulation to produce y.
        y = qammod(x,M);
        % Transmit signal through an AWGN channel.
%         ynoisy = awgn(y,SNR_dB,'measured');
        ynoisy = awgn(y,SNR_dB-10*log10(Ntx),'measured');
        % Create scatter plot from noisy data.
    %     scatterplot(ynoisy);
        % Demodulate ynoisy to recover the message.
        z=qamdemod(ynoisy,M);
        % Check symbol error rate.
        [num,~]= symerr(x,z);
        block = block + 1;
        Num = Num + num;
        if mod(block,1e3) == 0
            fprintf('There are %d error in %d blocks when SNR is %d.\n',Num,block,SNR_dB);
        end
    end
    
    SNR = [SNR SNR_dB];
    SER = [SER Num/(block*len)];
    fprintf('The SER is %f when SNR is %d.\n\n',SER(end),SNR_dB);
    SNR_dB = SNR_dB + 2;
    if(SER(end)<5e-5)
        break;
    end
end
result.SNRs = SNR;
result.SER = SER;
save(['SISO_norm_Ntx',num2str(Ntx),'.mat'],'result');
 figure
 semilogy(SNR,SER)
 grid on