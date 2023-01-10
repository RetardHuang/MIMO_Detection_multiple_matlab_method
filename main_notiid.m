clear
seed = 1;
rng(seed);
%% system parameters
SNR_dB_start = 35;
delta_SNR = 2;
sym_num = 100;
target_SER = 2e-4;
max_block = 2000;
target_errnum = 1000;


% MIMO parameters
TxRx.Ntx = 12;
TxRx.Nrx = 20;
TxRx.detector = 'EPD';  %'Limited_EP','EPD','Optimal_EP','AMP','LMMSE','LMMSE_AMP'
TxRx.ChannelCorrelation = 'Low'; %'Low','Medium','High'
if( TxRx.Ntx > TxRx.Nrx )
     error('Antenna number of Tx should be not more than Rx for UL-MIMO.')
end
Ntx = TxRx.Ntx;
Nrx = TxRx.Nrx;
Hc = zeros(Nrx,Ntx);
TxRx.Modulation_order = 4; %2,4,6,8
RM = 2^TxRx.Modulation_order;
[TxRx.Constellations,TxRx.ConsR] = modulation(TxRx);
Constellations = TxRx.ConsR;
TxRx.Es = sum(abs(TxRx.Constellations).^2)/2^(TxRx.Modulation_order);  
Iter = 1;

kappa = 100;   %condition munber of H
d2 = kappa.^(-2/Ntx*(0:Ntx-1));
d2 = d2/sum(d2) * Nrx;
d = sqrt(d2)';

filename = [TxRx.detector,'_Kappa',num2str(kappa),'_',num2str(Ntx),'T',num2str(Nrx),'R_',num2str(2^TxRx.Modulation_order ),'QAM_',...
    num2str(Iter),'_',num2str(sym_num),'symbols_',num2str(SNR_dB_start),'dB_start',num2str(target_errnum),'errors'];
% filename = 'Debug';
SER = [];
EsErr = [];
SNR_dB = SNR_dB_start;
SNRs = [];

while 1
    SNRs(end+1) = SNR_dB;
    SNR = 10^(SNR_dB/10);
    N0 = Ntx*TxRx.Es/SNR;
%     N0 = TxRx.Es/SNR;
    block_num = 0;
    sym_err = zeros(Iter,1);
    EstimationSize = 0;
    RealSize = 10;
    while block_num < max_block
        
        [x,sym_pos ]= tx(Ntx,sym_num,TxRx.Constellations);
        noise = sqrt(0.5)*(randn(Nrx,sym_num)+1i*randn(Nrx,sym_num));      
        
        detec_pos = zeros(Ntx,Iter,sym_num);
        for ss = 1:sym_num
           [V1,~] = qr(randn(Ntx,Ntx));
           [V2,~] = qr(randn(Ntx,Ntx)); 
           [U1,~] = qr(randn(Nrx,Nrx));
           [U2,~] = qr(randn(Nrx,Nrx));
           U1 = U1(:,1:Ntx);
           U2 = U2(:,1:Ntx);   
           Hc = sqrt(0.5*Ntx/Nrx)*(U1*diag(d)*V1 + 1i*U2*diag(d)*V2); 
           Hre=real(Hc);
           Him=imag(Hc);
           H = [Hre -1*Him;Him Hre]; %Real Value Channel
           HTH = H'*H;
           yc = Hc *x(:,ss) + noise(:,ss)*sqrt(N0);
           y = [real(yc); imag(yc)];
           HTy = H'*y;
           switch TxRx.detector
               case 'Limited_EP'
                   detec_pos(:,:,ss) = Limited_EP(TxRx,N0,Iter,HTH,HTy,Constellations); 
               case 'Limited_EP2'
                   detec_pos(:,:,ss) = Limited_EP2(TxRx,N0,Iter,HTH,HTy,Constellations); 
               case 'Optimal_EP'
                   detec_pos(:,:,ss) = Optimal_EP(TxRx,N0,Iter,HTH,HTy,Constellations); 
               case 'EPD'
                   detec_pos(:,:,ss) = EPD(TxRx,N0,Iter,HTH,HTy,Constellations,sym_pos(:,ss));
               case 'SD'
                   detec_pos(:,:,ss) = SD_Studer( TxRx,Hc,yc,RM,Iter,N0);
               case 'LMMSE'
                   detec_pos(:,:,ss) = LMMSE(TxRx,N0,HTH,HTy,Iter,TxRx.Es);
               case 'LP'
                   detec_pos(:,:,ss) = LP(TxRx,N0,Iter,Hc,yc,Constellations,sym_pos(:,ss));
               case 'AMP'
                   detec_pos(:,:,ss) = AMP(TxRx,y,H,N0,Iter,Constellations,sym_pos(:,ss));
               case 'LMMSE_AMP'
                   detec_pos(:,:,ss) =  LMMSE_AMP(TxRx,N0,H,y,Iter,TxRx.Es,Constellations);
               otherwise
                   error('detector is not supported.');
           end
           a = sum(sym_pos(:,ss) ~= squeeze(detec_pos(:,Iter,ss))); 
           b = a;
        end
        err = zeros(Iter,1);
        for ii = 1:Iter
            err(ii) = sum(sum(sym_pos ~= squeeze(detec_pos(:,ii,:))));
        end
        sym_err = sym_err + err;
        block_num = block_num + 1;
        if mod(block_num,100)==0
            fprintf('There are %d smymbols error when %d blocks have been transimitted at %d dB.\n\n',sym_err(ii),block_num,SNR_dB);
        end
        if sym_err > target_errnum
            break;
        end
    end
    SER(end+1,:) = sym_err/(block_num*Ntx*sym_num);
    EsErr(end+1,:) = EstimationSize/RealSize;
    fprintf('When SNR is %fdB the SER is %f.\n\n',SNR_dB,SER(end));
    SNR_dB = SNR_dB + delta_SNR;
    result.SER = SER;
    result.estierr = EsErr;
    result.SNRs = SNRs;
    result.TxRx = TxRx;
    save([filename,'.mat'], 'result');
    if(SER(end) < target_SER ||length(SNRs)>10)
        break;
    end
end
figure,
semilogy(result.SNRs,result.SER)
grid on