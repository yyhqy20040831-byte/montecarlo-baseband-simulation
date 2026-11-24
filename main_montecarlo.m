close all; clc;

here   = fileparts(mfilename('fullpath'));
outDir = fullfile(here, 'figures');
resDir = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist(resDir, 'dir'), mkdir(resDir); end

fprintf('=========== 1/7 transmitted waveforms ===========\n');
exp1_waveforms(outDir);

fprintf('=========== 2/7 decision variables ===========\n');
exp2_decision_vars(outDir);

fprintf('=========== 3/7 matched filter vs correlator ===========\n');
exp3_receivers(outDir);

fprintf('=========== 4/7 binary BER (Monte Carlo) ===========\n');
exp4_binary_ber(outDir, resDir);

fprintf('=========== 5/7 4PAM SER/BER (Monte Carlo) ===========\n');
exp5_4pam(outDir, resDir);

fprintf('=========== 6/7 Monte Carlo sample-size study ===========\n');
exp6_mc_size(outDir, resDir);

fprintf('=========== 7/7 power spectrum ===========\n');
exp7_psd(outDir);

fprintf('\nAll done. Figures in %s, data in %s\n', outDir, resDir);


function exp1_waveforms(outDir)
    Tb = 1; E = 1; sps = 32; Ts = Tb/sps;
    bits = [0 1 1 0 1 0 0 1];
    phi  = sqrt(1/Tb)*ones(1, sps);
    phi0 = sqrt(2/Tb)*[ones(1, sps/2), zeros(1, sps/2)];
    phi1 = sqrt(2/Tb)*[zeros(1, sps/2), ones(1, sps/2)];

    txAP = kron(sqrt(E)*(1-2*bits), phi);
    txOG = kron(sqrt(E)*(1-bits), phi0) + kron(sqrt(E)*bits, phi1);
    txOO = kron(sqrt(E)*bits, phi);

    % 4PAM: pair the bits and apply Gray mapping (2 bits/symbol)
    d = sqrt(2*E/5); lev = [-3*d, -d, d, 3*d];   % Gray labels: 00, 01, 11, 10
    pairs = reshape(bits, 2, []).';
    idx4 = zeros(size(pairs,1),1);
    for k = 1:size(pairs,1)
        b = pairs(k,:);
        if isequal(b,[0 0]), idx4(k)=1; elseif isequal(b,[0 1]), idx4(k)=2;
        elseif isequal(b,[1 1]), idx4(k)=3; else, idx4(k)=4; end
    end
    tx4 = kron(reshape(lev(idx4), 1, []), sqrt(1/(2*Tb))*ones(1, 2*sps));

    tBin = (0:numel(txAP)-1)*Ts;
    t4   = (0:numel(tx4)-1)*Ts;

    fig = figure('Visible','on','Color','w','Theme','light','Position',[80 80 900 640]);
    subplot(2,2,1); plot(tBin, txAP, 'LineWidth',1.3); grid on; ylim([-1.6 1.6]);
    title('Antipodal (2PAM)'); xlabel('t/T_b'); ylabel('s(t)');
    subplot(2,2,2); plot(tBin, txOG, 'LineWidth',1.3); grid on; ylim([-0.3 2.0]);
    title('Orthogonal'); xlabel('t/T_b'); ylabel('s(t)');
    subplot(2,2,3); plot(tBin, txOO, 'LineWidth',1.3); grid on; ylim([-0.3 1.6]);
    title('On-off'); xlabel('t/T_b'); ylabel('s(t)');
    subplot(2,2,4); plot(t4, tx4, 'LineWidth',1.3); grid on;
    title('4PAM (Gray, same bits)'); xlabel('t/T_b'); ylabel('s(t)');
    exportgraphics(fig, fullfile(outDir, 'fig1_waveforms.png'), ...
        'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('fig1_waveforms.png saved\n');
end


function exp2_decision_vars(outDir)
    Tb = 1; E = 1; sps = 8; Ts = Tb/sps; EbN0dB = 6;
    N0 = E/10^(EbN0dB/10); sigma = sqrt(N0/(2*Ts));
    Nbits = 4000;
    phi  = sqrt(1/Tb)*ones(1,sps);
    phi0 = sqrt(2/Tb)*[ones(1,sps/2), zeros(1,sps/2)];
    phi1 = sqrt(2/Tb)*[zeros(1,sps/2), ones(1,sps/2)];
    bits = randi([0 1], 1, Nbits);

    fig = figure('Visible','on','Color','w','Theme','light','Position',[60 60 900 640]);

    % Antipodal: histogram of the decision variable r
    tx = kron(sqrt(E)*(1-2*bits), phi);
    rx = tx + sigma*randn(size(tx));
    y = conv(rx, fliplr(phi))*Ts; r = y(sps:sps:end);
    subplot(2,2,1);
    histogram(r(bits==0), 40, 'FaceColor',[0 0.45 0.74]); hold on;
    histogram(r(bits==1), 40, 'FaceColor',[0.85 0.33 0.10]);
    xline(0,'--k'); grid on; title('Antipodal: decision variable r');
    xlabel('r'); ylabel('count'); legend('bit 0 sent','bit 1 sent','threshold');

    % Orthogonal: (r0, r1) scatter
    tx = kron(sqrt(E)*(1-bits), phi0) + kron(sqrt(E)*bits, phi1);
    rx = tx + sigma*randn(size(tx));
    y0 = conv(rx, fliplr(phi0))*Ts; y1 = conv(rx, fliplr(phi1))*Ts;
    r0 = y0(sps:sps:end); r1 = y1(sps:sps:end);
    subplot(2,2,2);
    plot(r0(bits==0), r1(bits==0), '.', 'Color',[0 0.45 0.74], 'MarkerSize',5); hold on;
    plot(r0(bits==1), r1(bits==1), '.', 'Color',[0.85 0.33 0.10], 'MarkerSize',5);
    plot([-1 2.5],[-1 2.5],'--k'); grid on; axis equal; xlim([-1 2.5]); ylim([-1 2.5]);
    title('Orthogonal: (r_0, r_1)'); xlabel('r_0'); ylabel('r_1');
    legend('bit 0 sent','bit 1 sent','decision boundary r_0 = r_1');

    % On-off keying: threshold at sqrt(E)/2
    tx = kron(sqrt(E)*bits, phi);
    rx = tx + sigma*randn(size(tx));
    y = conv(rx, fliplr(phi))*Ts; r = y(sps:sps:end);
    subplot(2,2,3);
    histogram(r(bits==0), 40, 'FaceColor',[0 0.45 0.74]); hold on;
    histogram(r(bits==1), 40, 'FaceColor',[0.85 0.33 0.10]);
    xline(sqrt(E)/2,'--k'); grid on; title('On-off: decision variable r');
    xlabel('r'); ylabel('count'); legend('bit 0 sent','bit 1 sent','threshold \surdE/2');

    % 4PAM: decision variable r with midpoint thresholds
    d = sqrt(2*E/5); lev = [-3*d,-d,d,3*d];
    sps4 = 16; Ts4 = 2*Tb/sps4; sig4 = sqrt(N0/(2*Ts4));
    phi4 = sqrt(1/(2*Tb))*ones(1,sps4);
    idx = randi([1 4],1,Nbits);
    tx = kron(lev(idx), phi4);
    rx = tx + sig4*randn(size(tx));
    y = conv(rx, fliplr(phi4))*Ts; r = y(sps4:sps4:end);
    subplot(2,2,4);
    histogram(r, 60, 'FaceColor',[0.49 0.18 0.56]); hold on;
    xline(-2*d,'--k'); xline(0,'--k'); xline(2*d,'--k');
    grid on; title('4PAM: decision variable r'); xlabel('r'); ylabel('count');
    legend('received samples','decision thresholds');

    exportgraphics(fig, fullfile(outDir, 'fig2_decision_vars.png'), ...
        'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('fig2_decision_vars.png saved\n');
end


function exp3_receivers(outDir)
    Tb = 1; E = 1; sps = 16; Ts = Tb/sps;
    Nbits = 12; EbN0dB = 6;
    phi = sqrt(1/Tb)*ones(1, sps);

    bits  = randi([0 1], 1, Nbits);
    N0    = E/10^(EbN0dB/10);
    sigma = sqrt(N0/(2*Ts));
    tx = kron(sqrt(E)*(1-2*bits), phi);
    rx = tx + sigma*randn(size(tx));

   
    yMF = conv(rx, fliplr(phi))*Ts;
    rMF = yMF(sps:sps:end).';

    
    R = reshape(rx, sps, Nbits).';
    rCOR = (R*phi.')*Ts;

    fprintf('max |r_MF - r_COR| = %.3e\n', max(abs(rMF - rCOR)));

    t = (0:numel(yMF)-1)*Ts;
    fig = figure('Visible','on','Color','w','Theme','light','Position',[80 80 900 420]);
    plot(t, yMF, 'Color',[0.7 0.7 0.7]); hold on;
    stem((1:Nbits)*Tb, rMF, 'o', 'Color',[0 0.45 0.74], 'LineWidth',1.2);
    stem((1:Nbits)*Tb+0.02, rCOR, 'x', 'Color',[0.85 0.33 0.10], 'LineWidth',1.2);
    yline(0, '--');
    grid on; xlabel('t/T_b'); ylabel('receiver output');
    legend('matched filter output (continuous)', 'MF sampled at t = kT_b', ...
           'correlator output', 'decision threshold 0');
    title(sprintf('Matched filter vs correlator, E_b/N_0 = %d dB (max diff %.1e)', ...
          EbN0dB, max(abs(rMF-rCOR))));
    exportgraphics(fig, fullfile(outDir, 'fig3_receivers.png'), ...
        'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('fig3_receivers.png saved\n');
end


function exp4_binary_ber(outDir, resDir)
    Tb    = 1;              % bit period (normalized)
    E     = 1;              % pulse energy of the "1" symbol (normalized)
    sps   = 8;              % samples per bit
    Ts    = Tb/sps;
    Nbits = 5e5;            % bits per Eb/N0 point (Monte Carlo sample size)
    EbN0dB = 0:1:8;
    EbN0   = 10.^(EbN0dB/10);
    qfun   = @(x) 0.5*erfc(x/sqrt(2));

    phi  = sqrt(1/Tb)*ones(1, sps);                              % whole-bit rectangle
    phi0 = sqrt(2/Tb)*[ones(1, sps/2), zeros(1, sps/2)];          % first half-bit
    phi1 = sqrt(2/Tb)*[zeros(1, sps/2), ones(1, sps/2)];          % second half-bit
    fprintf('Basis energies: phi=%.3f phi0=%.3f phi1=%.3f, <phi0,phi1>=%.3f\n', ...
            sum(phi.^2)*Ts, sum(phi0.^2)*Ts, sum(phi1.^2)*Ts, sum(phi0.*phi1)*Ts);

    berAP = zeros(size(EbN0dB));
    berOG = zeros(size(EbN0dB));
    berOO = zeros(size(EbN0dB));

    for ii = 1:numel(EbN0dB)
        N0      = E/EbN0(ii);
        sigma_n = sqrt(N0/(2*Ts));          % discrete AWGN std (PSD = N0/2)

        bits = randi([0 1], 1, Nbits);      % random binary source

        % ---------------- Antipodal ----------------
        tx = kron(sqrt(E)*(1-2*bits), phi);                 % bit0 -> +sqrt(E)
        rx = tx + sigma_n*randn(size(tx));                  % AWGN channel
        y  = conv(rx, fliplr(phi))*Ts;                      % matched filter
        r  = y(sps:sps:end);                                % sample at t = kTb
        bhat = r < 0;                                       % threshold 0
        berAP(ii) = mean(bhat ~= bits);

        % ---------------- Orthogonal ------------------
        tx = kron(sqrt(E)*(1-bits), phi0) + kron(sqrt(E)*bits, phi1);
        rx = tx + sigma_n*randn(size(tx));
        y0 = conv(rx, fliplr(phi0))*Ts;
        y1 = conv(rx, fliplr(phi1))*Ts;
        r0 = y0(sps:sps:end);
        r1 = y1(sps:sps:end);
        bhat = r1 > r0;                                     % compare two branches
        berOG(ii) = mean(bhat ~= bits);

        % ---------------- On-off keying --------------
        tx = kron(sqrt(E)*bits, phi);
        rx = tx + sigma_n*randn(size(tx));
        y  = conv(rx, fliplr(phi))*Ts;
        r  = y(sps:sps:end);
        bhat = r > sqrt(E)/2;                               % optimal threshold
        berOO(ii) = mean(bhat ~= bits);

        fprintf('Eb/N0 = %2d dB | AP %.3e | OG %.3e | OO %.3e\n', ...
                EbN0dB(ii), berAP(ii), berOG(ii), berOO(ii));
    end

    berAPth = qfun(sqrt(2*EbN0));
    berOGth = qfun(sqrt(EbN0));
    berOOth = qfun(sqrt(EbN0/2));

    fig = figure('Visible','on','Color','w','Theme','light','Position',[100 100 760 500]);
    semilogy(EbN0dB, berAPth, '-',  'Color',[0.85 0.33 0.10], 'LineWidth',1.4); hold on;
    semilogy(EbN0dB, berOGth, '-',  'Color',[0.00 0.45 0.74], 'LineWidth',1.4);
    semilogy(EbN0dB, berOOth, '-',  'Color',[0.47 0.67 0.19], 'LineWidth',1.4);
    semilogy(EbN0dB, berAP, '--o', 'Color',[0.85 0.33 0.10], 'LineWidth',1.2, ...
        'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerSize',5);
    semilogy(EbN0dB, berOG, '--s', 'Color',[0.00 0.45 0.74], 'LineWidth',1.2, ...
        'MarkerFaceColor',[0.00 0.45 0.74], 'MarkerSize',5);
    semilogy(EbN0dB, berOO, '--^', 'Color',[0.47 0.67 0.19], 'LineWidth',1.2, ...
        'MarkerFaceColor',[0.47 0.67 0.19], 'MarkerSize',5);
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('Bit error rate');
    legend('AP theory','OG theory','OO theory','AP Monte Carlo','OG Monte Carlo','OO Monte Carlo', ...
           'Location','southwest');
    title(sprintf('Binary signaling over AWGN, Monte Carlo with N = %g bits (solid: theory, dashed: simulation)', Nbits));
    exportgraphics(fig, fullfile(outDir, 'fig4_ber_binary.png'), ...
        'Resolution', 150, 'BackgroundColor', 'white');

    save(fullfile(resDir, 'ber_binary.mat'), 'EbN0dB', 'berAP', 'berOG', 'berOO', ...
         'berAPth', 'berOGth', 'berOOth', 'Nbits');
    T = table(EbN0dB(:), berAP(:), berAPth(:), berOG(:), berOGth(:), berOO(:), berOOth(:), ...
         'VariableNames', {'EbN0_dB','BER_AP_sim','BER_AP_theory','BER_OG_sim','BER_OG_theory','BER_OO_sim','BER_OO_theory'});
    writetable(T, fullfile(resDir, 'ber_binary.csv'));
    disp(T);
end


function exp5_4pam(outDir, resDir)
    Tb   = 1;                    % bit period; 4PAM symbol period = 2*Tb
    sps  = 16;                   % samples per 4PAM symbol
    Ts   = 2*Tb/sps;
    Nsym = 2.5e5;                % symbols per point (= 5e5 bits)
    EbN0dB = 0:1:10;
    EbN0  = 10.^(EbN0dB/10);
    qfun  = @(x) 0.5*erfc(x/sqrt(2));
    Eb    = 1;

    phi = sqrt(1/(2*Tb))*ones(1, sps);       % unit-energy symbol pulse
    d   = sqrt(2*Eb/5);
    lev = [-3*d, -d, d, 3*d];                % Gray: 00, 01, 11, 10

    ser = zeros(size(EbN0dB));
    ber = zeros(size(EbN0dB));

    for ii = 1:numel(EbN0dB)
        N0      = Eb/EbN0(ii);
        sigma_n = sqrt(N0/(2*Ts));

        idx  = randi([1 4], 1, Nsym);                    % symbol index 1..4
        A    = lev(idx);
        tx   = kron(A, phi);
        rx   = tx + sigma_n*randn(size(tx));
        y    = conv(rx, fliplr(phi))*Ts;
        r    = y(sps:sps:end);                           % = A + noise(N0/2)

        [~, ihat] = min(abs(r - lev.'));                 % nearest-level decision
        ser(ii) = mean(ihat ~= idx);

       
        b1map = [0 0 1 1];                               % first bit
        b2map = [0 1 1 0];                               % second bit
        bitErr = (b1map(ihat) ~= b1map(idx)) + (b2map(ihat) ~= b2map(idx));
        ber(ii) = sum(bitErr)/(2*Nsym);

        fprintf('Eb/N0 = %2d dB | 4PAM SER %.3e | BER %.3e\n', EbN0dB(ii), ser(ii), ber(ii));
    end

    serTh = 1.5*qfun(sqrt(0.8*EbN0));
    berTh = 0.75*qfun(sqrt(0.8*EbN0));
    ber2PAM = qfun(sqrt(2*EbN0));

    fig = figure('Visible','on','Color','w','Theme','light','Position',[100 100 760 500]);
    semilogy(EbN0dB, serTh, '-', 'Color',[0.49 0.18 0.56], 'LineWidth',1.4); hold on;
    semilogy(EbN0dB, berTh, '-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.4);
    semilogy(EbN0dB, ber2PAM, ':', 'Color',[0.2 0.2 0.2], 'LineWidth',1.6);
    semilogy(EbN0dB, ser, '--o', 'Color',[0.49 0.18 0.56], 'LineWidth',1.2, ...
        'MarkerFaceColor',[0.49 0.18 0.56], 'MarkerSize',5);
    semilogy(EbN0dB, ber, '--s', 'Color',[0.85 0.33 0.10], 'LineWidth',1.2, ...
        'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerSize',5);
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('Error rate');
    legend('4PAM SER theory','4PAM BER theory (Gray)','2PAM BER theory', ...
           '4PAM SER Monte Carlo','4PAM BER Monte Carlo','Location','southwest');
    title('4PAM vs 2PAM at the same bit energy E_b (solid: theory, dashed: simulation)');
    exportgraphics(fig, fullfile(outDir, 'fig5_ber_4pam.png'), ...
        'Resolution', 150, 'BackgroundColor', 'white');

    save(fullfile(resDir, 'ber_4pam.mat'), 'EbN0dB', 'ser', 'ber', 'serTh', 'berTh', 'ber2PAM');
    T = table(EbN0dB(:), ser(:), serTh(:), ber(:), berTh(:), ber2PAM(:), ...
        'VariableNames', {'EbN0_dB','SER_sim','SER_theory','BER_sim','BER_theory','BER_2PAM_theory'});
    writetable(T, fullfile(resDir, 'ber_4pam.csv'));
    disp(T);
end


function exp6_mc_size(outDir, resDir)
    E = 1; EbN0dB = 6; EbN0 = 10^(EbN0dB/10);
    N0 = E/EbN0; sigma = sqrt(N0/2);
    peTheory = 0.5*erfc(sqrt(2*EbN0)/sqrt(2));

    Nlist  = [1e3 1e4 1e5 1e6];
    Ntrial = 30;
    est = zeros(numel(Nlist), Ntrial);

    for i = 1:numel(Nlist)
        for j = 1:Ntrial
            bits = randi([0 1], 1, Nlist(i));
            r = sqrt(E)*(1-2*bits) + sigma*randn(1, Nlist(i));
            est(i,j) = mean((r < 0) ~= bits);
        end
        fprintf('N = %8g : mean %.4e, std %.4e, theory %.4e\n', ...
            Nlist(i), mean(est(i,:)), std(est(i,:)), peTheory);
    end

    fig = figure('Visible','on','Color','w','Theme','light','Position',[80 80 820 460]);
   
    semilogx(Nlist, est(:,1), 'o-', 'Color',[0.75 0.75 0.75], 'MarkerSize',4); hold on;
    for j = 2:Ntrial
        semilogx(Nlist, est(:,j), 'o-', 'Color',[0.75 0.75 0.75], ...
            'MarkerSize',4, 'HandleVisibility','off');
    end
    semilogx(Nlist, mean(est,2), 'o-', 'Color',[0 0.45 0.74], 'LineWidth',1.6, ...
        'MarkerFaceColor',[0 0.45 0.74]);
    yline(peTheory, '--', 'Color',[0.85 0.33 0.10], 'LineWidth',1.4);
    grid on; xlabel('Monte Carlo sample size N (bits)'); ylabel('estimated BER');
    legend('one Monte Carlo trial (30 grey lines)', ...
           'mean of 30 trials', ...
           'theory P_e = Q(\surd(2E_b/N_0))', 'Location','northeast');
    title(sprintf('Antipodal signaling, E_b/N_0 = %d dB: grey = 30 independent trials', EbN0dB));
    annotation('textbox', [0.60 0.53 0.36 0.19], ...
        'String', {'Grey lines = 30 independent Monte Carlo trials.', ...
                   'Each line is ONE run across the 4 sample sizes N;', ...
                   'the spread is the estimation variance, shrinking', ...
                   'as 1/sqrt(N) until it collapses onto the theory.'}, ...
        'FontSize', 9, 'Color', [0.30 0.30 0.30], ...
        'BackgroundColor', 'w', 'EdgeColor', [0.78 0.78 0.78], ...
        'Margin', 8, 'VerticalAlignment', 'top');
    exportgraphics(fig, fullfile(outDir, 'fig6_mc_size.png'), ...
        'Resolution', 150, 'BackgroundColor', 'white');

    T = table(Nlist(:), mean(est,2), std(est,0,2), repmat(peTheory,numel(Nlist),1), ...
        'VariableNames', {'N','mean_BER','std_BER','theory'});
    writetable(T, fullfile(resDir, 'mc_size.csv'));
    disp(T);
end


function exp7_psd(outDir)
    Tb = 1; sps = 32; Ts = Tb/sps; fs = 1/Ts;
    Nsym2 = 4096;                       % 2PAM symbols (= bits)
    Nsym4 = 2048;                       % 4PAM symbols (same bit count)

    b2 = randi([0 1], 1, Nsym2);
    tx2 = kron(1-2*b2, ones(1, sps));   % NRZ antipodal, unit amplitude

    d = sqrt(2/5); lev = [-3*d, -d, d, 3*d];
    i4 = randi([1 4], 1, Nsym4);
    tx4 = kron(lev(i4), ones(1, 2*sps)); % symbol period 2*Tb

    Nfft = 1024;
    [P2, f] = avgPeriodogram(tx2, fs, Nfft);
    [P4, ~] = avgPeriodogram(tx4, fs, Nfft);

    fig = figure('Visible','on','Color','w','Theme','light','Position',[80 80 820 460]);
    plot(f, 10*log10(P2/max(P2)), 'Color',[0 0.45 0.74], 'LineWidth',1.3); hold on;
    plot(f, 10*log10(P4/max(P4)), 'Color',[0.85 0.33 0.10], 'LineWidth',1.3);
    grid on; xlim([0 4]); ylim([-45 3]);
    xlabel('f / R_b'); ylabel('normalised PSD (dB)');
    legend('2PAM (symbol rate = R_b)','4PAM (symbol rate = R_b/2)','Location','northeast');
    title('Baseband PSD at the same bit rate: 4PAM halves the bandwidth');
    exportgraphics(fig, fullfile(outDir, 'fig7_psd.png'), ...
        'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('fig7_psd.png saved\n');
end


function [P, f] = avgPeriodogram(x, fs, Nfft)
    nseg = 2*floor(numel(x)/Nfft) - 1;
    w = 0.5 - 0.5*cos(2*pi*(0:Nfft-1).'/(Nfft-1));
    P = zeros(Nfft,1);
    for k = 1:nseg
        seg = x((k-1)*Nfft/2+1 : (k-1)*Nfft/2+Nfft).';
        X = fft(seg.*w, Nfft);
        P = P + abs(X).^2;
    end
    P = P/(nseg*sum(w.^2)*fs);
    f = (0:Nfft-1).'/Nfft*fs;
    P = P(1:Nfft/2+1); f = f(1:Nfft/2+1);
end
