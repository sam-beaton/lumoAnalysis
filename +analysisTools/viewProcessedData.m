function viewProcessedData(dataIn, info, params)

    % View pre-processed data - same as example in
    % NeuroDOT_Full_Processing_Script but created to tidy up main workflow
    % using modular function
    %
    % SLB 4/4/25
    
    figure('Position',[100 100 550 780])
    subplot(3,1,1); plot(dataIn(params.keep,:)'); 
    set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
    ylabel('log(\phi/\phi_0)') 
    m=max(max(abs(dataIn(params.keep,:))));
    subplot(3,1,2); imagesc(dataIn(params.keep,:),[-1,1].*m); 
    colorbar('Location','northoutside');
    xlabel('Time (samples)');ylabel('Measurement #')
    [ftmag,ftdomain] = fft_tts(squeeze(mean(dataIn(params.keep,:),1)),info.system.framerate); % Generate average spectrum
    subplot(3,1,3); semilogx(ftdomain,ftmag);
    xlabel('Frequency (Hz)');ylabel('|X(f)|');xlim([1e-3 1])
    
    nlrGrayPlots_180818(dataIn,info); % Gray Plot with synch points

end