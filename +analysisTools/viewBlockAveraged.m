function viewBlockAveraged(badata, params)

    % function to plot the block averaged data in channels which meet
    % particular criteria, outlined in params.keep.


    badata=bsxfun(@minus,badata,mean(badata,2));
    badataKeep = badata(params.keep,:);
    
    
    nSamples = size(badataKeep, 2);
    timeAxis = linspace(-params.dtPre, params.dtAfter, nSamples);
    
    figure('Position',[100 100 550 780])
    subplot(2,1,1); 
    plot(timeAxis, badataKeep'); 
    set(gca,'XLimSpec', 'tight');
    xlabel('Time (0.1s)');
    ylabel('log(\phi/\phi_0)'); 
    m=max(max(abs(badataKeep)));
    subplot(2,1,2); 
    imagesc(timeAxis, 1:size(badataKeep,1), badataKeep,[-1,1].*m); 
    colorbar('Location','northoutside');
    xlabel('Time (0.1s)');
    ylabel('Measurement #')

end