function model_plot_SFig14(obj_OAT,RelRange,Para_Name)
%%% Plot results of sensitivity analysis
figure;
for i = 1:size(obj_OAT,2)
    subplot(5,6,i)
    semilogx(RelRange, sqrt(obj_OAT(:,i))/sqrt(523),'k-','LineWidth',1.5);hold on;
    title(Para_Name(i));
    ylabel('RMSD');
    xlabel('Relative change')
    ylim([0.1 0.4]);
    xlim([0.2 5]);
    xticks([0.2 1 5]);
    xticklabels([0.2 1 5]);
    ax = gca;
    ax.LineWidth = 1.0;
    ax.TickLength = [0.040,0.050];
    box on;
end
end