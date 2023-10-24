function  time_trajectory_plot(Para_set, out , data)
%%% Plot setting
color_glu = [47 174 54]/255;
color_xyl = [232 66 64]/255;
color_eth = [4 51 145]/255;
color_od = [0 0 0];

yscale_left = [60, 60, 60, 60, 60, 80, 80, 80, 80, 80, 80, 80];
yscale_right = [40, 40, 40, 40, 40, 50, 50, 50, 50, 50, 50, 50];
tscale = [120, 120, 120, 120, 120, 120, 120, 120, 30 30 30 30];
tinterval = [30, 30, 30, 30, 30, 30, 30, 30, 10, 10, 10, 10];

title_list = [
    "Ratio YG:YX = 1:9"
    "Ratio YG:YX = 3:7"
    "Ratio YG:YX = 5:5"
    "Ratio YG:YX = 7:3"
    "Ratio YG:YX = 9:1"
    "YG0YX6"
    "YG0YX12"
    "YG0YX32"
    "YX0YG0"
    "YX0YG6"
    "YX0YG12"
    "YX0YG18"
];

%% Plot
figure;

for i = 1:length(data)

    subplot(4,3,i);
    
    hold on;
    
    T = data{i}.ydata(:,1);
    ymodel = Dolmodel(data{i}, Para_set);  % least-square fitting

    predlim = out.predlims{i};  % credible envelope
    obslim  = out.obslims{i};
    
    Glu_data = data{i}.ydata(:,2); % experimental data
    Xyl_data = data{i}.ydata(:,3);
    Eth_data = data{i}.ydata(:,4);
    OD_data  = data{i}.ydata(:,5);
    yyaxis left
    plot(T, Glu_data, 'square','Color', color_glu,'MarkerFaceColor',color_glu,'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',15);
    plot(T, Xyl_data, '^',     'Color', color_xyl,'MarkerFaceColor',color_xyl,'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',12);
    yyaxis right
    plot(T, Eth_data, 'v',     'Color', color_eth,'MarkerFaceColor',color_eth,'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',12);
    yyaxis left
    plot(T, OD_data,  'o',     'Color', color_od, 'MarkerFaceColor',color_od, 'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',12);

    Glu_model = ymodel(:,2);
    Xyl_model = ymodel(:,1);
    Eth_model = ymodel(:,3);
    OD_model  = ymodel(:,4);
    yyaxis left;    
    plot(T,Glu_model,'-','Color', color_glu, 'LineWidth',2,'DisplayName','Glucose');
    fill([T;T(end:-1:1)],[predlim{2}(1,:),predlim{2}(3,end:-1:1)],...
     color_glu,'EdgeColor','None','FaceAlpha',0.4');
    fill([T;T(end:-1:1)],[obslim{2}(1,:),obslim{2}(3,end:-1:1)],...
     color_glu,'EdgeColor','None','FaceAlpha',0.2');

    plot(T,Xyl_model,'-','Color', color_xyl, 'LineWidth',2,'DisplayName','Xylose');
    fill([T;T(end:-1:1)],[predlim{1}(1,:),predlim{1}(3,end:-1:1)],...
     color_xyl,'EdgeColor','None','FaceAlpha',0.4');
    fill([T;T(end:-1:1)],[obslim{1}(1,:),obslim{1}(3,end:-1:1)],...
     color_xyl,'EdgeColor','None','FaceAlpha',0.2');

    plot(T,OD_model, '-','Color', color_od,  'LineWidth',2,'DisplayName','OD600');
    fill([T;T(end:-1:1)],[predlim{4}(1,:),predlim{4}(3,end:-1:1)],...
     color_od,'EdgeColor','None','FaceAlpha',0.4');
    fill([T;T(end:-1:1)],[obslim{4}(1,:),obslim{4}(3,end:-1:1)],...
     color_od,'EdgeColor','None','FaceAlpha',0.2');

    ylim([0 yscale_left(i)]);
    yticks(0:20:yscale_left(i));
    ylabel('Glucose, Xylose, OD600')
    yyaxis right;
    plot(T,Eth_model,'-','Color', color_eth, 'LineWidth',2,'DisplayName','Ethanol');
    fill([T;T(end:-1:1)],[predlim{3}(1,:),predlim{3}(3,end:-1:1)],...
     color_eth,'EdgeColor','None','FaceAlpha',0.3');
    fill([T;T(end:-1:1)],[obslim{3}(1,:),obslim{3}(3,end:-1:1)],...
     color_eth,'EdgeColor','None','FaceAlpha',0.2');
    ylim([0 yscale_right(i)]);
    yticks(0:10:yscale_right(i));
    ylabel('Ethanol');
    xlim([0 tscale(i)]);
    xticks([0:tinterval(i):tscale(i)]);  
    xlabel('Time(h)');
    title(title_list(i));
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    ax.LineWidth = 1.0;
    ax.TickLength = [0.020,0.025];
    box on;
%    ax.Clipping = 'off';
end


end