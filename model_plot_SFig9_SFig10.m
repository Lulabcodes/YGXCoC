function model_plot_SFig9_SFig10(Para_set, data)
% Plot Supplementary Figure 9 and 10

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


%%% Model-Experiment comparison

opts = odeset('NonNegative',[1:9],'AbsTol',1e-7,'RelTol',1e-5);
Para_compo = Para_set([1:4,7:end]);
Para_tempo = Para_set([1:2,5:end]);


c = 0;
figure;
for i = 1:length(data)    
    switch data{i}.system
        case 1 % compositional DOL
            Params = Para_compo;
        case 2 % temporal DOL
            Params = Para_tempo;
    end

    T      = data{i}.ydata(:,1);
    Y0     = data{i}.y0;

    DelayCase = data{i}.delaycase;
    if DelayCase ~= 0
        DelayTime = abs(data{i}.delaytime);
        SecondStrainOD = data{i}.secondstrain;
    end
    % Solve ODE with given parameters
    switch DelayCase
        case -1 % add YG first
            [t1,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),0:0.1:DelayTime, Y0, opts);
            Y0d = y1(end,:);
            Y0d(4) = SecondStrainOD;    % add YX
            [t2,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),DelayTime:0.1:tscale(i), Y0d, opts);
            Y = [y1;y2(2:end,:)]; 
            tm = [t1;t2(2:end)];
        case 0 % add two strains simultaneously
            [tm,Y] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),0:0.1:tscale(i), Y0, opts);
        case 1 % add YX first
            [t1,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),0:0.1:DelayTime, Y0, opts);
            Y0d = y1(end,:);
            Y0d(5) = SecondStrainOD;    % add YG
            [t2,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),DelayTime:0.1:tscale(i), Y0d, opts);
            Y = [y1;y2(2:end,:)];  
            tm = [t1;t2(2:end)];
    end
    
    Glu_data = data{i}.ydata(:,2);
    Xyl_data = data{i}.ydata(:,3);
    Eth_data = data{i}.ydata(:,4);
    OD_data  = data{i}.ydata(:,5);
       

    Glu = Y(:,1);    
    Xyl = Y(:,2);
    Eth = Y(:,3);
    OD = Y(:,4)+Y(:,5);

    c = c+1;
    if i == 6
        c = c-5;  % plot compositional and temporal cases in two figures
        figure;
    end
    subplot(3,3,c);
    hold on;

    yyaxis left
    plot(T, Glu_data, 'square','Color', color_glu,'MarkerFaceColor',color_glu,'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',15);
    plot(T, Xyl_data, '^',     'Color', color_xyl,'MarkerFaceColor',color_xyl,'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',12);
    yyaxis right
    plot(T, Eth_data, 'v',     'Color', color_eth,'MarkerFaceColor',color_eth,'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',12);
    yyaxis left
    plot(T, OD_data,  'o',     'Color', color_od, 'MarkerFaceColor',color_od, 'MarkerEdgeColor','none','LineWidth',1,'MarkerSize',12);

    yyaxis left;    
    plot(tm,Glu,'-','Color', color_glu, 'LineWidth',2,'DisplayName','Glucose');
    plot(tm,Xyl,'-','Color', color_xyl, 'LineWidth',2,'DisplayName','Xylose');
    plot(tm,OD, '-','Color', color_od,  'LineWidth',2,'DisplayName','OD600');
    ylim([0 yscale_left(i)]);
    yticks(0:20:yscale_left(i));
    ylabel('Glucose, Xylose, OD600')
    yyaxis right;
    plot(tm,Eth,'-','Color', color_eth, 'LineWidth',2,'DisplayName','Ethanol');
    ylim([0 yscale_right(i)]);
    yticks(0:10:yscale_right(i));
    ylabel('Ethanol');

  %  legend('Location','northwest');
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
    ax.Clipping = 'off';
end
