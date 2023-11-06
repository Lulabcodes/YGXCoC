function model_plot_Fig6(Para_set)

% Model parameters
Para_compo = Para_set([1:4,7:end]);
Para_tempo = Para_set([1:2,5:end]);
opts = odeset('NonNegative',[1:9],'AbsTol',1e-7,'RelTol',1e-5);

figure;

%%% Plot Fig6a: Model prediction and validation of compositional DOL with
%%% varied initial YG fraction
% subplot(1,2,1); 
hold on;
title('Fig.6a');


RATIO_G = (0:0.01:1)';
initial_od = 10;

ini_Xyl = 47;
ini_Glu = 78;

Em = zeros(length(RATIO_G),1);
Tm = zeros(length(RATIO_G),1);

Pexp = [
0.1   36.8	0.296810658
0.3   45.6	0.366599598
0.4   46.9	0.368131868
0.5   47.0	0.369956744
0.7   43.3	0.356232003
0.9   41.1	0.333739837
]; % The data of validation experiments. each row shows 'YG fraction' 'Maximal ethanol production' and 'Maximal ethanol yield'

for i = 1:length(RATIO_G)
    OD_G = RATIO_G(i)/1 * initial_od;
    OD_X = initial_od - OD_G;
    
    [T,Y] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_compo),[0 150],[ini_Glu;ini_Xyl;0;OD_X;OD_G;0;0;1;1],opts);
    
    Em(i) = max(Y(:,3));
    Tm(i) = T(Y(:,3)==max(Y(:,3)));
end
% Plot model prediction and experimental validation
MaxEtOHyield_C = Em/(ini_Glu+ini_Xyl);

plot(RATIO_G, MaxEtOHyield_C,'k','LineWidth',2);
plot(Pexp(:,1),Pexp(:,3),'o','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','none');
ylim([0.1 0.4])
xlabel('Initial Y_G Fraction');
ax = gca;
ax.LineWidth = 1.0;
ax.TickLength = [0.020,0.025];
box on;

% Highlight the nearly optimal region
Nearly_Optimal = RATIO_G(Em>=max(Em)*0.98);
Optimal_Ratio = RATIO_G(Em == max(Em));
fill([Nearly_Optimal;Nearly_Optimal(end:-1:1)],[0.1*ones(length(Nearly_Optimal),1); 0.4*ones(length(Nearly_Optimal),1)],...
     'g','EdgeColor','None','FaceAlpha',0.3');

% Output information about optimal and nearly optimal conditions in command window
disp(['Yield reaches maximum ', num2str(max(Em)/(ini_Glu+ini_Xyl)), ' at YG fraction ', num2str(Optimal_Ratio)]);
disp(['Nearly optimal range is from ', num2str(min(Nearly_Optimal)),' to ', num2str(max(Nearly_Optimal)),...
    ' (yield > ', num2str(0.98*max(Em)/(ini_Glu+ini_Xyl)),')']);


%%% Plot Fig6b: Model prediction and validation of temporal DOL with varied
%%% delay time

ini_Glu = 71;
ini_Xyl = 41;


Tset = [1e-4,1:1:50]';
Tnum = length(Tset);

Emp = zeros(Tnum,1);

Pexpt = [
-32	    31.321055	0.283641539
-12	    31.17411	0.286896141
 -6	    30.135	    0.284713255
  0	    35.902391	0.344155770
  6	    44.094175	0.382663568
 12	    48.4034	    0.421744905
 18	    44.235645	0.373353690
 27	    47.87861	0.376423565
 40	    36.375	    0.314944438
];  


for i = 1:length(Tset)
    DT = Tset(i);
    [t1,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_tempo),[0 DT],[ini_Glu;ini_Xyl;0;9;0;0;0;1;1],opts);
    Y0 = y1(end,:);
    Y0(5) = 9;
    [t2,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_tempo),[t1(end), 120], Y0,opts);
    T = [t1(1:end-1);t2];
    Y = [y1(1:end-1,:);y2];

    Emp(i) = max(Y(:,3));
end

Emn = zeros(Tnum,1);

for i = 1:length(Tset)
    DT = Tset(i);    
    [t1,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_tempo),[0 DT],[ini_Glu;ini_Xyl;0;0;9;0;0;1;1],opts);
    Y0 = y1(end,:);
    Y0(4) = 9;
    [t2,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_tempo),[t1(end), 120], Y0,opts);
    T = [t1(1:end-1);t2];
    Y = [y1(1:end-1,:);y2];

    Emn(i) = max(Y(:,3));
end

% Plot model prediction and experimental validation
figure;
%subplot(1,2,2); 
hold on;
title('Fig.6b');
DelayT = [-Tset(end:-1:1);Tset];
MaxEtOHyield_T = [Emn(end:-1:1);Emp]/(ini_Glu+ini_Xyl);
plot(DelayT, MaxEtOHyield_T,'k','LineWidth',2);

plot(Pexpt(1:end-2,1),Pexpt(1:end-2,3),'^','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','b','MarkerEdgeColor','none');
plot(Pexpt(end-1:end,1),Pexpt(end-1:end,3),'o','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','none');
ylim([0.25 0.45]);
ylabel('Max yield (g/g)')
xlabel('Delay Time of Y_G(h)');

ax = gca;
ax.LineWidth = 1.0;
ax.TickLength = [0.020,0.025];
box on;

% Highlight the nearly optimal region
Nearly_Optimal = Tset(Emp>=max(Emp)*0.98);
Optimal_Time = Tset(Emp == max(Emp));
fill([Nearly_Optimal;Nearly_Optimal(end:-1:1)],[0.25*ones(length(Nearly_Optimal),1); 0.45*ones(length(Nearly_Optimal),1)],...
    'g','EdgeColor','None','FaceAlpha',0.3');

% Output information about optimal and nearly optimal conditions in command window
disp(['Yield reaches maximum ', num2str(max(Emp)/(ini_Glu+ini_Xyl)), ' at delay time ', num2str(Optimal_Time), 'h']);
disp(['Nearly optimal range is from ', num2str(min(Nearly_Optimal)),'h to ', num2str(max(Nearly_Optimal)),'h'...
    ' (yield > ', num2str(0.98*max(Emp)/(ini_Glu+ini_Xyl)),')']);