function pred = model_prediction(Para_set)
% calculate the maximal ethanol production with varied intial ratio and delay time based on the parameter set

Para_Compo = Para_set([1:4,7:end]);
Para_Tempo = Para_set([1:2,5:end]);

opts = odeset('NonNegative',[1:9],'AbsTol',1e-7,'RelTol',1e-5);

% calculate the maximal ethanol production with varied intial ratio
RATIO_G = (0:0.01:1)'; % Initial fraction of G
initial_od = 10;

ini_Glu_c = 78;
ini_Xyl_c = 47;

Em_comp = zeros(length(RATIO_G),1);


for i = 1:length(RATIO_G)
    OD_G = RATIO_G(i)/1 * initial_od;
    OD_X = initial_od - OD_G;
    Y0C = [ini_Xyl_c;ini_Glu_c;0;OD_X;OD_G;0;0;1;1];
    [~,Y] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_Compo),[0 120],Y0C,opts);

    Em_comp(i) = max(Y(:,3));
end


% calculate the maximal ethanol production with varied delay time
DTset = [-50:50]; % delay time
Tnum = length(DTset);

Em_temp = zeros(Tnum,1);

ini_Glu_t = 71;
ini_Xyl_t = 41;
OD_each = 9;
for j = 1: Tnum
    Stype = sign(DTset(j));
    DelayTime = abs(DTset(j));
    
    switch Stype
        case -1 % Adding YG first and YX later
            Y0T = [ini_Xyl_t, ini_Glu_t, 0, 0, OD_each, 0, 0, 1, 1];
            [t1,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_Tempo),0:0.1:DelayTime, Y0T, opts);
            Y0d = y1(end,:);
            Y0d(4) = OD_each;    % add YX
            [~,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_Tempo),t1(end):0.1:120, Y0d, opts);
            Y = [y1;y2(2:end,:)];  
        case 0 % Adding YG and YX simultaneously
            Y0T = [ini_Xyl_t, ini_Glu_t, 0, OD_each, OD_each, 0, 0, 1, 1];
            [~,Y] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_Tempo),0:0.1:120, Y0T, opts);
        case 1 % Adding YX first and YG later
            Y0T = [ini_Xyl_t, ini_Glu_t, 0, OD_each, 0, 0, 0, 1, 1];
            [t1,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_Tempo),0:0.1:DelayTime, Y0T, opts);
            Y0d = y1(end,:);
            Y0d(5) = OD_each;    % add YG
            [~,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Para_Tempo),t1(end):0.1:120, Y0d, opts);
            Y = [y1;y2(2:end,:)];  
    end

    Em_temp(j) = max(Y(:,3));
end

pred = [Em_comp;Em_temp];
end