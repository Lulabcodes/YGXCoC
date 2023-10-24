function  ymod = Dolmodel(data, Para_set)

opts = odeset('NonNegative',[1:9],'AbsTol',1e-7,'RelTol',1e-5);

% initialize objective function

% Parameters for compostional DOL or temporal DOL
switch data.system
    case 1 % compositional DOL
        Params = Para_set([1:4,7:end]);
    case 2
        Params = Para_set([1:2,5:end]);
end

% Sum of squares error

T       = data.ydata(:,1);
Y0     = data.y0;
% Solve ODE with given parameters
DelayCase = data.delaycase;
if DelayCase ~= 0
    DelayTime = abs(data.delaytime);
    SecondStrainOD = data.secondstrain;
end
switch DelayCase
    case -1
        [~,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),T(T<=DelayTime), Y0, opts);
        Y0d = y1(end,:);
        Y0d(4) = SecondStrainOD;    % add YX
        [~,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),T(T>=DelayTime), Y0d, opts);
        Y = [y1;y2(2:end,:)];         
    case 0
        [~,Y] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),T, Y0, opts);
    case 1
        [~,y1] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),T(T<=DelayTime), Y0, opts);
        Y0d = y1(end,:);
        Y0d(5) = SecondStrainOD;    % add YG
        [~,y2] = ode15s(@(t,y)Kinetic_Equation(t,y,Params),T(T>=DelayTime), Y0d, opts);
        Y = [y1;y2(2:end,:)];  
end
    
Y(:,4) = sum(Y(:,4:5),2);
ymod = Y(:,1:4);

end