function obj = objective_function(Para_set, data)
% sum-of-square error
obj = 0;

for i  = 1:length(data)
    ymodel = Dolmodel(data{i},Para_set);
    ydata = data{i}.ydata(:,[3,2,4,5]);  % adjust the columns of ydata to fit the order in ymodel
    obj = sum([obj, sum((ymodel-ydata).^2)]);
end

end
