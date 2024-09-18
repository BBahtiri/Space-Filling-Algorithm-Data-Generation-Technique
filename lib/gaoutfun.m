function [state,options,optchanged] = gaoutfun(options,state,flag)
optchanged = false;
format longEng
switch flag
    case 'init'
        disp(state.Best)
    case 'iter'
        sprintf('\nLast Improvement: \n')
        disp(state.LastImprovement)
        where_best = find(state.Score(:)==min(state.Best));
        if length(where_best) > 1
            where_best = where_best(end);
        end
        sprintf('\nBest: \n')
        disp(state.Best)
        sprintf('\nVariables for best solution: \n')
        disp(state.Population(where_best,:))
    case 'done'
        sprintf('\nLast Improvement: \n')
        disp(state.LastImprovement)
        where_best = find(state.Score(:)==min(state.Best));
        if length(where_best) > 1
            where_best = where_best(end);
        end
        sprintf('\nBest: \n')
        disp(state.Best)
        sprintf('\nVariables for best solution: \n')
        disp(state.Population(where_best,:))
end
end