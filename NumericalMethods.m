% Define the function
method = input('Enter the method to use (1 for Graphical, 2 for Incremental, 3 for Bisection, 4 for Regula-Falsi, 5 for Fixed-Point, 6 for Newton Rhapson, 7 for Secant Method): ');

root = NaN; % Initialize root
function x1 = g(x0, equation)
    % Example: Let's choose g(x) = x - f(x) / f'(x), similar to the Newton-Raphson method
    fx = evaluateFunction(equation, x0);
    fDashX = numerical_derivative(equation, x0);

    % Avoid division by zero
    if fDashX == 0
        error('Derivative is zero. Division by zero error.');
    end

    x1 = x0 - (fx / fDashX);
end

function df = numerical_derivative(func, x)
    % Numerical differentiation using central difference
    h = 1e-5; % Small step size
    df = (func(x + h) - func(x - h)) / (2 * h);
end

function displayResult(data, root)
    % Remove empty rows from data
    data = data(~cellfun('isempty', data(:, 1)), :);
    T = cell2table(data, 'VariableNames', {'Iteration', 'Xi', 'Xi+1', 'et%'});
    disp(T);
    disp(['Root found with Fixed-Point Iteration Method: ', num2str(root)]);
end

function f = evaluateFunction(equation, x)
    try
        f = equation(x);
    catch
        error('Error evaluating the function.');
    end
end
if method == 1
    string_equation = input('Enter the function (e.g., "exp(-x) - x"): ', 's');
    equation = str2func(['@(x)', string_equation]);
    startPoint = input('Enter the starting point (e.g., 0.01): ');
    endPoint = input('Enter the ending point (e.g., 0.01): ');
    increment = input('Enter the increment between numbers (e.g., 0.01): ');
    X=[];
    pointCounter=startPoint;
    Y=[];
    iteration = [];
    iter=0;
    while pointCounter<endPoint
        z = equation(pointCounter);
        iter=iter+1;
        iteration = [iteration;iter];
        X = [X;z];
        Y = [Y;pointCounter];
        pointCounter=pointCounter+increment;
    end
    root = fzero(equation, [startPoint, endPoint]);
    T = table(iteration,Y,X);
    disp(T);
    plot(Y,X);
    hold on;
    if ~isnan(root)
        plot(root,0,'r+', 'MarkerSize', 10, 'LineWidth', 1)
        title('Graphical Method');
        xlabel('X');
        ylabel('Y');
        grid on;
    end
    disp(['Root found with Graphical Method: ', num2str(root)]);
elseif method == 2
    string_equation = input('Enter the function (e.g., "exp(-x) - x"): ', 's');
    equation = str2func(['@(x)', string_equation]);
    startPoint = input('Enter the starting point (e.g., 0.01): ');
    endPoint = input('Enter the ending point (e.g., 0.01): ');
    increment = input('Enter the ΔX (e.g., 0.01): ');
    criterionStop = input('Enter stopping criterion (e.g., 0.0001): ');
    xl = startPoint
    xu = startPoint+increment;
    Xl=[];
    deltax = increment;
    DeltaX=[];
    Xu=[];
    iteration = [];
    iter=0;
    FXL = [];
    FXU = [];
    eA = [];
    fxl = 0;
    fxu = 0;
    ea = 1;
    fXufXl = [];
    xu1=0;
    X=[];
    pointCounter=startPoint;
    Y=[];
    iteration1 = [];
    iter1=0;
    while pointCounter<endPoint
        z = equation(pointCounter);
        iter=iter+1;
        iteration1 = [iteration1;iter1];
        X = [X;z];
        Y = [Y;pointCounter];
        pointCounter=pointCounter+0.01;
    end
    plot(Y,X);
    hold on;
    while criterionStop<ea
        iter=iter+1;
        iteration = [iteration;iter];
        fxl = equation(xl);
        fxu = equation(xu);
        fxufxl = fxl*fxu;
        ea = (abs(xu-xu1)/xu)*100;
        if ea<criterionStop
            root = xu;
        end
        DeltaX = [DeltaX;deltax];
        Xl = [Xl;xl];
        Xu = [Xu;xu];
        FXL = [FXL;fxl];
        FXU = [FXU;fxu];
        eA = [eA;ea];
        fXufXl = [fXufXl;fxufxl];
        xu1=xu;
        if fxufxl>0
            xl=xu;
            xu=xl+deltax;
        else
            deltax = deltax/10;
            xu = xl+deltax;
        end
    end
    T = table(iteration,Xl,DeltaX,Xu,FXL,FXU,eA,fXufXl);
    disp(T);
    hold on;
    if ~isnan(root)
        plot(root,0,'r+', 'MarkerSize', 10, 'LineWidth', 1)
        title('Incremental Method');
        xlabel('X');
        ylabel('Y');
        grid on;
    end
    disp(['Root found with Incremental Method: ', num2str(root)]);
elseif method == 3
    string_equation = input('Enter the function (e.g., "exp(-x) - x"): ', 's');
    equation = str2func(['@(x)', string_equation]);
    startPoint = input('Enter the starting point (e.g., 0.01): ');
    endPoint = input('Enter the ending point (e.g., 0.01): ');
    criterionStop = input('Enter stopping criterion (e.g., 0.0001): ');
    xl = startPoint;
    xu = endPoint;
    xr = (startPoint+endPoint)/2;
    Xl=[];
    Xu=[];
    Xr=[];
    iteration = [];
    iter=0;
    FXL = [];
    FXR = [];
    eA = [];
    fxl = 0;
    fxr = 0;
    ea = 1;
    fXrfXl = [];
    xr1=0;
    X=[];
    pointCounter=startPoint;
    Y=[];
    iteration1 = [];
    iter1=0;
    while pointCounter<endPoint
        z = equation(pointCounter);
        iter=iter+1;
        iteration1 = [iteration1;iter1];
        X = [X;z];
        Y = [Y;pointCounter];
        pointCounter=pointCounter+0.01;
    end
    plot(Y,X);
    hold on;
    while criterionStop<ea
        iter=iter+1;
        iteration = [iteration;iter];
        fxl = equation(xl);
        fxr = equation(xr);
        fxrfxl = fxl*fxr;
        ea = (abs(xr-xr1)/xr)*100;
        if ea<criterionStop
            root = xr;
        end
        Xr = [Xr;xr];
        Xl = [Xl;xl];
        Xu = [Xu;xu];
        FXL = [FXL;fxl];
        FXR = [FXR;fxr];
        eA = [eA;ea];
        fXrfXl = [fXrfXl;fxrfxl];
        xr1=xr;
        if fxrfxl>0
            xl=xr;
            xr = (xl+xu)/2;
        else
            xu = xr;
            xr = (xl+xu)/2;
        end
    end
    T = table(iteration,Xl,Xr,Xu,FXL,FXR,eA,fXrfXl);
    disp(T);
    hold on;
    if ~isnan(root)
        plot(root,0,'r+', 'MarkerSize', 10, 'LineWidth', 1)
        title('Bisection Method');
        xlabel('X');
        ylabel('Y');
        grid on;
    end
    disp(['Root found with Bisection Method: ', num2str(root)]);
elseif method == 4
    string_equation = input('Enter the function (e.g., "exp(-x) - x"): ', 's');
    equation = str2func(['@(x)', string_equation]);
    startPoint = input('Enter the Xl (e.g., 0.01): ');
    endPoint = input('Enter the Xu (e.g., 0.01): ');
    criterionStop = input('Enter stopping criterion (e.g., 0.0001): ');
    xl = startPoint;
    xu = endPoint;
    xr = 0;
    Xl=[];
    Xu=[];
    Xr=[];
    iteration = [];
    iter=0;
    FXL = [];
    FXR = [];
    FXU = [];
    eA = [];
    fxl = 0;
    fxr = 0;
    fxu = 0;
    ea = 1;
    fXrfXl = [];
    xr1=0;
    X=[];
    pointCounter=startPoint;
    Y=[];
    iteration1 = [];
    iter1=0;
    while pointCounter<endPoint
        z = equation(pointCounter);
        iter=iter+1;
        iteration1 = [iteration1;iter1];
        X = [X;z];
        Y = [Y;pointCounter];
        pointCounter=pointCounter+0.01;
    end
    plot(Y,X);
    hold on;
    while criterionStop<ea
        iter=iter+1;
        iteration = [iteration;iter];
        fxl = equation(xl);
        fxu = equation(xu);
        xr = ((xu*fxl)-(xl*fxu))/(fxl-fxu);
        fxr = equation(xr);
        fxrfxl = fxl*fxr;
        ea = (abs(xr-xr1)/xr)*100;
        if ea<criterionStop
            root = xr;
        end
        Xr = [Xr;xr];
        Xl = [Xl;xl];
        Xu = [Xu;xu];
        FXL = [FXL;fxl];
        FXU = [FXU;fxu];
        FXR = [FXR;fxr];
        eA = [eA;ea];
        fXrfXl = [fXrfXl;fxrfxl];
        xr1=xr;
        if fxrfxl>0
            xu = xr;
        else
            xl = xr;
        end
    end
    T = table(iteration,Xl,Xu,Xr,FXL,FXU,FXR,eA,fXrfXl);
    disp(T);
    hold on;
    if ~isnan(root)
        plot(root,0,'r+', 'MarkerSize', 10, 'LineWidth', 1)
        title('Regula-Falsi Method');
        xlabel('X');
        ylabel('Y');
        grid on;
    end
    disp(['Root found with Regula-Falsi Method: ', num2str(root)]);
elseif method == 5
    string_equation = input('Enter the function (e.g., "exp(-x)-x"): ', 's');
     try
        equation = str2func(['@(x)' string_equation]);
        x0 = input('Enter the initial guess: ');
        tol = input('Enter the allowed error: ');

        max_iterations = 1000;
        iterations = 0;
        error = inf;
        data = cell(max_iterations, 4);

        while error > tol && iterations < max_iterations
            iterations = iterations + 1;
            x1 = g(x0, equation);  % Fixed-point iteration: x1 = g(x0)

            if iterations == 1
                error = inf;
            else
                error = abs((x1 - str2double(data{iterations-1, 3})) / x1) * 100;
            end

            data{iterations, 1} = iterations;
            data{iterations, 2} = num2str(x0, '%.9f');
            data{iterations, 3} = num2str(x1, '%.9f');
            data{iterations, 4} = num2str(error, '%.9f');

            x0 = x1;
        end

        root = x0;
        displayResult(data, root);

        % Plotting the function with the root marked
        x_vals_fixedpoint = linspace(root - 1, root + 1, 1000);
        y_vals_fixedpoint = arrayfun(@(x) equation(x), x_vals_fixedpoint);

        figure;
        plot(x_vals_fixedpoint, y_vals_fixedpoint, 'LineWidth', 2);
        hold on;
        plot(root, equation(root), 'ro', 'MarkerSize', 10);
        title('Simple Fixed-Point Method');
        xlabel('X');
        ylabel('Y');
        grid on;
        hold off;
    catch ME
        disp(ME.message);
   
    end

elseif method == 6
    string_equation = input('Enter the function (e.g., "exp(-x) - x"): ', 's');
    equation = str2func(['@(x)', string_equation]);
    
    % Define symbolic variable x
    syms x;
    % Convert the string to a symbolic function
    equation_sym = str2sym(string_equation);
    % Compute the symbolic derivative
    equationNew_sym = diff(equation_sym, x);
    % Convert the symbolic derivative to a function handle
    equationNew = matlabFunction(equationNew_sym);

    startPoint = input('Enter the Xi (e.g., 0.01): ');
    endPoint = input('Enter the limit of Xi(e.g., 0.01):');
    criterionStop = input('Enter stopping criterion (e.g., 0.0001): ');
    xi = startPoint;
    Xi=[];
    iteration = [];
    iter=0;
    FXI = [];
    FFXI = [];
    FXU = [];
    eT=[];
    et=1;
    fxi = 0;
    ffxi = 0;
    xi1=0;
    X=[];
    pointCounter=startPoint;
    Y=[];
    iteration1 = [];
    iter1=0;
    root = fzero(equation, [startPoint, endPoint]);
    while pointCounter<endPoint
        z = equation(pointCounter);
        iter1=iter1+1;
        iteration1 = [iteration1;iter1];
        X = [X;z];
        Y = [Y;pointCounter];
        pointCounter=pointCounter+0.01;
    end
    plot(Y,X);
    hold on;
    while criterionStop<et
        iter=iter+1;
        iteration = [iteration;iter];
        fxi = equation(xi);
        ffxi = equationNew(xi);
        et = (abs(xi-root)/root)*100;
        if et<criterionStop
            root = xi;
        end
        Xi = [Xi;xi];
        FXI = [FXI;fxi];
        FFXI = [FFXI;ffxi];
        eT = [eT;et];
        xi1=xi;
        xi = xi-(fxi/ffxi);
    end
    T = table(iteration,Xi,FXI,FFXI,eT);
    disp(T);
    hold on;
    if ~isnan(root)
        plot(root,0,'r+', 'MarkerSize', 10, 'LineWidth', 1)
        title('Newton Raphson Method');
        xlabel('X');
        ylabel('Y');
        grid on;
    end
    disp(['Root found with Newton Raphson Method: ', num2str(root)]);
elseif method == 7
    string_equation = input('Enter the function (e.g., "exp(-x) - x"): ', 's');
    equation = str2func(['@(x)', string_equation]);
    startPoint = input('Enter the X0 (e.g., 0.01): ');
    endPoint = input('Enter the X1 (e.g., 0.01): ');
    criterionStop = input('Enter stopping criterion (e.g., 0.0001): ');
    x0 = startPoint;
    x1 = endPoint;
    xi = 0;
    et = 1;
    eT = [];
    X0=[];
    X1=[];
    Xi = [];
    iteration = [];
    iter=0;
    FX0 = [];
    FX1 = [];
    %FXI = [];
    fx0 = 0;
    fx1 = 0;
    fxi = 0;
    xr1=0;
    pointCounter=startPoint;
    Y=[];
    X=[];
    iteration1 = [];
    iter1=0;
    root = fzero(equation, [startPoint, endPoint]);
    while pointCounter<endPoint
        z = equation(pointCounter);
        iter=iter+1;
        iteration1 = [iteration1;iter1];
        X = [X;z];
        Y = [Y;pointCounter];
        pointCounter=pointCounter+0.01;
    end
    plot(Y,X);
    hold on;
    while criterionStop<et
        iter=iter+1;
        iteration = [iteration;iter];
        fx0 = equation(x0);
        fx1 = equation(x1);
        %fxi = equation(xi);
        xi = x1-((fx1*(x0-x1))/(fx0-fx1));
        
        et = (abs(xi-root)/root)*100;
        if et<criterionStop
            root = xi;
        end
        Xi = [Xi;xi];
        X1 = [X1;x1];
        X0 = [X0;x0];
        FX0 = [FX0;fx0];
        FX1 = [FX1;fx1];
        %FXI = [FXI;fxi];
        eT = [eT;et];
        x0=x1;
        x1=xi;
    end
    T = table(iteration,X0,X1,Xi,FX0,FX1,eT);
    disp(T);
    hold on;
    if ~isnan(root)
        plot(root,0,'r+', 'MarkerSize', 10, 'LineWidth', 1)
        title('Secant Method');
        xlabel('X');
        ylabel('Y');
        grid on;
    end
    disp(['Root found with Secant Method: ', num2str(root)]);
else

end