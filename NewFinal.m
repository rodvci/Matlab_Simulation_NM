while true
    % Prompt the user to choose a root-finding method
    disp('Choose a root-finding method:');
    disp('1. GRAPHICAL METHOD');
    disp('2. INCREMENTAL METHOD');
    disp('3. BISECTION METHOD');
    disp('4. FALSE POSITION METHOD');
    disp('5. SIMPLE FIXED-POINT METHOD');
    disp('6. NEWTON-RAPHSON METHOD');
    disp('7. SECANT METHOD');
    method = input('Enter the number corresponding to your choice: ');

    switch method
        case 1
            % Graphing Method
            disp('GRAPHICAL METHOD');
            % Define the function
            disp('Enter your equation in terms of x (use MATLAB syntax ex. x^3+x^2-4*x+1):');
            eq_str = input('f(x) = ', 's');
            % Replace 'e' or 'Euler' with the numerical value
            eq_str = strrep(eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
            eq_str = strrep(eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value
            f = str2func(['@(x)' eq_str]);

            % Input initial value for x and step size
            valX = input('Enter the initial x value: ');
            step = 0.2; % Step size for graphing method

            % Initialize variables
            currX = valX + step;
            currY = f(currX);
            prevY = f(valX);
            signChange = false;

            % Store iteration information
            iteration_info = [];
            iterations = 0;

            while ~signChange
                iterations = iterations + 1;
                % Store current X and Y values
                iteration_info(end+1, :) = [currX, currY];

                if prevY * currY < 0
                    % Root finding using false position method
                    tolerance = 0.001;
                    lowerBound = currX - step;
                    upperBound = currX;
                    root = 0;

                    while abs(upperBound - lowerBound) > tolerance
                        root = (lowerBound * f(upperBound) - upperBound * f(lowerBound)) / (f(upperBound) - f(lowerBound));
                        fA = f(lowerBound);
                        fC = f(root);

                        if fC == 0.0
                            break;
                        elseif fA * fC < 0
                            upperBound = root;
                        else
                            lowerBound = root;
                        end
                    end

                    rootVal = root;
                    signChange = true; % Exit the loop as root is found
                else
                    prevY = currY;
                    currX = currX + step;
                    currY = f(currX);
                end
            end

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the iteration information in a table
            fprintf('Iteration information:\n');
            iteration_table = array2table(iteration_info, 'VariableNames', {'X', 'Y'});
            disp(iteration_table);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the root
            fprintf('Root found at:\n');
            disp(rootVal);

            % Plot the function and the root
            clf;
            fplot(f, [valX-1, currX+1]);
            hold on;
            plot(rootVal, f(rootVal), 'ro', 'MarkerSize', 8); % Plot the root
            title('Graphing Method Visual');
            xlabel('x');
            ylabel('f(x)');
            grid on;
            legend('Function', 'Root');

        case 2
            % Incremental Method
            disp('INCREMENTAL METHOD');
            % Define the function
            disp('Enter your equation in terms of x (use MATLAB syntax ex. x^3+x^2-4*x+1):');
            eq_str = input('f(x) = ', 's');
            % Replace 'e' or 'Euler' with the numerical value
            eq_str = strrep(eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
            eq_str = strrep(eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value
            f = str2func(['@(x)' eq_str]);

            % Input range start [a], deltaX and tolerance
            a = input('Enter the xL: ');
            deltaX = input('Enter the deltaX: ');
            tol = 0.001; % Set tolerance to 0.001

            % Initialize array to store roots
            roots = zeros(1, 1000); % Preallocate for a maximum of 1000 roots

            % Counter for the number of roots found
            root_count = 0;

            % Store iteration information
            iteration_info = [];

            % Initialize variables for the incremental search method
            xL = a;
            xU = xL + deltaX;
            iterations = 0;
            prevXU = 0;
            Error = abs(xU - prevXU) / abs(xU) * 100;

            while Error > tol
                iterations = iterations + 1;
                fxL = f(xL);
                fxU = f(xU);
                positive = fxL * fxU;

                % Calculate error
                if iterations == 1
                    currentError = NaN;
                else
                    currentError = abs(xU - prevXU) / abs(xU) * 100;
                end

                % Store iteration information
                iteration_info(end+1, :) = [iterations, xL, deltaX, xU, fxL, fxU, currentError];

                if Error < tol && fxL * fxU == 0
                    root_count = root_count + 1;
                    roots(root_count) = xU;
                    break;
                end

                if fxL * fxU < 0
                    % xL = xL;
                    deltaX = deltaX / 10.0;
                else
                    xL = xU;
                end

                prevXU = xU;
                xU = xL + deltaX;
                Error = abs(xU - prevXU) / abs(xU) * 100;
            end

            % Final root
            root_count = root_count + 1;
            roots(root_count) = xU;

            % Trim the excess preallocated space
            roots = roots(1:root_count);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the iteration information in a table
            fprintf('Iteration information:\n');
            iteration_table = array2table(iteration_info, 'VariableNames', {'Iteration', 'xl', 'deltaX', 'xu', 'f(xl)', 'f(xu)', 'Error'});
            disp(iteration_table);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the roots
            fprintf('Roots:\n');
            disp(roots);

            % Plot the roots
            clf;
            fplot(f, [a-1, xU+1]);
            hold on;
            plot(roots, arrayfun(f, roots), 'ro', 'MarkerSize', 8); % Using arrayfun to apply function to each root
            title('Incremental Method Graph Visual');
            xlabel('x');
            ylabel('f(x)');
            grid on;
            legend('Function', 'Roots');

        case 3
            % Bisection Method
            disp('BISECTION METHOD');
            % Define the function
            disp('Enter your equation in terms of x (use MATLAB syntax ex. x^3+x^2-4*x+1):');
            eq_str = input('f(x) = ', 's');
            % Replace 'e' or 'Euler' with the numerical value
            eq_str = strrep(eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
            eq_str = strrep(eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value
            f = str2func(['@(x)' eq_str]);

            % Input range [a, b] and tolerance
            a = input('Enter the xL: ');
            b = input('Enter the xU: ');
            tol = 0.001; % Set tolerance to 0.001

            % Initialize array to store roots
            roots = zeros(1, 1000); % Preallocate for a maximum of 1000 roots

            % Counter for the number of roots found
            root_count = 0;

            % Store iteration information
            iteration_info = [];

            % Initialize variables for the bisection method
            xL = a;
            xU = b;
            iteration = 0;
            prevXR = 0;

            while (xU - xL) > tol
                iteration = iteration + 1;
                xR = (xL + xU) / 2.0;
                fxL = f(xL);
                fxR = f(xR);
                fxU = f(xU);
                positive = fxL * fxU;

                % Format numbers for display (not strictly necessary for calculations)
                decimalnumA = xL;
                decimalnumB = xR;
                decimalnumC = xU;
                decimalnumFA = fxL;
                decimalnumFB = fxR;
                decimalnumFC = fxU;

                % Calculate error
                if iteration == 1
                    Error = NaN;
                else
                    Error = abs((xR - prevXR) / xR) * 100;
                end

                % Store iteration information
                iteration_info(end+1, :) = [iteration, decimalnumA, decimalnumB, decimalnumC, decimalnumFA, decimalnumFB, decimalnumFC, Error];

                if Error < tol && fxL * fxU == 0
                    root_count = root_count + 1;
                    roots(root_count) = xR;
                    break;
                end

                if fxL * fxR < 0
                    xU = xR;
                else
                    xL = xR;
                end

                prevXR = xR;
            end

            % Final root
            root_count = root_count + 1;
            roots(root_count) = xR;

            % Trim the excess preallocated space
            roots = roots(1:root_count);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the iteration information in a table
            fprintf('Iteration information:\n');
            iteration_table = array2table(iteration_info, 'VariableNames', {'Iteration', 'xL', 'xR', 'xU', 'f(xL)', 'f(xR)', 'f(xU)', 'Error'});
            disp(iteration_table);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the roots
            fprintf('Roots:\n');
            disp(roots);

            % Plot the roots
            clf;
            fplot(f, [a-1, b+1]);
            hold on;
            plot(roots, arrayfun(f, roots), 'ro', 'MarkerSize', 8); % Using arrayfun to apply function to each root
            title('Bisection Method Graph Visual');
            xlabel('x');
            ylabel('f(x)');
            grid on;
            legend('Function', 'Roots');

        case 4
            % Regula Falsi Method
            disp('FALSE POSITION METHOD');
            % Define the function
            disp('Enter your equation in terms of x (use MATLAB syntax ex. x^3+x^2-4*x+1):');
            eq_str = input('f(x) = ', 's');
            % Replace 'e' or 'Euler' with the numerical value
            eq_str = strrep(eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
            eq_str = strrep(eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value
            f = str2func(['@(x)' eq_str]);

            % Input range [a, b] and tolerance
            a = input('Enter the xL: ');
            b = input('Enter the xU: ');
            tol = 0.001; % Set tolerance to 0.001

            % Initialize array to store roots
            roots = zeros(1, 1000); % Preallocate for a maximum of 1000 roots

            % Counter for the number of roots found
            root_count = 0;

            % Store iteration information
            iteration_info = [];

            % Initialize variables for the Regula Falsi method
            xL = a;
            xU = b;
            iteration = 0;
            prevXU = 0;

            while (xU - xL) > tol
                iteration = iteration + 1;
                fxL = f(xL);
                fxU = f(xU);

                % Calculate xR using Regula Falsi formula
                xR = (xL * fxU - xU * fxL) / (fxU - fxL);
                fxR = f(xR);

                % Calculate error
                if iteration == 1
                    Error = NaN;
                else
                    Error = abs((xR - prevXU) / xR) * 100;
                end

                % Store iteration information
                iteration_info(end+1, :) = [iteration, xL, xU, xR, fxL, fxU, fxR, Error];

                if Error < tol && fxL * fxU == 0
                    root_count = root_count + 1;
                    roots(root_count) = xR;
                    break;
                end

                if fxL * fxR < 0
                    xU = xR;
                else
                    xL = xR;
                end

                prevXU = xR;
            end

            % Final root
            root_count = root_count + 1;
            roots(root_count) = xR;

            % Trim the excess preallocated space
            roots = roots(1:root_count);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the iteration information in a table
            fprintf('Iteration information:\n');
            iteration_table = array2table(iteration_info, 'VariableNames', {'Iteration', 'xL', 'xU', 'xR', 'f(xL)', 'f(xU)', 'f(xR)', 'Error'});
            disp(iteration_table);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the roots
            fprintf('Roots:\n');
            disp(roots);

            % Plot the roots
            clf;
            fplot(f, [a-1, b+1]);
            hold on;
            plot(roots, arrayfun(f, roots), 'ro', 'MarkerSize', 8); % Using arrayfun to apply function to each root
            title('False Position Method Graph Visual');
            xlabel('x');
            ylabel('f(x)');
            grid on;
            legend('Function', 'Roots');

        case 5
            %FIXED POINT ITERATION
            disp('SIMPLE FIXED POINT-METHOD');
        % Define the function
        disp('Enter your equation in terms of x (use MATLAB syntax ex. x^3+x^2-4*x+1):')
        eq_str = input('f(x) = ', 's');

        % Replace 'e' or 'Euler' with the numerical value
        eq_str = strrep(eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
        eq_str = strrep(eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value

        f = str2func(['@(x)' eq_str]);

        % Define the function to differentiate
        differentiateFunction = @(xValue) (f(xValue + 0.000001) - f(xValue)) / 0.000001;

        % Define the fixed-point iteration function g(x)
        g = @(x) x - f(x) / differentiateFunction(x);

        % Input initial guess and tolerance
        x0 = input('Enter the initial guess x0: ');
        tol = 0.001; % Set tolerance to 0.001

        % Initialize array to store iteration data
        iteration_info = [];

        % Initialize variables for the fixed-point iteration method
        iterations = 0;
        error = Inf;

        while error > tol
            iterations = iterations + 1;
            x1 = g(x0); % Fixed-point iteration: x1 = g(x0)

            % Calculate error
            if iterations == 1
                prevX0 = NaN;
            else
                prevX0 = iteration_info(end, 2);
            end
            error = abs(x1 - prevX0) / abs(x1) * 100;

            % Format numbers for display (not strictly necessary for calculations)
            decimalnumX0 = x0;
            decimalnumX1 = x1;
            decimalnumError = error;

            % Store iteration information
            iteration_info(end+1, : ) = [iterations, decimalnumX0, decimalnumX1, decimalnumError];

            x0 = x1;
        end

        % Print a separator line
        fprintf('----------------------------------------\n');

        % Display the iteration information in a table
        fprintf('Iteration information:\n');
        iteration_table = array2table(iteration_info, 'VariableNames', {'Iteration', 'x0', 'x1', 'Error'});
        disp(iteration_table);

        % Print a separator line
        fprintf('----------------------------------------\n');

        % Display the root
        fprintf('Root:\n');
        fprintf('x = %f\n', x1);

        % Plot the function and the root
        clf;
        xL = min(iteration_info(:, 2)); % Minimum x0 value from iterations
        xU = max(iteration_info(:, 2)); % Maximum x0 value from iterations
        fplot(f, [xL-1, xU+1]);
        hold on;
        plot(x1, f(x1), 'ro', 'MarkerSize', 8); % Plot the root
        title('Fixed Point Iteration Method');
        xlabel('x');
        ylabel('f(x)');
        grid on;
        legend('Function', 'Root');


                case 6
                    % Newton-Raphson Method
                    disp('NEWTON-RAPHSON METHOD');
                    % Define the function and its derivative
                    disp('Enter your equation in terms of x (use MATLAB syntax ex. x^3+x^2-4*x+1):');
                    eq_str = input('f(x) = ', 's');
                    disp('Enter the derivative of your equation (use MATLAB syntax ex. 3*x^2+2*x-4):');
                    d_eq_str = input('f''(x) = ', 's');
                    % Replace 'e' or 'Euler' with the numerical value
                    eq_str = strrep(eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
                    eq_str = strrep(eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value
                    d_eq_str = strrep(d_eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
                    d_eq_str = strrep(d_eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value
                    f = str2func(['@(x)' eq_str]);
                    df = str2func(['@(x)' d_eq_str]);

                    % Input initial guess and tolerance
                    x0 = input('Enter the initial guess: ');
                    tol = 0.001; % Set tolerance to 0.001

                    % Initialize variables for the Newton-Raphson method
                    iterations = 0;
                    x_old = x0;
                    error = tol + 1;

                    % Store iteration information
                    iteration_info = [];

                    while error > tol
                        iterations = iterations + 1;
                        x_new = x_old - f(x_old) / df(x_old);
                        error = abs((x_new - x_old) / x_new) * 100;

                        % Store iteration information
                        iteration_info(end+1, :) = [iterations, x_old, x_new, error];

                        x_old = x_new;
                    end

                    % Print a separator line
                    fprintf('----------------------------------------\n');
                    % Display the iteration information in a table
                    fprintf('Iteration information:\n');
                    iteration_table = array2table(iteration_info, 'VariableNames', {'Iteration', 'x_old', 'x_new', 'Error'});
                    disp(iteration_table);

                    % Print a separator line
                    fprintf('----------------------------------------\n');
                    % Display the root
                    fprintf('Root found at:\n');
                    disp(x_new);

                    % Plot the function and the root
                    clf;
                    fplot(f, [x0-1, x0+1]);
                    hold on;
                    plot(x_new, f(x_new), 'ro', 'MarkerSize', 8); % Plot the root
                    title('Newton-Raphson Method Graph Visual');
                    xlabel('x');
                    ylabel('f(x)');
                    grid on;
                    legend('Function', 'Root');

        case 7
            % Secant Method
            disp('SECANT METHOD');
            % Define the function
            disp('Enter your equation in terms of x (use MATLAB syntax ex. x^3+x^2-4*x+1):');
            eq_str = input('f(x) = ', 's');
            % Replace 'e' or 'Euler' with the numerical value
            eq_str = strrep(eq_str, 'e', num2str(exp(1))); % 'e' to numerical value
            eq_str = strrep(eq_str, 'Euler', num2str(exp(1))); % 'Euler' to numerical value
            f = str2func(['@(x)' eq_str]);

            % Input initial guesses and tolerance
            x0 = input('Enter the first initial guess (x0): ');
            x1 = input('Enter the second initial guess (x1): ');
            tol = 0.001; % Set tolerance to 0.001

            % Initialize variables for the Secant method
            iterations = 0;
            error = tol + 1;

            % Store iteration information
            iteration_info = [];

            while error > tol
                iterations = iterations + 1;
                fx0 = f(x0);
                fx1 = f(x1);

                x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
                error = abs((x_new - x1) / x_new) * 100;

                % Store iteration information
                iteration_info(end+1, :) = [iterations, x0, x1, x_new, fx0, fx1, error];

                x0 = x1;
                x1 = x_new;
            end

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the iteration information in a table
            fprintf('Iteration information:\n');
            iteration_table = array2table(iteration_info, 'VariableNames', {'Iteration', 'x0', 'x1', 'x_new', 'f(x0)', 'f(x1)', 'Error'});
            disp(iteration_table);

            % Print a separator line
            fprintf('----------------------------------------\n');
            % Display the root
            fprintf('Root found at:\n');
            disp(x_new);

            % Plot the function and the root
            clf;
            fplot(f, [x0-1, x1+1]);
            hold on;
            plot(x_new, f(x_new), 'ro', 'MarkerSize', 8); % Plot the root
            title('Secant Method Graph Visual');
            xlabel('x');
            ylabel('f(x)');
            grid on;
            legend('Function', 'Root');
            
        otherwise
            disp('Invalid choice. Please enter a number between 1 and 7.');
    end

    % Ask if the user wants to calculate again
    choice = input('Do you want to calculate again? [y]es/[n]o: ', 's');
    if ~strcmpi(choice, 'y') && ~strcmpi(choice, 'yes')
        disp('EXIT.....');
        break; % Exit the loop if the user does not want to calculate again
    end
end