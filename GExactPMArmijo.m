function [X,f,time,outiter,nfev,info] = GExactPMArmijo(n,X)

% Start timing

tic;

% Initial Parameters

tol        = 10^(-4);
ftol       = 10^(-4);

maxoutiter = 500;
maxinneriter = 1000;
stopit     = 2;
alphamin   = 10^(-10);
alphamax   = 10^(10);
stpmin     = 10^(-10);
sigma1     = 0.1;
sigma2     = 0.9;

% Counters

outiter = 0;
nfev    = 0;

% Compute the function value

[f] = evalf(X);
nfev = nfev + 1;

% Compute the gradient

[g] = evalg(X);

normg = norm(g,'fro');

normgX0 = normg;

% Define alpha

alpha = min(alphamax, max(alphamin, 1.0 / normg ) );

fprintf('------------------------------------------------\n')
fprintf('       GExactPM employing Armijo stepsize       \n')
fprintf('------------------------------------------------\n')
fprintf('Number of variables : %i \n',n)
fprintf('Tolerance for convergence: %.0e \n\n',tol)

while (1)
    
    % Print information
    
    if ( mod(outiter,10) == 0 )
        fprintf('\n')
        fprintf('%-5s   %-8s  %-8s\n','it','f','|x-xprev|/|xprev|')
    end
    if ( outiter == 0 )
        fprintf('%5d   %5.2e  %1s\n',outiter,f,'-')
    else
        fprintf('%5d   %5.2e  %8.2e\n',outiter,f,normXXprev)
    end
    
    % Project X - alpha * g onto the feasible set
    
    Xprev = X;
    
    Y = X - alpha * g;
    
    Y = ( Y + Y' ) / 2.0d0;
    
    % Compute the exact projection of Y
     
    [W] = spec_proj(Y);
    
    W = ( W + W' ) / 2.0;
    
    % Define the search direction
    
    d = W - X;
    
    % --------------------------------
    % Stopping criteria
    % --------------------------------
    
    % Test (unconstrained) optimality
    
    if ( normg / normgX0 <= 10^(-6) )
        info = 0;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('CPU time(s): %.1f \n',time)
        return
    end
    
    % Check the progress on the iterates
    
    if ( outiter > 0 && normXXprev <= tol  )
        stop = stop + 1;
    else
        stop = 0;
    end
    
   if ( stop == stopit )
        
        info = 1;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('CPU time(s): %.1f \n',time)
        return
    end
        
    % Test whether the number of iterations is exhausted
    
    if ( outiter == maxoutiter )
        info = 2;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('The number of maximum iterations was reached.\n')
        fprintf('CPU time(s): %.1f \n',time)
        
        return
    end
    
    % Test whehter the norm of d is too small
    
    if ( norm( d, 'fro' ) <= 10^(-8) )
        info = 4;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('CPU time(s): %.1f \n',time)
        return
    end
    
    % --------------------------------
    
    % Increment outiter
    
    outiter = outiter + 1;
    
    % Compute the step size
    
    gtd = g(:)' * d(:);
    
    stp = 1.0;
    
    while (1)
        
        Xtrial = X + stp * d;
        
        [ftrial] = evalf(Xtrial);
        nfev = nfev + 1;
        
        if ( ftrial <= f + ftol * stp * gtd )
            break
        end
        
        if ( stp <= stpmin ) 
            break
            disp('Warning: stp = stpmin in the backtracking procedure')
        end
			
        stpq = ( ( gtd / ( ( f - ftrial) / stp + gtd ) ) / 2.0 ) * stp;

        if ( stpq >= sigma1 * stp && stpq <= sigma2 * stp )
            stp = stpq;
        else
            stp = stp / 2.0;
        end
        
    end
    
    % Update X
    
    Xprev = X;
    
    X = Xtrial;
    
    % Compute norm(X,Xprev)/norm(Xprev)
        
    normXXprev = norm( X - Xprev,'fro' )/norm( Xprev,'fro' );
    
    % Compute the function value

    f = ftrial;

    % Compute the gradient
    
    gprev = g;

    [g] = evalg(X);

    normg = norm(g,'fro');

    % Define alpha

    s = X - Xprev;
    y = g - gprev;
    
    a = s(:)' * s(:);
    b = s(:)' * y(:);
    
    if ( b <= 10^(-12) )
        alpha = alphamax;
    else
        alpha = min(alphamax, max(alphamin, a / b ) );
    end
        
end