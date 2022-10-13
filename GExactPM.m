function [X,f,time,outiter,info] = GExactPM(n,X,L)

% Start timing

tic;

% Initial Parameters

tol        = 10^(-4);
maxoutiter = 1000;
stopit     = 2;

% Counters

outiter = 0;

fprintf('------------------------------------------------\n')
fprintf('      GExactPM employing constant stepsize      \n')
fprintf('------------------------------------------------\n')
fprintf('Number of variables : %i \n',n)
fprintf('Tolerance for convergence: %.0e \n\n',tol)

while (1)
    
    % Compute the function value
    
    [f] = evalf(X);
    
    % Compute the gradient
    
    [g] = evalg(X);
    
    normg = norm(g,'fro');
    
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
    
    % --------------------------------
    % Stopping criteria
    % --------------------------------
    
    % Test (unconstrained) optimality
    
    if ( normg <= 10^(-6) )
        info = 0;
        time = toc;
        
        % Print information
        
        [f] = evalf(X);
        
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
        
        [f] = evalf(X);
        
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
        
        [f] = evalf(X);
        
        fprintf('\n')
        fprintf('The number of maximum iterations was reached.\n')
        fprintf('CPU time(s): %.1f \n',time)
        
        return
    end
    
    % --------------------------------
    
    % Increment outiter
    
    outiter = outiter + 1;
    
    % Define alpha
    
    alpha = 0.9999 * ( 1.0  / L );
    
    % Compute the next iterate
    
    Xprev = X;
    
    Y = X - alpha * g;
    
    Y = ( Y + Y' ) / 2.0d0;
    
    % Compute the exact projection of Y
    
    [W] = spec_proj(Y);
    
    X = ( W + W' ) / 2.0;
    
    % Compute norm(X,Xprev)/norm(Xprev)
        
    normXXprev = norm( X - Xprev,'fro' )/norm( Xprev,'fro' );
end