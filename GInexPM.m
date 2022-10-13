function [X,f,time,outiter,info] = GInexPM(n,X,L,erropt)

% Start timing

tic;

% Initial Parameters

tol          = 10^(-4);
a            = 1;
k            = 1;
maxoutiter   = 1000;
maxinneriter = 1000;
stopit       = 2;

% Counters

outiter  = 0;
easyProj = 0;

fprintf('------------------------------------------------\n')
fprintf('       GInexPM employing constant stepsize      \n')
fprintf('------------------------------------------------\n')
fprintf('Number of variables : %i \n',n)
fprintf('Tolerance for convergence: %.0e \n\n',tol)

while (1)
    
    % Compute the function value
    
    [f] = evalf(X);
    
    % Compute the gradient
    
    [g] = evalg(X);
    
    normg = norm(g,'fro');
    
    if ( outiter == 0 )
        normgX0 = normg;
    end
    
    % Print information
    
    if ( mod(outiter,10) == 0 )
        fprintf('\n')
        fprintf('%-5s   %-8s  %-5s %2s  %-8s\n','it','f','k','IS','|x-xprev|/|xprev|')
    end
    if ( outiter == 0 )
        fprintf('%5d   %5.2e  %1s      %1s       %1s\n',outiter,f,'-','-','-')
    else
        fprintf('%5d   %5.2e  %-5d %2d    %8.2e\n',outiter,f,k,infoProj,normXXprev)
    end
    
    % --------------------------------
    % Stopping criteria
    % --------------------------------
    
    % Test (unconstrained) optimality
    
    if ( normg / normgX0 <= 10^(-6) )
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
        
        [f] = evalf(X);
        
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
        
        [f] = evalf(X);
        
        fprintf('\n')
        fprintf('The number of maximum iterations was reached.\n')
        fprintf('CPU time(s): %.1f \n',time)
        
        return
    end
    
    % --------------------------------
    
    % Increment outiter
    
    outiter = outiter + 1;
    
    % Define gamma
        
    anablaf2 = a / normg^2;
    
    gamma(2) = min( anablaf2 / 2 , 0.49995 );
    gamma(1) = anablaf2 - gamma(2);
    
    if ( erropt == 1 ||  erropt == 2 || erropt == 3 || erropt == 4 )
        gamma(3) = 0.1;
    end
    
    if ( erropt == 5 )
        gamma(3) = 0.0;
    end
    
    if ( erropt == 6 )
        gamma(1) = anablaf2;
        gamma(2) = 0.0;
        gamma(3) = 0.0;
    end
    
    if ( erropt == 7 )
        gamma(1) = 0.0;
        gamma(2) = min( anablaf2 , 0.49995 );
        gamma(3) = 0.0;
    end
    
    if ( erropt == 8 )
      gamma(1) = 0.0;
      gamma(2) = 0.0;
      gamma(3) = 0.0;
    end
    
    % Define alpha
    
    alpha = 0.9999 * ( 1.0 - 2.0 * gamma(3) ) / L;
    
    % Compute the next iterate
    
    Xprev = X;
    
    Y = X - alpha * g;
    
    Y = sparse( ( Y + Y' ) / 2.0d0 );
    
    % Set the rank-k
    
    if ( outiter > 1 )
        if ( ( infoProj == 0 || infoProj == 1 ) && outiterProj == 1 )
            easyProj = easyProj + 1;
            if ( easyProj == 5 )
                k = max( k - 1, 1 );
                easyProj = 0;
            end
        else
            easyProj = 0;
        end

    end
    
    % Compute the inexact projection of Y relative to X

    [W,k,outiterProj,infoProj] = SpecDp(n,k,X,Y,gamma,erropt);

    
    X = ( W + W' ) / 2.0;
    
    % Compute norm(X,Xprev)/norm(Xprev)
        
    normXXprev = norm( X - Xprev,'fro' )/norm( Xprev,'fro' );
    
    % Update a
    
    if ( outiter == 1 )
        a = 2.0 - 1.0 / ( log(2) );
    else
        a = 1.0 / ( log( outiter ) ) - 1.0 / ( log( outiter + 1 ) );
    end
    
    a = 10^2 * a;
end