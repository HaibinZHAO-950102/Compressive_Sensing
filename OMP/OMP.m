function [x, Out] = OMP(A,y,varargin)

    % Test for number of required parametres
    if (nargin-length(varargin)) ~= 2
        error('Wrong number of required parameters');
    end
    %--------------------------------------------------------------
    % Set parameters to their defaults
    %--------------------------------------------------------------
    opts.tol = 1e-3;
    opts.maxIter = 1e3;
    %--------------------------------------------------------------
    % Set parameters to user specified values
    %--------------------------------------------------------------
    if (rem(length(varargin),2)==1)
        error('Options should be given in pairs');
    else
        for i=1:2:(length(varargin)-1)
            switch upper(varargin{i})
                case 'STOPTOLERANCE'
                    opts.tol = varargin{i+1};
                case 'MAXITER'
                    opts.maxit = varargin{i+1};
                otherwise
                    error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
    [m_A, n] = size(A);
    [m_y, n_y] = size(y);
    if n_y ~= 1
        error(['y must be a column vector']);
    end
    if m_A ~= m_y
        error(['A and y must have same # of row'])
    end
    % --------------------------------------------------------------
    x = zeros(n,1);
    Omega = [];
    A_Omega = [];
    r = y;
    iter = 1;
    varepsilon = 1e-3; % step 1

    while true
        [~, max_index] = max(abs(A'*r));
        Omega = [Omega max_index];
        A_Omega = [A_Omega A(:,max_index)]; % step 2
        x_k = A_Omega\y; % step 3
        r = y - A_Omega*x_k; % step 3
        iter = iter + 1; % step 4
        if (iter > opts.maxIter) || (norm(r) <= opts.tol) || (norm(A'*r, inf) <= varepsilon)
            break; % ��ֹ����
        end
    end  
    x(Omega) = x_k; %step 5
    Out.iter = iter;
end