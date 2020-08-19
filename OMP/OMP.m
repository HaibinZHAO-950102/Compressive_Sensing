function [x, Out] = OMP(A,y,varargin)
% This is a l1 minimization solver for orthogonal matching pursuit: 
%   minimize     ||x||_1
%   subject to   y = Ax,
% -----------------------------------------------------------------
% Author: Beier ZHU
%         Tsinghua University
% -----------------------------------------------------------------
% 
% =========================Required inputs=========================  
% A -- an m x n matrix
% y -- an m x 1 vector
% =========================Optional inputs=========================
% 'maxIter' -- maximum number of iterations
% 'StopTolerance' -- stopping tolerance
% ===========================Outputs===============================
% x -- last iterate (hopefully an approximate solution)
% Out.iter -- # of iterations taken

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
            break; % ÖÕÖ¹Ìõ¼þ
        end
    end  
    x(Omega) = x_k; %step 5
    Out.iter = iter;
end