function kfcs(sigma_mu, sigma_sr, name)
printfigure = 0;


load phi
load basis_ml

            LOAD = ['load ',name];
            eval(LOAD);
            

            x = 0 : dx : 10;
            nx = length(x);

            Dt = 0.1;
            m_rh = m_rh(:,1:Dt/dt:end);
            nt = size(m_rh,2);

            % random select of sensors
            number = 12;
            S = zeros(nt,number);
            S(:,1) = ceil(rand(nt,1)*64);
            for t = 1 : nt
                k = 2;
                while k <= number
                    temp = ceil(rand()*64);
                    if abs(S(t,:) - temp) ~= 0
                        S(t,k) = temp;
                        k = k + 1;
                    end
                end
                S(t,:) = sort(S(t,:));
            end


            % CS Basis
            Phi_cs = zeros(nt, number, 64);
            for t = 1 : nt
                for i = 1 : number
                    Phi_cs(t,i,S(t,i)) = 1;
                end
            end

            z = zeros(64, nt);
            z(:,1) = MT * m_rh(:,1);


            k = 0.1;


            E = sparse(2:nx,1:nx-1,1,nx,nx);
            P = Dt*k/dx^2;
            A = -2*speye(nx)+(E+E'); 
            A(1,1) = -1; % Neumann Boundary
            A(nx,nx) = -1; % Neumann Boundary
            D = speye(nx) - P*A;

            t = 0 : Dt : 20;
            u = zeros(nt, nx);
            u(:,round(3/dx+1)) = 0.1 * sin(t - pi / 4) / dx;
            u(:,round(5/dx+1)) = -0.2 * sin(t) / dx;
            u(:,round(7/dx+1)) = 0.01 * t / dx;


            Phi_KF = zeros(length(p_index),nx);
            for i = 1 : length(p_index)
                Phi_KF(i,p_index(i)) = 1;
            end


            H = zeros(64, nx);
            for i = 1 : 64
                H(i,p_index(i)) = 1;
            end



            f_e_integrated_1 = zeros(nx, nt);
            f_e_integrated_1(:,1) = f_sr(:,1);
            Ce = speye(nx) * 0;
            Cw = phi * eye(51) * phi' * sigma_sr;  % Systemrauschen
            
            for t = 2 : nt
                h = squeeze(Phi_cs(t,:,:)) * THETA;
                y = m_rh(S(t,:),t);

                %dynamic weight
                e = 0.01;
                Cv = eye(length(p_index)) * 1; 
                for i = 1 : size(S,2)
                    Cv(S(t,i),S(t,i)) = e;
                end

                fp = D \ f_e_integrated_1(:, t-1) + D^-1 * u(t-1,:)' * Dt;
                Cp = D \ (Ce + Cw) / D';
                K = Cp * H' / (H * Cp * H' + Cv);


                % iterative
                zi = z(:,t-1);
                for iterative = 1 : 10
                    cvx_begin quiet;
                        variable a(64,1);
                        minimize(norm(a,1));
                        subject to;
                            norm(h * (zi + a) - y) <= e;
                    cvx_end;

                    zi_new = zi + a;
                    ypi = THETA * zi_new;

                    fi = (eye(size(K,1)) - K * H) * fp + K * ypi;
                    zi_kf = MT * Phi_KF * fi;

                    if norm(zi_kf - zi) < 1e10
                        break
                    else
                        zi = zi_kf;
                    end
                end
                z(:, t) = zi_kf;
                f_e_integrated_1(:,t) = fi;

%                 ['t = ',num2str(t),',    iterativ = ',num2str(iterative)]
                Ce = (speye(size(K,1)) - K * H) * Cp;
            end

            f_e_kfcs = f_e_integrated_1;

            name = ['f_e_kfcs_',name];
            save(name,'f_e_kfcs')
            

