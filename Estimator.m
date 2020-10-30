classdef Estimator
    properties
        RisPhaseProfile
        Error_Squared      %collects the results
        config      % contains the config
        randomSeed  % scalar for seed of random number generator
        risElementLoc %contias the location of ris elements
        IV    %stores some initial calculations for channel parameters
        FIMch % fisher information matrix for channel parameters
        FIMpo %fisher information matrix for positional parametrs
        PEB %peb vector
        Jacob
    end
    
    methods
        function obj = Estimator(config,risPhaseProfile,randomSeedin)
            % constructor: input the config file and the Ris phase profile
            obj.config = config;
            obj.randomSeed = randomSeedin;
            obj.RisPhaseProfile = risPhaseProfile; % this is the RIS phase profile
            obj.Error_Squared = obj.errinit;
            obj.risElementLoc = obj.getRisElLoc;
            obj.IV = obj.getInits;
        end
        
        function Err_Squared = errinit(obj)
            % initializing the error matrices to zeros
            Z = zeros(obj.config.NoiseSampleNum,obj.config.PointNum);
            Err_Squared.PEB = Z;
            Err_Squared.tau_b = Z;
            Err_Squared.tau_r = Z;
            Err_Squared.phi_az = Z;
            Err_Squared.phi_el = Z;
            Err_Squared.DT = Z;
        end
        
        function RisElLoc = getRisElLoc(obj)
            % calculates a 3XM matrix containing ris element locations
            Lris = (obj.config.Mc-1)*obj.config.risElementDist;
            risXvec = linspace(0,Lris,obj.config.Mc)-Lris/2;
            risZvec = risXvec;
            [risX, risZ] = meshgrid(risXvec,risZvec);
            RisElLoc= [risX(:),zeros(length(risX(:)),1),risZ(:)];
        end
        
        function IV = getInits(obj)
            % calculates some initial values like channel parameters and
            % Power, noise variance , and poistions
            
            IV.sqEs = sqrt(obj.config.TxPower/obj.config.Nsc); %quare root symbol energy
            IV.sqN = sqrt(obj.config.Noise_Factor*obj.config.NPSD*obj.config.Df/2);%square root noise power
            IV.snr = obj.config.TxPower/(obj.config.Nsc*obj.config.Noise_Factor*obj.config.NPSD*obj.config.Df);
            
            IV.uePos = SupPoint(obj.config.xyz(1,:),obj.config.xyz(2,:),obj.config.xyz(3,:)); %ue positions
            IV.BsPos = SupPoint(obj.config.bsPos); %Bs positon
            IV.BsUeVec = SupPoint(IV.uePos.X-IV.BsPos.X, IV.uePos.Y-IV.BsPos.Y, IV.uePos.Z-IV.BsPos.Z);% vector from Bs to Ue
            IV.F = ifft(eye(obj.config.Nsc),obj.config.N_F); %IFFT matrix
            IV.options = optimoptions('fminunc','Display','off', 'OptimalityTolerance', 1e-10);%use for fminunc
            
            
            IV.Gam = zeros(obj.config.N_F_tilde^2,obj.config.T); %the vector that represents the all ris pahse shifts
            for tc =1:size(obj.RisPhaseProfile,2)
                G = ifft2(reshape(obj.RisPhaseProfile(:,tc),obj.config.Mc,obj.config.Mc),obj.config.N_F_tilde,obj.config.N_F_tilde);
                IV.Gam(:,tc) = G(:);
            end
            
            IV.gb = sqrt((obj.config.lambda./(4*pi*IV.BsUeVec.abs(:))).^2).*exp(2*pi*1j*rand(size(IV.BsUeVec.abs(:)))); % los gain
            IV.gr = sqrt((obj.config.lambda./(4*pi*IV.BsPos.abs(:))).^2 .* (obj.config.lambda./(4*pi*IV.uePos.abs(:))).^2).*exp(2*pi*1j*rand(size(IV.BsUeVec.abs(:)))); % reflected path gatin
            IV.tau_b = IV.BsUeVec.abs(:)./obj.config.c; % los ToA
            IV.tau_r = (IV.BsPos.abs(:)+IV.uePos.abs(:))./obj.config.c; % NLOS TOA
            IV.k = IV.uePos.getKv(obj.config.lambda); % wavenumber vector
            IV.ur = exp(-1j*IV.k*obj.risElementLoc.') * obj.RisPhaseProfile; %  u_r
            
            
            IV.kaz = IV.uePos.getKDazv(obj.config.lambda);  %d k /d \phi_az
            IV.kel = IV.uePos.getKDelv(obj.config.lambda);
            IV.urAz = (exp(-1j*IV.k*obj.risElementLoc.').* (-1j*IV.kaz*obj.risElementLoc.')) * obj.RisPhaseProfile;
            IV.urEl = (exp(-1j*IV.k*obj.risElementLoc.').* (-1j*IV.kel*obj.risElementLoc.')) * obj.RisPhaseProfile;
        end
        
        function obj = estimate(obj)
            % this is the mainfunction that performs the estimation
            rng(obj.randomSeed)
            maxFinder = @(m1,f1, f2, f3) m1+(f2-f3)/(2*(f3+f2-2*f1));% quadrature approx, f1= f(m1), f2= f(m1-1), f3= f(m1+1)
            urfuc = @(phi) exp(-1j*2*pi/obj.config.lambda*[sin(phi(1))*cos(phi(2)),sin(phi(1))*sin(phi(2)),cos(phi(1))]*obj.risElementLoc.')*obj.RisPhaseProfile; %calculates the ur vector for a given angle
            delta_tau = @(t1,t2)  min([abs(t1-t2),abs(t1-t2+1/obj.config.Df),abs(t1-t2-1/obj.config.Df)]); % calculates the time error
            
            tic
            for ix =1:obj.config.PointNum % this for loop goes through all the ue positions
                fprintf('Estimation %.2f %% done! \n', (ix-1)/obj.config.PointNum *100);
                for inr = 1:obj.config.NoiseSampleNum % for different noise realizations
                    
                    Dt = (1/obj.config.Df)*rand;
                    db = exp(-1j*2*pi*mod(obj.config.Df*(obj.IV.tau_b+Dt)*[0:obj.config.Nsc-1],1));
                    dr = exp(-1j*2*pi*mod(obj.config.Df*(obj.IV.tau_r+Dt)*[0:obj.config.Nsc-1],1));
                    % received signal without noise
                    Y_NL = obj.IV.gb(ix)*obj.IV.sqEs*db(ix,:).'*ones(1,obj.config.T) + obj.IV.gr(ix)*obj.IV.sqEs*((dr(ix,:).'*ones(1,obj.config.T)).*(ones(obj.config.Nsc,1)*obj.IV.ur(ix,:)));
                    Y = Y_NL+obj.IV.sqN*(randn(size(Y_NL))+1j*randn(size(Y_NL)));%received signal
                    
                    %% Estimation tau_b
                    y_hat= sum(Y,2);
                    [y0,k_tilde] = max(vecnorm(ifft(y_hat,obj.config.N_F),2,2));
                    y_hat_fun = @(tau) -1*vecnorm(obj.IV.F(k_tilde,:)*(y_hat .* exp(-1j*2*pi*mod(obj.config.Df*tau*[0:obj.config.Nsc-1].',1))))/y0;
                    delta_tilde = fminunc(y_hat_fun,0,obj.IV.options);
                    Est_taub = (k_tilde-1)/(obj.config.Df*obj.config.N_F)-delta_tilde;
                    %% Estimation tau_r
                    % removing LOS signal
                    gHatB = sum(y_hat .* exp(1j*2*pi*mod(obj.config.Df*Est_taub*[0:obj.config.Nsc-1].',1)))/(obj.config.Nsc*obj.config.T);
                    Yt = Y -  gHatB * repmat(exp(-1j*2*pi*mod(obj.config.Df*Est_taub*[0:obj.config.Nsc-1].',1)),1,obj.config.T);
                    %%
                    [y0,k_tilde] =  max(vecnorm(ifft(Yt,obj.config.N_F),2,2));
                    y_hat_fun = @(tau) -1*max(vecnorm(obj.IV.F(k_tilde,:)*(Yt .* repmat(exp(-1j*2*pi*mod(obj.config.Df*tau*[0:obj.config.Nsc-1].',1)),1,obj.config.T)),2,2)/y0);
                    delta_tilde = fminunc(y_hat_fun,0,obj.IV.options);
                    Est_taur = (k_tilde-1)/(obj.config.Df*obj.config.N_F)-delta_tilde;
                    %% estimating AOD
                    Ytt = Yt .* repmat(exp(1j*2*pi*mod(obj.config.Df*Est_taur*[0:obj.config.Nsc-1].',1)),1,obj.config.T);
                    y_phi = sum(Ytt);
                    % IFFT search
                    h_Gam = obj.IV.Gam*y_phi'./(y_phi*y_phi'); % finding the best gain as in
                    Err = vecnorm(h_Gam*y_phi-obj.IV.Gam,2,2).^2;
                    [~,mindx] = min(Err);
                    [l_tilde,m_tilde] = ind2sub([obj.config.N_F_tilde,obj.config.N_F_tilde],mindx);
                    % quadratic refinement
                    ErrS = reshape(Err,obj.config.N_F_tilde,obj.config.N_F_tilde);
                    mrF1 = ErrS(l_tilde,m_tilde);
                    mrF2 = ErrS(mod(l_tilde-2,obj.config.N_F_tilde)+1,m_tilde);
                    mrF3 = ErrS(mod(l_tilde,obj.config.N_F_tilde)+1,m_tilde);
                    mcF1 = ErrS(l_tilde,m_tilde);
                    mcF2 = ErrS(l_tilde,mod(m_tilde-2,obj.config.N_F_tilde)+1);
                    mcF3 = ErrS(l_tilde,mod(m_tilde,obj.config.N_F_tilde)+1);
                    l_hat = maxFinder(l_tilde, mrF1, mrF2, mrF3);
                    m_hat = maxFinder(m_tilde, mcF1, mcF2, mcF3);
                    
                    % calculating (normalized)wavenumbre vector based on  l_hat and m_hat
                    k3 = - obj.config.lambda/obj.config.risElementDist*(l_hat-1)/obj.config.N_F_tilde;
                    if k3<-1 % this is because of 2*pi amibiguity
                        k3=2+k3;
                    end
                    
                    k1 =  - obj.config.lambda/obj.config.risElementDist*(m_hat-1)/obj.config.N_F_tilde;
                    if k1<-1
                        k1=2+k1;
                    end
                    if(k1^2+k3^2>1)
                        dk = sqrt(1/(k1^2+k3^2));
                        k1=k1*dk;
                        k3=k3*dk;
                    end
                    k_bar = [k1,real(sqrt(1-k1^2-k3^2)),k3];
                    
                    %% refinement based on quasi-newton algorithm
                    squ = y_phi*y_phi';
                    y0 = vecnorm((urfuc(k_bar)*y_phi'/squ)*y_phi-urfuc(k_bar)); %normalization factor
                    likelihood_phi= @(phi)vecnorm((urfuc(phi)*y_phi'/squ)*y_phi-urfuc(phi))/y0;
                    phi_tilde = fminunc(likelihood_phi,[acos(k_bar(3)),atan2(k_bar(2),k_bar(1))],obj.IV.options);
                    k_bar = [sin(phi_tilde(1))*cos(phi_tilde(2)),sin(phi_tilde(1))*sin(phi_tilde(2)),cos(phi_tilde(1))];
                    
                    %%
                    DiffTau = min(abs(Est_taur-Est_taub),abs(Est_taur+(1/obj.config.Df)-Est_taub));%to cover for an wrap arround becoause of Dt\approx 1/DF
                    d_fun = @(d) abs(vecnorm(k_bar*d)-vecnorm(k_bar*d-obj.IV.BsPos.getV(1,1))+obj.IV.BsPos.abs-DiffTau*obj.config.c);
                    d_tilde = fminunc(d_fun,1,obj.IV.options);
                    Est_position = d_tilde*k_bar; % approx position
                    
                    %% calculating Errors
                    Est_Dt = Est_taub- vecnorm(Est_position-obj.IV.BsPos.getV(1,1))/obj.config.c;
                    Etb = delta_tau(Est_Dt, Dt);
                    E = Est_position-obj.IV.uePos.getV(1,ix);
                    obj.Error_Squared.PEB(inr,ix) = E*E';
                    obj.Error_Squared.DT(inr,ix) =  Etb^2;
                    obj.Error_Squared.tau_b(inr,ix) =  delta_tau(obj.IV.tau_b(ix)+Dt,Est_taub)^2;
                    obj.Error_Squared.tau_r(inr,ix) =  delta_tau(obj.IV.tau_r(ix)+Dt,Est_taur)^2;
                    obj.Error_Squared.phi_az(inr,ix) =  (atan2(Est_position(2),Est_position(1))-atan2(obj.IV.k(ix,2),obj.IV.k(ix,1)))^2;
                    obj.Error_Squared.phi_el(inr,ix) =  (acos(Est_position(3)/vecnorm(Est_position))-acos(obj.IV.k(ix,3)/vecnorm(obj.IV.k(ix,:))))^2;
                end
                
                timePassed=toc;
                fprintf('Estimation finishes about %.2f minutes! \n', timePassed/ix*(obj.config.PointNum-ix)/60);
                
            end
            fprintf('Total time = %.2f minutes', timePassed/60)
        end
        
        function Jc = getJacobian(obj)
            % returens a cell of the length of number of points each with a 8X8 jacobian matrix
            J = cell(8,8);
            J{1,1} = 1/obj.config.c*(obj.IV.uePos.X-obj.IV.BsPos.X)./(obj.IV.BsUeVec.abs); %d\tau_{b}/d\p1
            J{1,2} = 1/obj.config.c*(obj.IV.uePos.Y-obj.IV.BsPos.Y)./(obj.IV.BsUeVec.abs); %d\tau_{b}/d\p2
            J{1,3} = 1/obj.config.c*(obj.IV.uePos.Z-obj.IV.BsPos.Z)./(obj.IV.BsUeVec.abs); %d\tau_{b}/d\p3
            J{2,1} = 1/obj.config.c*(obj.IV.uePos.X)./(obj.IV.uePos.abs); %d\tau_{r}/d\p1
            J{2,2} = 1/obj.config.c*(obj.IV.uePos.Y)./(obj.IV.uePos.abs); %d\tau_{r}/d\p2
            J{2,3} = 1/obj.config.c*(obj.IV.uePos.Z)./(obj.IV.uePos.abs); %d\tau_{r}/d\p3
            J{3,1} = -obj.IV.uePos.Y./(obj.IV.uePos.XY2); %d\phi_{az}/d\p1
            J{3,2} = obj.IV.uePos.X./(obj.IV.uePos.XY2); %d\phi_{az}/d\p1
            J{4,1} = obj.IV.uePos.X.*obj.IV.uePos.Z./(obj.IV.uePos.abs2.*obj.IV.uePos.XY); %d\phi_{el}/d\p1
            J{4,2} = obj.IV.uePos.Y.*obj.IV.uePos.Z./(obj.IV.uePos.abs2.*obj.IV.uePos.XY); %d\phi_{el}/d\p2
            J{4,3} = (-obj.IV.uePos.abs2 + obj.IV.uePos.Z.*obj.IV.uePos.Z)./(obj.IV.uePos.abs2.*obj.IV.uePos.XY); %d\phi_{el}/d\p3
            Jc = cell(size(obj.IV.uePos.X));
            for ir=1:size(obj.IV.uePos.X,1)
                for ic=1:size(obj.IV.uePos.X,2)
                    Jc{ir,ic}=[J{1,1}(ir,ic),J{1,2}(ir,ic),J{1,3}(ir,ic),-1,0,0,0,0;...
                        J{2,1}(ir,ic),J{2,2}(ir,ic),J{2,3}(ir,ic),-1,0,0,0,0;...
                        J{3,1}(ir,ic),J{3,2}(ir,ic),0,0,0,0,0,0;...
                        J{4,1}(ir,ic),J{4,2}(ir,ic),J{4,3}(ir,ic),0,0,0,0,0;...
                        0,0,0,0,1,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1];
                end
            end
        end
        
        function FIMch = getFIMCh(obj)
            db = exp(-1j*2*pi*mod(obj.config.Df*obj.IV.tau_b*[0:obj.config.Nsc-1],1));
            dr = exp(-1j*2*pi*mod(obj.config.Df*obj.IV.tau_r*[0:obj.config.Nsc-1],1));
            d = -1j*2*pi*[0:obj.config.Nsc-1]*obj.config.Df; %this is to calculate the extra factors in d/d\tau terms
            FIMch = cell(1,obj.config.PointNum);
            for ip = 1:obj.config.PointNum
                fprintf('Calculating CRB %.2f %% done! \n', (ip-1)/obj.config.PointNum*100);
                FIMch{1,ip} = zeros(8);
                for it = 1:obj.config.T
                    for in = 1:obj.config.Nsc
                        u = [
                            d(in)*obj.IV.gb(ip)*db(ip,in),...
                            d(in)*obj.IV.gr(ip)*dr(ip,in)*obj.IV.ur(ip,it),...
                            obj.IV.gr(ip)*dr(ip,in)*obj.IV.urAz(ip,it),...
                            obj.IV.gr(ip)*dr(ip,in)*obj.IV.urEl(ip,it),...
                            db(ip,in),...
                            1j*db(ip,in),...
                            dr(ip,in)*obj.IV.ur(ip,it),...
                            1j*dr(ip,in)*obj.IV.ur(ip,it),...
                            ];
                        FIMch{1,ip} = FIMch{1,ip} + 2*obj.IV.snr* real(u'*u);
                    end
                end
            end
        end
        
        function obj = PEBcalc(obj)
            %calculates fisher information for positional parameters and
            %PEB
            obj.Jacob = obj.getJacobian;
            obj.FIMch = obj.getFIMCh;
            obj.FIMpo = cell(1,obj.config.PointNum);
            obj.PEB = zeros(1,obj.config.PointNum);
            for ip = 1:obj.config.PointNum
                Jcb = obj.Jacob{1,ip};
                obj.FIMpo{1,ip} = Jcb.' * obj.FIMch{1,ip} * Jcb;
                
                fim = obj.FIMpo{1,ip};% using Schur’s complement
                A = fim(1:4,1:4);
                B = fim(1:4,5:8);
                C = fim(5:8,1:4);
                D = fim(5:8,5:8);
                eqFim = A-B*inv(D)*C;
                FimInv = inv(eqFim);
                obj.PEB(1,ip) = sqrt(trace(FimInv(1:3,1:3)));
            end
        end
        
        function crb = getCrb(obj,FIM,parNum)
            % calculatest the crb bound for a given FIM (which is a cell) and a
            % parameter,
            crb=zeros(size(FIM));
            for ic = 1:length(FIM)
                fim=FIM{ic};% using Schur’s complement
                A = fim(1:4,1:4);
                B = fim(1:4,5:8);
                C = fim(5:8,1:4);
                D = fim(5:8,5:8);
                eqFim = A-B*inv(D)*C;
                FimInv = inv(eqFim);
               crb(ic) = sqrt(FimInv(parNum,parNum)); 
            end
        end
        
    end
end


