function config = GetConfig()
%returns a structure of the configuration
config.secFactor = 1e9; % time scaling factor 1e9 means that time is in 
%nano seconds and frequency in GHz
config.c = 3e8/config.secFactor; %light speed m/ns
%% ofdm
config.f = 30e9/config.secFactor; %frequency Ghz
config.Df= 120e3/config.secFactor;%bandwidth Ghz
config.Nsc = 3e3;% Subcarrier number
config.T = 256;% Subcarrier number
config.lambda = config.c/config.f;
%% geometry
config.bsPos = [5,5,0]; % Bs position
config.risPos = [0,0,0]; % RIS position
config.Mc = 64; % number of RIS elements in a row
config.risElementDist = .5*config.lambda; % RIS element distance
config.TxPower = 100e-3; % transmit power
config.NPSD = 1e-3*10^(-174/10)*config.secFactor; % noise spectral density
config.Noise_Factor = 10^(8/10); % noise factor
% ------------------------
config.NoiseSampleNum = 100; % number of iterations (this number is 5000 in the paper)
config.N_F = 2^12; % IFFT length for delay estimation
config.N_F_tilde=2^8; % IFFT length for angle estimation
config.PointNum=30; % number of points that is  considered in FIG 1
config.r_vec = 10.^(linspace(0,log10(35),config.PointNum));
config.xyz = [-config.r_vec/sqrt(2);config.r_vec/sqrt(2);-10*ones(1,length(config.r_vec))];%contians xyz corrdinates of the points

end

