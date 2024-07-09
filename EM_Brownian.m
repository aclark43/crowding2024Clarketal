%1/18/15 MR, modified by JI
% Simulate eye movements as a Brownian process
%
%   D: standard diffusion coefficient (the double of Kuang et al, 2012)
%   fc:  sampling frequency (Hertz)
%   N:  number of simulated samples in the eye movement trace
%   nTraces: number of traces to generate, added by Janis
%   biasVars: 3-vector [angle (degrees), speed, speed variance]
%
%   Xe, Ye:  eye movement sequence

function [Xe, Ye]  = EM_Brownian(D,fc,N, nTraces, biasVars)
if ~exist('nTraces', 'var')
    nTraces = 1; 
end

% simplified simulation
K = sqrt(2*D/fc); % standard deviation of jump size at each step % JI: added the 2* so that <r^2> = 4Dt
Xe = K*randn(nTraces,N); 
Ye = K*randn(nTraces,N);

%% biased version added by Janis, based on Murat Aytekin's writeup
if exist('biasVars', 'var')
    % basically what this does is change the mean location of the gaussian
    % at each time based on a bias speed and direction. In the unbiased
    % case, the mean is 0. To estimate the diffusion coefficient in the
    % biased case, just remove the mean displacement at each time, then use
    % one of the methods to estimate the diffusion coefficients in the
    % unbiased case.
    theta = biasVars(1);
    mu = biasVars(2);
    vr = biasVars(3);
    
    biasSpeedX = cosd(theta) * (mu + sqrt(vr) * randn(nTraces, N));
    biasSpeedY = sind(theta) * (mu + sqrt(vr) * randn(nTraces, N));
    
    Xe = Xe + biasSpeedX / fc;
    Ye = Ye + biasSpeedY / fc;
end

Xe = cumsum(Xe, 2);
Ye = cumsum(Ye, 2);

%% old, complicated version?
%{
% more complex simulation
 dt = 1/fc;
 Xeye = zeros (2, N);

 %  Stepsize is normal. I added 2 here to make the D come right
 s = 2* sqrt ( D * dt ) * randn ( 1, N - 1 );

 %  Direction is random.
 a = randn ( 2, N - 1 );
 v = s ./ sqrt ( sum ( a.^2 ) );
 b = spdiags ( v', 0, N-1, N-1 );
 dx(1:2,1:N-1) = a * b;

%  Each position is the sum of the previous steps.
 Xeye(:,2:N) = cumsum ( dx(:,1:N-1), 2 );

 Xe = Xeye(1,:)';   Ye = Xeye(2,:)';
%}
end 


% % % %1/18/15 MR
% % % % Simulate eye movements as a Brownian process
% % % %
% % % %   D: standard diffusion coefficient (the double of Kuang et al, 2012)
% % % %   fc:  sampling frequency (Hertz)
% % % %   N:  number of simulated samples in the eye movement trace
% % % %
% % % %   Xe, Ye:  eye movement sequence
% % % % this code was taken from
% % % % : http://people.sc.fsu.edu/~jburkardt/m_src/brownian_motion_simulation/brownian_motion_simulation.m
% % % 
% % % % Bias = [Angle Amplitude StDeviation];
% % % % example
% % % % Angle = 25; % in deg
% % % % Amplitude = 10; % in arcmin (n deg over the entire duration of the drift period)
% % % % StDeviation = 3; % standard deviation in arcmin around the angular bias (smaller std tighter distribution of values around the bias direction)
% % % 
% % % function [Xe Ye]  = EM_Brownian(D,fc,N, addBias,Bias)
% % % 
% % % % older version
% % % % simplified simulation
% % % % bias is a vector in polar coordinates (amplitude and direction)
% % % % Bias = [angle amp] % angle in deg  and amp in arcmin/s
% % % 
% % % % Dx = 10;
% % % % Dy = 10;
% % % % R = [cosd(0) -sind(0); sind(0) cosd(0)];
% % % % % original simulation in cartesian coordinates
% % % % K = sqrt((2*Dx)/fc);
% % % % Xe = cumsum(K*randn(1,N)+0.02); 
% % % % K = sqrt((2*Dy)/fc);
% % % % Ye = cumsum(K*randn(1,N)+0.02); 
% % % % Xe = R*[Xe; Ye];
% % % % Ye = Xe(2,:);
% % % % Xe = Xe(1,:);
% % % 
% % % % Random walk in polar coordinates
% % %  
% % %  dt = 1/fc;
% % %  Xeye = zeros (2, N);
% % %  Yeye = zeros (2, N);
% % % 
% % %  %  Stepsize is normal. I added 4 here to make the D come right
% % %  % the variance of the random walk is given by dispsq = 2*2D*dt
% % %  s = sqrt ( 4* D * dt ) * randn ( 1, N - 1 );
% % %  
% % %  %  Direction is random.
% % %  if addBias
% % %      ang = normrnd(deg2rad(Bias(1)),Bias(3), 2,N-1);
% % %  else
% % %      ang = randn ( 2, N - 1 );
% % %  end
% % %  v = s ./ sqrt ( sum ( ang.^2 ) );
% % %  
% % %  b = repmat(v,2,1);
% % % 
% % %  dx(1:2,1:N-1) = (ang.* b);
% % %  if addBias
% % %       % add a drift
% % %      [bx by] = pol2cart(deg2rad(Bias(1)),Bias(2)/fc);
% % %      dx = dx + repmat([bx; by], 1, N-1);
% % %  end
% % % %scatter(dx(1,:),dx(2,:))
% % % 
% % % %  Each position is the sum of the previous steps.
% % %  Xeye(:,2:N) = cumsum ( dx(:,1:N-1), 2 );
% % % 
% % %  Xe = Xeye(1,:);   Ye = Xeye(2,:);
% % % 
% % %  