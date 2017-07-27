function  [S] = kalman_adaptive(xp, z, S, dt, m)
% kalman_adaptive: Adaptive Kalman filter algorithm for NaveGo INS/GPS system, 
% Ding version.
%
% INPUT:
%  xp, 21x1 a posteriori state vector (old).
%   z, 6x1 innovations vector.
%  dt, time period. 
%   S, data structure with at least the following fields:
%       F,  21x21 state transition matrix.
%       H,   6x21 observation matrix.
%       Q,  12x12 process noise covariance.
%       R,   6x6  observation noise covariance.
%       Pp, 21x21 a posteriori error covariance.
%       G,  21x12 control-input matrix.      
%
% OUTPUT:
%    S, the following fields are updated:
%       xi, 21x1 a priori state vector (new).
%       xp, 21x1 a posteriori state vector (new).
%       A,  21x21 state transition matrix.
%       K,  21x6  Kalman gain matrix.
%       Qd, 21x6  discrete process noise covariance.
%       Pi, 21x21 a priori error covariance.
%       Pp, 21x21 a posteriori error covariance.  
%       C,   6x6  innovation (or residual) covariance.
%
%   Copyright (C) 2014, Rodrigo Gonzalez, all rights reserved.
%
%   This file is part of NaveGo, an open-source MATLAB toolbox for
%   simulation of integrated navigation systems.
%
%   NaveGo is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License (LGPL)
%   version 3 as published by the Free Software Foundation.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this program. If not, see
%   <http://www.gnu.org/licenses/>.
%
% Reference:
%			Ding, 2007.
%
%           Dan Simon. Optimal State Estimation. Chapter 5. John Wiley 
% & Sons. 2006.   
%
% Version: 004
% Date:    2017/05/10
% Author:  Rodrigo Gonzalez <rodralez@frm.utn.edu.ar>
% URL:     https://github.com/rodralez/navego

persistent D
persistent i

[p,q] = size(S.R);
RS = p*q;

if isempty(D)    

    D = zeros(m,RS);
    i = 1;
    S.alpha = 1;
    S.Qd = (S.G * S.Q * S.G') .* dt;
    S.xi = xp;
    S.Pi = S.Pp;
end

I = eye(max(size(S.F)));
S.xp = xp;

% Discretization of continous-time system
S.A =  expm(S.F * dt);          % "Exact" expression
% S.A = I + (S.F * dt);         % Approximated expression
% S.Qd = (S.G * S.Q * S.G') .* dt;

%%

% Step 1, update the a priori covariance matrix Pi
S.Pi = (S.A * S.Pp * S.A') + S.Qd;

% Step 2, update Kalman gain
S.C = (S.R + S.H * S.Pi * S.H');
S.K = (S.Pi * S.H') / (S.C) ;

% Step 3, update the a posteriori state xp
S.xi = S.A * S.xp;
S.xp = S.xi + S.K * (z - S.H * S.xi);

i = i + 1;
if (i > m)
    i = 1;
end

% Step 4, update the a posteriori covariance matrix Pp
J = (I - S.K * S.H);
S.Pp = J * S.Pi * J' + S.K * S.R * S.K';    % Joseph stabilized version     
% S.Pp = (I - S.K * S.H) * S.Pi ;           % Alternative implementation
S.Pp =  0.5 .* (S.Pp + S.Pp');


%% DING 2007

% Innovation covariance
d = z - (S.H * S.xi);
IN = d * d';
D(i, :) = reshape(IN, 1, RS);
IC_m = mean(D);
IC = diag(diag(reshape(IC_m, p, q)));

S.alpha = abs( trace(IC - S.R) / trace(S.H * S.Pi * S.H') ); 

S.Qd = (S.Qd * sqrt(S.alpha));
end
