function  [S] = kalman_adaptive_2(xp, z, S, dt, m)
% kalman_adaptive_2: Adaptive Kalman filter algorithm for NaveGo INS/GPS
% system, Mohamed version.
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
%			Mohamed, 1999.
%
%           Dan Simon. Optimal State Estimation. Chapter 5. John Wiley
% & Sons. 2006.
%
% Version: 004
% Date:    2017/05/10
% Author:  Rodrigo Gonzalez <rodralez@frm.utn.edu.ar>
% URL:     https://github.com/rodralez/navego

persistent X
persistent D
persistent i

[pp,qp] = size(S.Pp);
PS = pp*qp;

[pr,qr] = size(S.R);
RS = pr*qr;

if isempty(D)
    
    D = zeros(m,RS);
    X = zeros(m,PS);
    i = 1;
    
    S.Qd = (S.G * S.Q * S.G') .* dt;
    S.C = S.R;
    
    S.xi = xp;
    S.Pi = S.Pp;
end


%% MOHAMED 1999

% Innovation covariance
dz = z - (S.H * S.xi);
% d = diag(diag(dz * dz'));
d = (dz * dz');
D(i, :) = reshape(d, 1, RS);
D_m = mean(D);
Di = reshape(D_m, pr, qr);

S.R = Di + S.H * S.Pi * S.H'; % Mohamed (1999), Eq. 18
S.R = diag(diag(S.R));
S.alpha = diag(S.R);
S.C = -Di;

%% DING 2007

labmda = abs( trace(Di - S.R) / trace(S.H * S.Pi * S.H') ); 

S.Qd = (S.Qd * sqrt(labmda));

%% 

I = eye(max(size(S.F)));
S.xp = xp;

% Discretization of continous-time system
S.A =  expm(S.F * dt);          % "Exact" expression
% S.A = I + (S.F * dt);         % Approximated expression
% S.Qd = (S.G * S.Q * S.G') .* dt;

% Step 1, update the a priori covariance matrix Pi
S.Pi = (S.A * S.Pp * S.A') + S.Qd;

% Step 2, update Kalman gain
% S.C = (S.R + S.H * S.Pi * S.H');
S.K = (S.Pi * S.H') / (S.C) ;

% Step 3, update the a posteriori state xp
S.xi = S.A * S.xp;
S.xp = S.xi + S.K * (z - S.H * S.xi);

Pp_old = S.Pp;

% Step 4, update the a posteriori covariance matrix Pp
J = (I - S.K * S.H);
S.Pp = J * S.Pi * J' + S.K * S.R * S.K';    % Joseph stabilized version
% S.Pp = (I - S.K * S.H) * S.Pi ;           % Alternative implementation
S.Pp =  0.5 .* (S.Pp + S.Pp');

%% MOHAMED 1999
% Qk Update

% dx = S.xp - S.xi;
% x = dx * dx';
% X(i, :) = reshape(x, 1, PS);
% X_m = mean(X);
% Xi = reshape(X_m, pp, qp);

% S.Qd = (Xi + S.Pp - (S.A * Pp_old * S.A'));   % Mohamed (1999), Eq. 23
% S.Qd = S.K * Di * S.K' * dt;                  % Mohamed (1999), Eq. 5

%%

i = i + 1;
if (i > m)
    i = 1;
end

end
