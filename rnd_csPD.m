%%%%%%%%%%%%%%%%%%%%%%%
function A = rnd_csPD(n,imagScale)
% Devuelve A compleja simetrica, no hermitiana, definida positiva
% imagScale controla la “fuerza” de la parte imaginaria (≈0.1 -- 0.5)

if nargin<2, imagScale = 0.3; end

% --- Parte real (SPD) ---
C = randn(n);                 % aleatoria real
R = C.'*C + n*eye(n);         % SPD y con autovalores ≥ n

% --- Parte imaginaria (simetrica, quizá grande) ---
S0 = randn(n);  S = (S0+S0.')/2;      % simetrica real
S  = imagScale * S / norm(S,2);       % control de magnitud

A = R + 1i*S;                 % matriz resultante
end