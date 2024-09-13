function [parameters, expectedinfomatrix, lvalue] = recoverParameter(data, years)
    disp('alpha, r, sigma'); % Mostra l'ordine dei parametri
    n = length(data) - 1;
    % 'n' corrisponde al numero di cambiamenti, dato che è discreto.
    % Nota che abbiamo sempre n+1 punti nel nostro campione.
    dt = years / (n + 1); % Dimensione di ogni "passo".
    parameters = zeros(1, 3);
    expectedinfomatrix = zeros(3, 3);
    S0 = 0; % Inizializzazione per la somma di tutti R(t(i-1))
    S1 = 0; % Inizializzazione per la somma di tutti R(t(i))
    S00 = 0; % Inizializzazione per la somma di tutti R(t(i-1))^2
    S01 = 0; % Inizializzazione per la somma di tutti R(t(i-1))*R(t(i))
    for i = 1:n % Calcola ogni somma
        S0 = S0 + data(i);
        S1 = S1 + data(i + 1);
        S00 = S00 + data(i) * data(i);
        S01 = S01 + data(i) * data(i + 1);
    end
    ahat = -(1 / dt) * log((n * S01 - S0 * S1) / (n * S00 - S0^2)); % Assegna un valore ad alpha-hat
    rhat = 1 / (n * (1 - exp(-ahat * dt))) * (S1 - exp(-ahat * dt) * S0); % Assegna un valore a r-hat
    sigmapart = 0; % Inizializzazione per la somma nell'estimatore di sigma
    for j = 1:n % Calcola la somma nell'estimatore di sigma
        sigmapart = sigmapart + (data(j + 1) - data(j) * exp(-ahat * dt) - rhat * (1 - exp(-ahat * dt)));
    end
    sigma2 = 2 * ahat * sigmapart / (n * (1 - exp(-2 * ahat * dt))); % Assegna un valore a sigma2-hat
    parameters(1) = ahat;
    parameters(2) = rhat;
    parameters(3) = sqrt(sigma2);
    A = exp(-ahat * dt); % Per semplificare (vedi dimostrazione per teorema 3.2 e 3.3)
    C = sigma2 / (2 * ahat) * (1 - exp(-2 * ahat * dt)); % Per semplificare (vedi dimostrazione)
    dCdalpha = -(sigma2 * A^2 * (-2 * ahat * dt + exp(2 * ahat * dt) - 1)) / (2 * ahat^2); % Derivata per semplificare (vedi dimostrazione)
    d2La2 = -1 / C * S00 * dt^2 * A^2; % Derivata per semplificare (vedi dimostrazione)
    d2Lb2 = -n / C * (rhat * dt * A)^2; % Derivata per semplificare (vedi dimostrazione)
    d2Lc2 = -dCdalpha * n / (2 * C^2); % Derivata per semplificare (vedi dimostrazione)
    d2Lab = rhat * (dt * A)^2 * S0 / C; % Derivata per semplificare (vedi dimostrazione)
    d2Lalpha2 = -(d2La2 + d2Lb2 + d2Lc2 + 2 * d2Lab); % Seconda derivata completa (vedi dimostrazione)
    d2Lr2 = (2 * n * ahat * (1 - exp(-ahat * dt))^2) / (sigma2 * (1 - exp(-2 * ahat * dt))); % Seconda derivata completa (vedi dimostrazione)
    d2Lsigma22 = (2 * n) / sigma2; % Seconda derivata completa (vedi dimostrazione)
    % Le seguenti sette linee completano la matrice di informazione negativa aspettata
    expectedinfomatrix(1, 1) = d2Lalpha2;
    expectedinfomatrix(2, 2) = d2Lr2;
    expectedinfomatrix(3, 3) = d2Lsigma22;
    expectedinfomatrix(2, 1) = -2 * ahat * dt * (S0 - n * rhat) / (sigma2 * (1 + exp(ahat * dt)));
    expectedinfomatrix(1, 2) = expectedinfomatrix(2, 1);
    expectedinfomatrix(3, 1) = -n / (sqrt(sigma2) * ahat) * (1 - exp(-2 * ahat * dt) * (2 * ahat * dt + 1)) / (1 - exp(-2 * ahat * dt));
    expectedinfomatrix(1, 3) = expectedinfomatrix(3, 1);
    lvalue = -n / 2 * log(sigma2 / (2 * ahat) * (1 - exp(-2 * ahat * dt))) - n / 2 * log(2 * pi) - (sigma2 / ahat * (1 - exp(-2 * ahat * dt)))^(-1) * sigmapart;
    % Assegna un valore al log-likelihood basato sugli stimatori
end

