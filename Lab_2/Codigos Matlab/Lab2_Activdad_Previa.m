% Parámetros
Frecuencia_Base = 1; % Frecuencia base
Valores_Alpha = [0, 0.25, 0.75, 1]; % Factores de roll-off
t = linspace(0, 5/Frecuencia_Base, 1000);  % Tiempo
f = linspace(-2*Frecuencia_Base, 2*Frecuencia_Base, 1000); % Frecuencia

% Respuesta al impulso he(t)
for Alpha = Valores_Alpha
    fd = Alpha * Frecuencia_Base;
    he_t = 2*Frecuencia_Base * (sin(2*pi*Frecuencia_Base*t)./(2*pi*Frecuencia_Base*t)) .* (cos(2*pi*fd*t)./(1 - (4*fd*t).^2));
    he_t(t==0) = 2*Frecuencia_Base;  % Valor en t = 0

    figure;
    plot(t, he_t);
    title(['Respuesta al impulso, \alpha = ' num2str(Alpha)]);
    xlabel('Tiempo t');
    ylabel('he(t)');
    grid on;
end

% Respuesta en frecuencia He(f)
for Alpha = Valores_Alpha
    H = zeros(size(f));
    for i = 1:length(f)
        Frecuencia_absoluta = abs(f(i));
        if Frecuencia_absoluta < Frecuencia_Base * (1 - Alpha)
            H(i) = 1;
        elseif Frecuencia_absoluta <= Frecuencia_Base * (1 + Alpha)
            H(i) = 0.5 * (1 + cos(pi/(2*Alpha*Frecuencia_Base) * (Frecuencia_absoluta - Frecuencia_Base * (1 - Alpha))));
        else
            H(i) = 0;
        end
    end

    figure;
    plot(f, H);
    title(['Respuesta en frecuencia, \alpha = ' num2str(Alpha)]);
    xlabel('Frecuencia f');
    ylabel('He(f)');
    grid on;
end
