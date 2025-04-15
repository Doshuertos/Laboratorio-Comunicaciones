% Parámetros principales
num_bits = 10^4;           % Número de bits
Rs = 1e3;                  % Tasa de símbolos (1 kHz)
sps = 8;                   % Muestras por símbolo
Fs = Rs * sps;             % Frecuencia de muestreo 
rolloff = 0.75;               % Factor de roll-off 
span = 10;                 % Span del filtro en símbolos
snr_dB = 20;               % Relación señal/ruido en dB

% Generación de bits aleatorios
bits = randi([0, 1], 1, num_bits);

% Codificación NRZ-L (0 -> -1, 1 -> +1)
symbols = 2 * bits - 1;

% Sobremuestreo (interpolación)
symbolsUp = upsample(symbols, sps);

% Crear y aplicar filtro Raised Cosine (coseno alzado)
rrcFilter = rcosdesign(rolloff, span, sps, 'normal');
filteredSignal = filter(rrcFilter, 1, symbolsUp);

% Calcular potencia de la señal para ajustar el nivel de ruido
filteredSignalPower = sum(abs(filteredSignal).^2) / length(filteredSignal);
noisePower = filteredSignalPower / (10^(snr_dB / 10));
noise = sqrt(noisePower) * randn(size(filteredSignal));

% Agregar ruido gaussiano blanco (AWGN)
receivedSignal = filteredSignal + noise;

% Crear vector de tiempo para graficar la señal (opcional)
t = (0:length(receivedSignal)-1) / Fs;


% Mostrar diagrama de ojo
figure;
eyediagram(receivedSignal, 2 * sps);  % Dos símbolos por ventana de ojo
title(['Diagrama de ojo para pulso coseno alzado (α = ', num2str(rolloff), ') y 1KHz de Fs']);
xlabel('Tiempo');
ylabel('Amplitud');