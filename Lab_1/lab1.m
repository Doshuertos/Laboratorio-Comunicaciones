clc
clear all
close all

% Parámetros iniciales
fc = 1000; % frecuencia de la señal sinosoidal
fm = 100000; % frecuencia de muestreo 
tm = 1/fm; % periodo de muestreo
tiempo_max = 0.01;%aqui va un 1 segun guia pero va 0.01 para que se vea bien
ls = 200; % numero de muestras son 200

% Generar vector de tiempo
Tiempo = (0:ls-1)*tm;

% Señal sinosoidal (portadora)
y = sin(2*pi*fc*Tiempo);

% Modulación PAM
fs = 5000; % frecuencia de muestreo
ts = 1/fs; % periodo de muestreo
tau = 0.5*ts; % duracion del pulso
d = tau/ts; % ciclo del trabajo

r = floor(ts/tm); % Número de muestras entre dos instantes de muestreo
s = floor(tau/tm); % Número de muestras que dura el pulso activo
disp(r)
% ----- Señal natural ------  PAM MUESTRO NATURAL
Vector_natural_muestral = zeros(1,length(Tiempo)); %Vector generado de 0
for i=1:length(y) % se genera un "tren de pulsos"
    if mod(i,r)==0
        Vector_natural_muestral(i:i+s) = 1;
    end
end
Vector_natural_muestral = Vector_natural_muestral(1:length(Tiempo));
Muestreo_Natural = y.*Vector_natural_muestral;

% ----- Señal Instantanea ------ PAM MUESTRO INSTANTANEO
Vector_instantaneo_muestral = zeros(1,length(Tiempo));
for i=1:length(y)
    if mod(i,r)==0
       Vector_instantaneo_muestral(i:i+s) = y(i);
    end
end
Muestreo_Instantaneo = Vector_instantaneo_muestral(1:length(Tiempo));
%--------- Transformadas ---------

%fft Y original
Transformada_Sinusoydal=fft(y);
Transformada_Muestreo_Natural = fft(Muestreo_Natural);
Transformada_Muestreo_Instantaneo = fft(Muestreo_Instantaneo);
Tamano_Transformada = length(y); 
f = linspace(0, fm/2, Tamano_Transformada/2+1); % Eje de frecuencias

% Calcular la magnitud del espectro (solo parte positiva)
Espectro_Sinusoidal = abs(Transformada_Sinusoydal(1:Tamano_Transformada/2+1));
Espectro_Muestreo_Natural = abs(Transformada_Muestreo_Natural(1:Tamano_Transformada/2+1));
Espectro_Muestreo_Instantaneo = abs(Transformada_Muestreo_Instantaneo(1:Tamano_Transformada/2+1));

%---------- PCM -------------
% Parámetros PCM
N = 200; % Número de bits para PCM
Amplitud_maxima  = max(Muestreo_Instantaneo);
Amplitud_minima  = min(Muestreo_Instantaneo);
Cuantizacion = (Amplitud_maxima - Amplitud_minima) / (2^N - 1);
PCM_Cuantizado = round((Muestreo_Instantaneo-Amplitud_minima)/Cuantizacion) * Cuantizacion + Amplitud_minima; % Cuantizar la señal
Error_Cuantizacion = Muestreo_Instantaneo - PCM_Cuantizado;% calcular el error
erro_medio_de_Cuantizacion = mean(abs(Error_Cuantizacion)); %Mean : devuelve un valor numerico de promedio , abs: valor abs
disp(['Error Cuantizacion :', num2str(erro_medio_de_Cuantizacion)]); % num2str : devuelve un areglo numerico en un areglo de caracteres
% ------ Graficas ---------

% --------- PAM, PAM NAT, PAM INST ---------

figure;
subplot(4,1,1);
plot(Tiempo, y);
title('Señal de Pulso');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% y_pam_natural 
subplot(4,1,2);
plot(Tiempo, Muestreo_Natural);
title('PAM Muestreo Natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% y_pam_instantaneo
subplot(4,1,3);
plot(Tiempo, Muestreo_Instantaneo);
title('PAM Muestreo Instantáneo');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% ----- Pam inst , Pam nat y pam juntas ----------

figure;
plot(Tiempo, y, 'y', 'LineWidth', 1.5); 
hold on;
plot(Tiempo, Muestreo_Natural, 'c', 'LineWidth', 1.5);
plot(Tiempo, Muestreo_Instantaneo, 'm', 'LineWidth', 1.5);
hold off;
title('Superposición de Señales');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('Señal Original', 'PAM Natural', 'PAM Instantáneo');
grid on;

%---------- FFT ------------

figure;
subplot(3,1,1);
plot(f, Espectro_Sinusoidal);
title('Espectro de la Señal Senoidal');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

subplot(3,1,2);
plot(f, Espectro_Muestreo_Natural);
title('Espectro del Muestreo Natural');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

subplot(3,1,3);
plot(f, Espectro_Muestreo_Instantaneo);
title('Espectro del Muestreo Instantáneo');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;
%--------- FFT JUNTAS ---------
figure;
plot(f, Espectro_Sinusoidal, 'b', 'LineWidth', 1.5); hold on;
plot(f, Espectro_Muestreo_Natural, 'r', 'LineWidth', 1.5);
plot(f, Espectro_Muestreo_Instantaneo, 'g', 'LineWidth', 1.5);
hold off;

title('Comparación de Transformadas de Fourier');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
legend('Original', 'Muestreo Natural', 'Muestreo Instantáneo');
grid on;


%------- PCM ----------

figure;
plot(Tiempo, y, 'b', 'LineWidth', 1.5);
hold on;
plot(Tiempo, Muestreo_Instantaneo, 'r', 'LineWidth', 1.5);
stem(Tiempo, PCM_Cuantizado, 'm', 'Marker', 'o', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Original, Señal PAM Instantánea y Señal PAM Cuantificada (PCM)');
legend('Señal Original', 'Señal PAM Instantánea', 'Señal PAM Cuantificada (PCM)');
grid on;
figure;
plot(Tiempo, Error_Cuantizacion, 'k--', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Error de Cuantización');
title('Error de Cuantización para la Señal PAM Cuantificada (PCM)');
grid on;
