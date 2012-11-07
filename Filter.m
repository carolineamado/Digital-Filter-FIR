%FIR COM FASE LINEAR - EQUAÇÃO T4 - BANDPASS

clear all
close all

%Ordem do filtro (Ímpar, valor dado)
L = 31; 

%Frequência de corte inferior (Valor Dado)
fcl = 7e3;

%Frequência de corte superior (Valor Dado)
fc2 = 28e3;

%Frequência de amostragem (Valor Dado)
Fs = 70e3;

%Frequência de corte normalizadas, pois a função firls exige isso
freq2 = 2*fc2/Fs;
freq1 = 2*fcl/Fs;

%Vetor com as frequências de corte
freq = [0 freq1/2 freq1 freq2 (1-freq1/2) 1];

%Vetor Amplitudes para interpolarização, requerido pela função do matlab
amp = [0 0 1 1 0 0];


% Calculo dos coeficientes de h
h = firls(L,freq,amp,'h');

%Grafico dos coeficientes de h discretos
figure(1)
subplot(2,1,1)
stem(0:L,h)
title('Coeficientes (simetria impar)');
xlabel('índice n'); ylabel('valor do coeficiente');

%Intervalo omega discreto que serão os pontos utilizado no gráfico
omegadis = 0:0.01:pi;
for l = 1:length(omegadis)
   for n= 0:(L-1)/2
    valores(n+1) = 2*h(n+1)*sin(omegadis(l)*(n - (L/2)));
    end
    soma = sum(valores);
    H(l)= exp(-i*((omegadis(l)*(L/2)) + (pi/2)) )*soma;
end

figure(3)
plot(omegadis/pi,abs(H),'-b');grid minor
title('Resposta de Frequência');
xlabel('Omega normalizado'); ylabel('Magnitude');

%Quantização

%Criação de uma variavel chamada htemp para se manipular os coeficientes do
%filtro
htemp = h;

%Parte em que é determinado o bit de sinal e armazenado em hsinal
for m=1:length(htemp)
    
 
    if htemp(m) < 0        
        
        hsinal(m)= 1;
        
        htemp(m) = -1*htemp(m);
    
    else
        
        hsinal(m) = 0;
    end  
   
     % trasnformação em binário
    for s = 2:8
            htemp(m) = htemp(m)*2;
            
            if htemp(m) < 1  
            
                hbin(m,s) = 0;
            
            end 
            if htemp(m) > 1
                
                hbin(m,s) = 1;
                
                htemp(m) = htemp(m) - fix(htemp(m));
            end
            
            
    end
    
end


hquan = zeros(1,length(h));

for a = 1:length(h)
    
    %Passando os números de binario para decimal
    for m = 2:8 
        hquan(a) = hquan(a)+ hbin(a,m)*(2^(-m+1));
    end
    
    if hsinal(a) == 1 %Analizando o sinal dos coeficientes
        hquan(a) = hquan(a)*(-1);
    end
end

figure(1)
subplot(2,1,2)
stem(0:L,hquan)
title('Coeficientes Quantizados');
xlabel('índice n'); ylabel('valor do coeficiente');


% Calculo e plot de H3 com coeficientes de 8 bits equação 8.12
for l = 1:length(omegadis)
    for n = 0:(L-1)/2
        valores(n+1) = 2*hquan(n+1)*sin(omegadis(l)*(n - (L/2)));
    end
    soma = sum(valores);
    Hquan(l) = exp(-i* ( (omegadis(l)*(L/2)) + (pi/2)) )*soma;
end

figure(3)
plot(omegadis/pi,abs(Hquan),'-b');grid minor
title('Resposta de Frequência Quantizado');
xlabel('Omega normalizado'); ylabel('Magnitude');

%Gráficos
figure
subplot(2,1,1)
plot(omegadis/pi,abs(H),'-b');grid minor
hold

plot(omegadis/pi,abs(Hquan),':r'); grid minor
title('Módulo de H (azul- Módulo de Função H e vermelho- H quantizado)');
xlabel('Omega normalizado'); ylabel('Magnitude');

subplot(2,1,2)
plot(omegadis/pi,angle(H),'-b');grid minor
hold

plot(omegadis/pi,angle(Hquan),':r'); grid minor
title('Fase de H');
xlabel('Omega normalizado'); ylabel('Fase');


