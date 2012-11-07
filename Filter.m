%Digital Filter FIR with linear phase - EQUATION T4 - BANDPASS

clear all
close all

%Order of the filter (Odd, value given)
L = 31; 

%Low-frequency cutoff (Value Given)
fcl = 7e3;

%Upper cut-off frequency (Given Value)
fc2 = 28e3;

%Sampling frequency (Given Value)
Fs = 70e3;

%Cutoff frequency standard because it requires by the function firls
freq2 = 2*fc2/Fs;
freq1 = 2*fcl/Fs;

%Vector with cutoff frequencies
freq = [0 freq1/2 freq1 freq2 (1-freq1/2) 1];

%Vector Amplitudes to interpolate, required by the matlab function
amp = [0 0 1 1 0 0];


%Calculation of coefficients h
h = firls(L,freq,amp,'h');

%Graph discrete coefficients h
figure(1)
subplot(2,1,1)
stem(0:L,h)
title('Ratios (odd symmetry)');
xlabel('index n'); ylabel('coefficient');

%Interval omega discrete points that will be used in the graph
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
title('Frequency Response');
xlabel('Omega standard'); ylabel('Magnitude');

%Quantization

%Creating a variable called htemp to manipulate the filter coefficients

htemp = h;

%Part where is determined that the sign bit and stored in hsinal
for m=1:length(htemp)
    
 
    if htemp(m) < 0        
        
        hsinal(m)= 1;
        
        htemp(m) = -1*htemp(m);
    
    else
        
        hsinal(m) = 0;
    end  
   
     % Converting in binary
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
title('Quantized coefficients');
xlabel('index n'); ylabel('cofficient');


% Calculate and plot H3 8-bit coefficients with equation 8.12
for l = 1:length(omegadis)
    for n = 0:(L-1)/2
        valores(n+1) = 2*hquan(n+1)*sin(omegadis(l)*(n - (L/2)));
    end
    soma = sum(valores);
    Hquan(l) = exp(-i* ( (omegadis(l)*(L/2)) + (pi/2)) )*soma;
end

figure(3)
plot(omegadis/pi,abs(Hquan),'-b');grid minor
title('Frequency Response quantized');
xlabel('Omega standard'); ylabel('Magnitude');

%Grafics
figure
subplot(2,1,1)
plot(omegadis/pi,abs(H),'-b');grid minor
hold

plot(omegadis/pi,abs(Hquan),':r'); grid minor
title('Module H (blue-Function Module H and H red-quantized)');
xlabel('Omega standard'); ylabel('Magnitude');

subplot(2,1,2)
plot(omegadis/pi,angle(H),'-b');grid minor
hold

plot(omegadis/pi,angle(Hquan),':r'); grid minor
title('Phase of H');
xlabel('Omega standard'); ylabel('Phase');


