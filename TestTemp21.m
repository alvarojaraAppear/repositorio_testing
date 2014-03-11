% - Implementation of algorithms for QRS detection from ECG signals using
% TMS320C6713 processor platform; Geoffrey Green
%

close all; clear all; clc

%% DATA 

procesada = xlsread('DATA.xlsx',1,'d7000:d15500'); % Sheet 1: 63684
% procesada = xlsread('DATA.xlsx',2,'d1:d10000')'; % Sheet 2 62443
% procesada = xlsread('DATA.xlsx',3,'d1000:d32050'); % Sheet 3: 58180
% procesada = xlsread('DATA.xlsx',4,'d9000:d32050'); % Sheet 4: 54851

data = max(procesada) - procesada;

for i = 1 : length(data)        % se revisa que no falten datos en el vector
    if isnan(data(i)) == 1
        data(i) = data(i - 1);  % si faltan se completa con el dato anterior 
    end
    
end

fs = 83;
T = length(data)/fs;
t = 0 : 1/(length(data)/T) : (T - (1/(length(data)/T))) ;

plot(t, data), grid on, title('Data original')
figure

%% Pan Tompkins

newData = data;                                       % Solo para seguir la notacion utilizada en las versiones de prueba

% Inicialización de variables

lengthLPArray = 12;                                   % Arreglo que almacena los 12 estados previos
LPArray = zeros(1, lengthLPArray);                    % LPArray = 0; Inicialmente todos los estados previos = 0
LPOutput = zeros(1, length(newData));                 % LPOutput = 0;
LPOutputPrev = [0 0];                                 % LPOutputPrev(1) = estado previo (n - 1) ; LPOutputPrev(2) = estado previo previo (n - 2)

lengthHPArray = 32;                                   % Arreglo que almacena los 32 estados previos
HPArray = zeros(1, lengthHPArray);                    % HPArray = 0; Inicialmente todos los estados previos = 0
HPOutput = zeros(1, length(newData));                 % HPOutput = 0;
HPOutputPrev = 0;                                     % HPOutputPrev = salida previa (n - 1), inicialmente = 0    

lengthDerivateArray = 4;                              % Arreglo que almacena los 4 valores utilizados para derivar un dato
DerivateArray = zeros(1, lengthDerivateArray);        % DerivateArray = 0;
DerivateOutput = zeros(1, length(newData));           % LPOutput = 0;

SquareOutput = zeros(1, length(newData));             % Arreglo que almacena el resultado de elevar cada dato al cuadrado 

lengthIntegrationArray = 28;                          % Ventana de integracion de 30 datos                    
IntegrationArray = zeros(1, lengthIntegrationArray);  % IntegrationArray = 0;
IntegrationOutput = zeros(1, length(newData));        % IntegrationOutput = 0;

lengthArrayWin = 83*20;                                 % fs = 83; 5sg = 83*5 = 415
ArrayWin = zeros(1, lengthArrayWin);
normArrayWin = zeros(1, lengthArrayWin);
maximoArrayWin = 1;

lengthDataRealWin = lengthArrayWin;
dataRealWin = zeros(1, lengthDataRealWin);
tempDataRealWin = zeros(1, lengthDataRealWin);
maximoDataRealWin = 1;

datoAnteriorFiltro = 0;
threshold = 0;
largoMaximosLocales = lengthArrayWin;
maximosLocales = zeros(1, largoMaximosLocales);
ondaR = zeros(1, largoMaximosLocales);

bpm = 0;
bpmPrevio = 0;
k = 60/(lengthArrayWin/fs)

arrayBPM = zeros(1,length(data));

% Se aplica el algoritmo a UN dato y se almacena el resultado en la
% variabla resultado. Se hace ésto con todos lo s datos, uno a la vez

valorFinal = 0;

for n = 1 : length(newData)                                                      
    
    tic
    
    %% --> LowPass Filter: y(n) = 2*y(n - 1) - y(n - 2) + (1/32)*(x(n) - 2*x(n - 6) + x(n - 12))
    
     LPOutput(n) = 2*LPOutputPrev(1) - LPOutputPrev(2) + (1/32)*(newData(n) - 2*LPArray(6) + LPArray(12)); % Se filtra el dato n
    
     for index = (lengthLPArray - 1) : -1 : 1         % Se actualiza el arreglo para el sgte. dato 
         temp = LPArray(index);
         LPArray(index + 1) = temp;
     end
     LPArray(1) = newData(n);
     
     LPOutputPrev(2) = LPOutputPrev(1);               % Se actualizan los estados anteriores para el sgte. dato
     LPOutputPrev(1) = LPOutput(n);
    
     %% --> HighPass Filter: y(n) = y(n - 1) - (1/32)*x(n) + x(n - 16) - x(n - 17) + (1/32)*x(n - 32)
    
     HPOutput(n) = HPOutputPrev - (1/32)*LPOutput(n) + HPArray(16) - HPArray(17) + (1/32)*HPArray(32); % Se filtra el dato n

     for index = (lengthHPArray - 1) : -1 : 1         % Se actualiza el arreglo para el sgte. dato
         temp = HPArray(index);
         HPArray(index + 1) = temp;
     end
     HPArray(1) = LPOutput(n);
     
     HPOutputPrev = HPOutput(n);                      % Se actualiza el estado anterior para el sgte. dato

     %% --> Derivate: y(n) = (1/8)*( 2*x(n) + x(n - 1) - x(n - 3) - 2*x(n - 4) )

     DerivateOutput(n) = (1/8)*( 2*HPOutput(n) + DerivateArray(2) - DerivateArray(3) - 2*DerivateArray(4) );  % Se deriva el dato n     
     
     for index = (lengthDerivateArray - 1) : -1 : 1    % Se actualiza el arreglo para el sgte. dato
         temp = DerivateArray(index);
         DerivateArray(index + 1) = temp;
     end
     DerivateArray(1) = HPOutput(n);

     %% --> Squaring Point to point
     
     SquareOutput(n) = (DerivateOutput(n)^2);          % Se eleva al cuadrado el dato n

     %% --> Integration: window size = 30;
     
     IntegrationOutput(n) = SquareOutput(n);           % Se integra el dato n 
     
     for index = 2 : lengthIntegrationArray
         IntegrationOutput(n) = IntegrationOutput(n) + IntegrationArray(index);
     end
     for index = (lengthIntegrationArray - 1) : -1 : 1 % Se actualiza la ventana para el sgte. dato
         temp = IntegrationArray(index);
         IntegrationArray(index + 1) = temp;
     end
     IntegrationArray(1) = SquareOutput(n);
     
     %% Beat Detection
     
     % Se ordena el arreglo que contiene los datos de la ventana a analizar
     for index = 2 : 1 : lengthArrayWin                                           
         temp = ArrayWin(index);     
         ArrayWin(index - 1) = temp;
  
         temp = dataRealWin(index); 
         dataRealWin(index - 1) = temp;          
     end  
     ArrayWin(lengthArrayWin) = IntegrationOutput(n);     
     dataRealWin(lengthArrayWin) = newData(n);
     
     % Falta Suavizar la señal!!!
     
     % Busqueda del maximo de la ventana 
     for index = 1 : 1 : lengthArrayWin             
           
         if ArrayWin(index) >= maximoArrayWin;   
             maximoArrayWin = ArrayWin(index);         
         end         
         if dataRealWin(index) >= maximoDataRealWin;   
             maximoDataRealWin = dataRealWin(index); 
         end        
     end
     
     % Ajuste de amplitud  [0 - 1]
        
     for index = 1 : lengthArrayWin    
         normArrayWin(index) = ArrayWin(index)/maximoArrayWin; 
         tempDataRealWin(index) = dataRealWin(index)/maximoDataRealWin; 
     end   
     
     % Umbral inicial
     % threshold = mean(tempDataRealWin);
     threshold = 0.7;
     
     % Busqueda  de maximos locales  
     ventana = 18;                                           % Promedio entre caso normal (25 ptos. por lado) y peor caso (11 puntos por lado)
     
     for index = 1 + ventana : (lengthArrayWin - ventana)
           
         % Ventana de comparacion
         flag = 1;
         
         for ind = index - ventana : 1  : index + ventana
             if ( normArrayWin(index) < normArrayWin(ind) ) 
                 flag = 0;
             end
         end
         
         if flag == 1;
             ondaR(index) = normArrayWin(index);
             primeraOndaR = index;
             bpm = bpm + 1;
         end
             
     end
     
     %pause
     pause(0.0120);
     clf  
     bpm
     bpmPrevio
     
     bpm = (bpm*0.1 + bpmPrevio*0.9);
     bpmPrevio = bpm;
     bpm = bpm*k
     arrayBPM(n) = bpm;
     
     if bpm >= 220 
         titulo = 'Original - PanTompkins; bpm > 220 ';
     else
        bpm = uint64(bpm);
        titulo = strcat('Original - PanTompkins; bpm = ', num2str(bpm)); 
     end
     plot(tempDataRealWin, 'b--'), title(titulo), grid on
     hold on, plot(normArrayWin, 'r')
     hold on, plot(ondaR, 'ok'), grid on
     
     maximoArrayWin = ArrayWin(1);
     datoAnteriorFiltro = 0;
     maximosLocales = zeros(1, largoMaximosLocales);
     ondaR = zeros(1, largoMaximosLocales);
     
     bpm = 0;
     
     toc
     
    % Fin Beat Detection
     
end

disp('final')
arrayBPM = arrayBPM(lengthArrayWin:length(arrayBPM));
figure, plot(arrayBPM), grid on, title(strcat('Promedio: ', num2str(mean(arrayBPM)), '; std: ', num2str(std(arrayBPM))))

% resultado = conv(IntegrationOutput, gausswin(17));     % Se suaviza el resultado para buscar los máximos (Pan - Tompkins v1.1)
% resultado = IntegrationOutput;
% 
% figure
% plot(data/max(data), 'b--'), title('FINAL'), grid on
% hold on, plot(resultado/max(resultado), 'r')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     