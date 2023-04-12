function [receiverCoord, iterationsNumber] = pseudorangeMethod(satellite, delayTime)
%%--------------------Проверка аргументов функции-----------------------------
% usage = strjoin({"Input satellite matrix must have 4 columns (x, y, z, deltaTau),",
%              "that means func need a number of satellite coords and", 
%              "divergence of time scales for each satellite.", 
%              "\nAlso this func requires at least 4 rows of arguments", 
%              "that means you need to input data about at least 4 satellites.",
%              "\n delayTime argument is signal acceptance delay time for each", 
%              "satellite. It must be a vector of arguments, which equals to", 
%              "number of satellite rows."});
  
[satelliteRows, satelliteCols] = size(satellite);
if satelliteCols ~= 4
    error("Wrong nuber of columns in satellite argument. Input satellite matrix must have 4 columns (x, y, z, deltaTau).");
end % if

if satelliteRows < 4
    error("Wrong nuber of rows in satellite argument. Number must be equal or greater than 4");
end % if

[delayTimeRows, delayTimeCols] = size(delayTime);
if delayTimeRows ~= satelliteRows ||  delayTimeCols ~= 1
    error("Wrong nuber of rows in delayTime argument. Number of rows must be equal to satellite number of rows");
end % if
  
%% --------------------- Начальные параметры ---------------------------------
receiverCoord = [0; 0; 0; 0]; % Начальные координаты НАП и начальное расхождение времени (deltaTau) 
c = 2.99792458e+08; % Скорость света в вакууме, [м/с]

%% ---------------------- Рассчёт псевдодальностей ----------------------------
P = delayTime .* c;  % вектор псевдодальностей исходя из времени задержки приёма сигнала
%% ---------------------- Оcновной цикл ---------------------------------------
accuracy = 1;  % Значение точности текущего решения
accuracyRequired = 0.00000001;  % Требуемая точность
iterationsNumber = 0 % Счётчик итераций цикла
while accuracy > accuracyRequired
%% ------------------ Рассчёт оценок псевдодальностей -----------------------
% вектор псевдодальностей по МНК
    circumflexP = sqrt((receiverCoord(1) - satellite(:, 1)).^2  + ...
        (receiverCoord(2) - satellite(:, 2)).^2  + ...
        (receiverCoord(3) - satellite(:, 3)).^2) + satellite(:, 4) .* c;
    deltaP = P - circumflexP;  % ошибки оценок псевдодальностей
%% ---------------- Рассчёт матрицы направляющих косинусов ----------------
% расстояния до навигационного спутника
    R = sqrt((receiverCoord(1) - satellite(:, 1)).^2  + ...
        (receiverCoord(2) - satellite(:, 2)).^2  + ...
        (receiverCoord(3) - satellite(:, 3)).^2);
    H = [];   % матрица направляющих косинусов
    
    H(:, 1) = (receiverCoord(1) - satellite(:, 1)) ./ R(:);
    H(:, 2) = (receiverCoord(2) - satellite(:, 2)) ./ R(:);
    H(:, 3) = (receiverCoord(3) - satellite(:, 3)) ./ R(:);
    H(:, 4) = 1;
%% ---------------------- Коррекция текущего решения ----------------------
    prevReceiverCoord = receiverCoord;  %Сохраняем предыдущее решение
%% -------------------- Оценка ошибки текущего решения --------------------
    deltaNavigationReceiver = inv(transpose(H) * H) * transpose(H) * deltaP;
    receiverCoord = receiverCoord + deltaNavigationReceiver;
  
%% Проверка достижения точности
  accuracy = sqrt((receiverCoord(1) - prevReceiverCoord(1)) ^ 2 + ...
      (receiverCoord(2) - prevReceiverCoord(2)) ^ 2 + ...
      (receiverCoord(3) - prevReceiverCoord(3)) ^ 2);
%% Увеличиваем счётчик итераций на 1       
  iterationsNumber = iterationsNumber + 1;
end % while
receiverCoord = receiverCoord(1:3);
end % function

