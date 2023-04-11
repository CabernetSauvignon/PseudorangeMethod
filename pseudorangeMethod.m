function receiverCoord = pseudorangeMethod(satellite, delayTime)

  %% --------------------- Начальные параметры ---------------------------------
  receiverCoord = [0; 0; 0; 0]; % Начальные координаты НАП и начальное расхождение времени (deltaTau) 
  c = 2.99792458e+08; % Скорость света в вакууме, [м/с]

  %%---------------------- Рассчёт псевдодальностей ----------------------------
  P = delayTime .* c;  % вектор псевдодальностей исходя из времени задержки приёма сигнала

  %%---------------------- Оcновной цикл ---------------------------------------
  accuracy = 1;  % Значение точности текущего решения
  while accuracy > 0.00000001
    
    %%------------------ Рассчёт оценок псевдодальностей -----------------------
    % вектор псевдодальностей по МНК
    circumflexP = sqrt((receiverCoord(1) .- satellite(:, 1)).^2  .+ 
                       (receiverCoord(2) .- satellite(:, 2)).^2  .+ 
                       (receiverCoord(3) .- satellite(:, 3)).^2) .+ satellite(:, 4) .* c;
    deltaP = [P .- circumflexP];  % ошибки оценок псевдодальностей
  
    %%-------------- Рассчёт матрицы направляющих косинусов---------------------
    % расстояния до навигационного спутника
    R = sqrt((receiverCoord(1) .- satellite(:, 1)).^2  .+ 
             (receiverCoord(2) .- satellite(:, 2)).^2  .+ 
             (receiverCoord(3) .- satellite(:, 3)).^2);
    H = [];   % матрица направляющих косинусов
    H(:, 1) = (receiverCoord(1) .- satellite(:, 1)) ./ R(:);
    H(:, 2) = (receiverCoord(2) .- satellite(:, 2)) ./ R(:);
    H(:, 3) = (receiverCoord(3) .- satellite(:, 3)) ./ R(:);
    H(:, 4) = 1;
    
    %%------------------- Коррекция текущего решения----------------------------
    prevReceiverCoord = receiverCoord;  %Сохраняем предыдущее решение
    % Оценка ошибки текущего решения
    deltaNavigationReceiver = inv(transpose(H) * H) * transpose(H) * deltaP;
    receiverCoord = receiverCoord + deltaNavigationReceiver;
  
  %% Проверка достижения точности
  accuracy = sqrt((receiverCoord(1) - prevReceiverCoord(1)) ^ 2 + 
         (receiverCoord(2) - prevReceiverCoord(2)) ^ 2 + 
         (receiverCoord(3) - prevReceiverCoord(3)) ^ 2);
  endwhile

  receiverCoord = receiverCoord(1:3);
end

