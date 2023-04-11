clear all;

%% ------------------------------------------------------------------------
c = 2.99792458e+08; % Скорость света в вакууме, [м/с]
% ---------------------- Координаты НАП для проверки ----------------------
X_nav_receiver = 2.757469638336251e+06; % координата навигационного приемника (НАП), [м]
Y_nav_receiver = 1.616185590852821e+06; % координата навигационного приемника (НАП), [м]
Z_nav_receiver = 5.501093579768415e+06; % координата навигационного приемника (НАП), [м]
%% ------------------------------------------------------------------------
X1 = 2.195937759194386e+07; % координата 6 НС, [м]
Y1 = -1.677765131669631e+06; % координата 6 НС, [м]
Z1 = 1.289116411925952e+07; % координата 6 НС, [м]

X2 = 7.792103569304634e+06; % координата 7 НС, [м]
Y2 = -9.096796457026156e+06; % координата 7 НС, [м]
Z2 = 2.247747106121401e+07; % координата 7 НС, [м]

X3 = -1.012699752380525e+07; % координата 8 НС, [м]
Y3 = -1.215900544590689e+07; % координата 8 НС, [м]
Z3 = 1.993296089085508e+07; % координата 8 НС, [м]

X4 = 1.901730264632814e+07; % координата 9 НС, [м]
Y4 = -1.262741554696167e+07; % координата 9 НС, [м]
Z4 = 1.129322395574378e+07; % координата 9 НС, [м]

X5 = -1.911654053901512e+06; % координата 15 НС, [м]
Y5 = 1.593531304313390e+07; % координата 15 НС, [м]
Z5 = 1.984113231665855e+07; % координата 15 НС, [м]

X6 = -2.230243011054459e+06; % координата 17 НС, [м]
Y6 = 1.922436416999018e+07; % координата 17 НС, [м]
Z6 = 1.659128674288339e+07; % координата 17 НС, [м]

X7 = 6.406348881553141e+06; % координата 18 НС, [м]
Y7 = 2.462733100637249e+07; % координата 18 НС, [м]
Z7 = 1.251121457619099e+06; % координата 18 НС, [м]

X8 = -1.041100304385226e+07; % координата 24 НС, [м]
Y8 = 2.867677014174020e+06; % координата 24 НС, [м]
Z8 = 2.308987787737002e+07; % координата 24 НС, [м]
satelliteXYZ = [X1, Y1, Z1; 
                X2, Y2, Z2;
                X3, Y3, Z3;
                X4, Y4, Z4;
                X5, Y5, Z5;
                X6, Y6, Z6;
                X7, Y7, Z7;
                X8, Y8, Z8];
%% ------------------------------------------------------------------------
tau1 = 0.069504426668593; % время задержки, [с] в приеме радионавигационного сигнала от 6 НС
tau2 = 0.069033498355289; % время задержки, [с] в приеме радионавигационного сигнала от 7 НС
tau3 = 0.079220185830764; % время задержки, [с] в приеме радионавигационного сигнала от 8 НС
tau4 = 0.074647673609317; % время задержки, [с] в приеме радионавигационного сигнала от 9 НС
tau5 = 0.069368087242531; % время задержки, [с] в приеме радионавигационного сигнала от 15 НС
tau6 = 0.071379412978226; % время задержки, [с] в приеме радионавигационного сигнала от 17 НС
tau7 = 0.078998330236525; % время задержки, [с] в приеме радионавигационного сигнала от 18 НС
tau8 = 0.073409889758624; % время задержки, [с] в приеме радионавигационного сигнала от 24 НС
tau = [tau1; tau2; tau3; tau4; tau5; tau6; tau7; tau8];
%% ------------------------------------------------------------------------
deltaTau1 = 0.0; % расхождение шкал времени НАП с 6 НС
deltaTau2 = 0.0; % расхождение шкал времени НАП с 7 НС
deltaTau3 = 0.0; % расхождение шкал времени НАП с 8 НС
deltaTau4 = 0.0; % расхождение шкал времени НАП с 9 НС
deltaTau5 = 0.0; % расхождение шкал времени НАП с 15 НС
deltaTau6 = 0.0; % расхождение шкал времени НАП с 17 НС
deltaTau7 = 0.0; % расхождение шкал времени НАП с 18 НС
deltaTau8 = 0.0; % расхождение шкал времени НАП с 24 НС
deltaTau = [deltaTau1; deltaTau2; deltaTau3; deltaTau4; deltaTau5; deltaTau6; deltaTau7; deltaTau8];
%% ------------------------------------------------------------------------
satelliteXYZdeltaTau = [satelliteXYZ deltaTau];
%% ------------------------------------------------------------------------
navigationReceiver = [0; 0; 0; 0]; % Начальные координаты НАП и начальное расхождение времени (deltaTau) 

tic
%%------Рассчёт псевдодальностей----
P = tau .* c; % вектор псевдодальностей

%% Основной цикл
accuracy = 1;
while accuracy > 0.00000001
  %% Рассчёт оценок псевдодальностей  
    LeastSquareP = sqrt((navigationReceiver(1) .- satelliteXYZdeltaTau(:, 1)).^2  .+ 
                     (navigationReceiver(2) .- satelliteXYZdeltaTau(:, 2)).^2  .+ 
                     (navigationReceiver(3) .- satelliteXYZdeltaTau(:, 3)).^2) .+ satelliteXYZdeltaTau(:, 4) .* c;

  deltaP = [P .- LeastSquareP];
  
  %% Рассчёт матрицы направляющих косинусов
    R = sqrt((navigationReceiver(1) .- satelliteXYZdeltaTau(:, 1)).^2  .+ 
                (navigationReceiver(2) .- satelliteXYZdeltaTau(:, 2)).^2  .+ 
                (navigationReceiver(3) .- satelliteXYZdeltaTau(:, 3)).^2);
    H = [];
    H(:, 1) = (navigationReceiver(1) .- satelliteXYZdeltaTau(:, 1)) ./ R(:);
    H(:, 2) = (navigationReceiver(2) .- satelliteXYZdeltaTau(:, 2)) ./ R(:);
    H(:, 3) = (navigationReceiver(3) .- satelliteXYZdeltaTau(:, 3)) ./ R(:);
    H(:, 4) = 1;
  
  %% Коррекция текущего решения
  deltaNavigationReceiver = inv(transpose(H) * H) * transpose(H) * deltaP;
  previousNavigationReceiver = navigationReceiver;
  navigationReceiver = navigationReceiver + deltaNavigationReceiver;
  
  %% Проверка достижения точности
  accuracy = sqrt((navigationReceiver(1) - previousNavigationReceiver(1)) ^ 2 + 
                  (navigationReceiver(2) - previousNavigationReceiver(2)) ^ 2 + 
                  (navigationReceiver(3) - previousNavigationReceiver(3)) ^ 2);
endwhile

clc;
toc
disp("navigationReceiver"), disp(navigationReceiver)
disp("accuracy"), disp(accuracy)
disp("Рассчитанные координаты НАП"), disp(navigationReceiver(1:3))
disp("Координаты НАП для проверки"), disp(X_nav_receiver), disp(Y_nav_receiver), disp(Z_nav_receiver)

function receiverCoord = pseudorangeMethod(satellite, delayTime)


  navigationReceiver = [0; 0; 0; 0]; % Начальные координаты НАП и начальное расхождение времени (deltaTau) 

  %%------Рассчёт псевдодальностей----
  P = []; % вектор псевдодальностей
  for i = 1:8
    P(1, i) = c * tau(i);
  endfor

  %% Основной цикл
  accuracy = 1;
  steps = 0;
  while accuracy > 0.00000001
    %% Рассчёт оценок псевдодальностей
    [numRowsSatellite,numColsSatellite] = size(satelliteXYZdeltaTau);
    LeastSquareP = [];
    deltaP = [];
    for i = 1:numRowsSatellite
      LeastSquareP(i) = sqrt((navigationReceiver(1) - satelliteXYZdeltaTau(i, 1)).^2  + 
                      (navigationReceiver(2) - satelliteXYZdeltaTau(i, 2)).^2  + 
                      (navigationReceiver(3) - satelliteXYZdeltaTau(i, 3)).^2) + satelliteXYZdeltaTau(i, 4) .* c;
    endfor

    deltaP = [P .- LeastSquareP];
  
    %% Рассчёт матрицы направляющих косинусов
    H = [];
    R = [];
    for i = 1:numRowsSatellite
      R(i) = sqrt((navigationReceiver(1) - satelliteXYZdeltaTau(i, 1)).^2  + 
                 (navigationReceiver(2) - satelliteXYZdeltaTau(i, 2)).^2  + 
                 (navigationReceiver(3) - satelliteXYZdeltaTau(i, 3)).^2);
                
      H(i, 1) = (navigationReceiver(1) - satelliteXYZdeltaTau(i, 1)) / R(i);
      H(i, 2) = (navigationReceiver(2) - satelliteXYZdeltaTau(i, 2)) / R(i);
      H(i, 3) = (navigationReceiver(3) - satelliteXYZdeltaTau(i, 3)) / R(i);
      H(i, 4) = 1;
    endfor
  
    %% Коррекция текущего решения
    %navigationReceiver = navigationReceiver + transpose(H) * transpose(deltaP);
    deltaNavigationReceiver = inv(transpose(H) * H) * transpose(H) * transpose(deltaP);
    previousNavigationReceiver = navigationReceiver;
    navigationReceiver = navigationReceiver + deltaNavigationReceiver;
  
    %% Проверка достижения точности
    accuracy = sqrt((navigationReceiver(1) - previousNavigationReceiver(1)) .^ 2 + 
                    (navigationReceiver(2) - previousNavigationReceiver(2)) .^ 2 + 
                    (navigationReceiver(3) - previousNavigationReceiver(3)) .^ 2);
  
    steps++
  endwhile

end

