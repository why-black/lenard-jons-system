function plot_thermodynamics()
    try
        files = {
            'system_temperature.txt', 'Temperature'
            'system_pressure.txt', 'Pressure'
            'system_kinetic_energy.txt', 'KineticEnergy'
            'system_potential_energy.txt', 'PotentialEnergy'
            'system_total_energy.txt', 'TotalEnergy'
        };
        
        % Очищаем фигуру перед построением
        figure('Position', [100, 100, 1200, 800]);
        hold on;
        
        % Для каждого файла
        for i = 1:length(files)
            filename = files{i,1};
            if ~exist(filename, 'file')
                error('Файл %s не найден', filename);
            end
            
            % Чтение данных с обработкой ошибок
            try
                data = load(filename);
                
                % Проверка размера данных
                if size(data, 2) < 2
                    error('Файл %s содержит менее 2 столбцов', filename);
                end
                
                % Сортировка по времени (первый столбец)
                data = sortrows(data, 1);
                
                time = data(:,1);
                values = data(:,2);
                
                % Создание субграфиков
                subplot(3, 2, i);
                plot(time, values, 'LineWidth', 1.5);
                
                % Настройка подписей
                switch files{i,2}
                    case 'Temperature'
                        title('Температура системы');
                        ylabel('T');
                    case 'Pressure'
                        title('Давление системы');
                        ylabel('P');
                    case 'KineticEnergy'
                        title('Кинетическая энергия');
                        ylabel('E_k');
                    case 'PotentialEnergy'
                        title('Потенциальная энергия');
                        ylabel('E_p');
                    case 'TotalEnergy'
                        title('Полная энергия');
                        ylabel('E_{total}');
                end
                
                xlabel('Время');
                grid on;
                
            catch ME
                warning('Ошибка обработки %s: %s', filename, ME.message);
            end
        end
        
        % Дополнительный график - фазовый портрет
        subplot(3, 2, 6);
        
        % Объединение данных (пример для T-P)
        try
            temp_data = load('system_temperature.txt');
            pres_data = load('system_pressure.txt');
            
            % Сортировка и обрезка по минимальной длине
            min_len = min(size(temp_data, 1), size(pres_data, 1));
            time = temp_data(1:min_len, 1);
            temperature = temp_data(1:min_len, 2);
            pressure = pres_data(1:min_len, 2);
            
            % Разделение на фазы
            phase = ones(size(time));
            phase(time > max(time)/2) = 2;
            
            scatter(temperature, pressure, 10, phase, 'filled');
            title('Фазовый портрет: T-P');
            xlabel('Температура'); ylabel('Давление');
            colormap(jet(2));
            colorbar('Ticks', [1,2], 'TickLabels', {'Этап I', 'Этап II'});
            grid on;
            
        catch
            warning('Не удалось построить фазовый портрет');
        end
        
        saveas(gcf, 'thermodynamics_plot.png');
        disp('График сохранен как thermodynamics_plot.png');
        
    catch ME
        disp('Ошибка:');
        disp(ME.message);
    end
end