function plot_rdf()
    % Данные для жидкостей
    phases = {
        'liquid', 'b', '-', 'ρ=0.75σ^{-3}'
        'compressed', 'r', '--', 'ρ=1.05σ^{-3}'
    };
    
    figure;
    hold on;
    grid on;
    
    for i = 1:size(phases, 1)
        filename = ['rdf_' phases{i,1} '.csv'];
        % Явно указываем имена переменных при чтении
        data = readtable(filename, 'ReadVariableNames', true);
        
        % Используем индексы столбцов вместо имен
        plot(data{:,1}, data{:,2}, ...  % Столбец 1 = r, столбец 2 = g_r
            'Color', phases{i,2}, ...
            'LineStyle', phases{i,3}, ...
            'DisplayName', phases{i,4}, ...
            'LineWidth', 2);
    end
    
    title('Функция радиального распределения');
    xlabel('r (σ)');
    ylabel('g(r)');
    legend show;
    legend('Location', 'northeast'); % Явно указываем положение легенды
    
    saveas(gcf, 'rdf_comparison.png');
end