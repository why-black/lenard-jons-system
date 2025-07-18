using System;
using System.IO;
using System.Globalization;

class MD_Simulation
{
    // Параметры системы
    int N = 1372;
    int M = 20000;
    int nx = 7, ny = 7, nz = 7;
    double alat = 1.75;
    double dt = 0.005;
    double T_ext = 2.8;
    double T0 = 0.01;
    double vel_dT = 0.2;
    double dr_s = 0.3;
    double r_cut = 2.5;
    double L_x, L_y, L_z;
    double time = 0;
    // Для итерационного алгоритма
    double[] v_x = Array.Empty<double>();
    double[] v_y = Array.Empty<double>();
    double[] v_z = Array.Empty<double>();
    double[] v_xo = Array.Empty<double>();
    double[] v_yo = Array.Empty<double>();
    double[] v_zo = Array.Empty<double>();
    double[] b_x = Array.Empty<double>();
    double[] b_y = Array.Empty<double>();
    double[] b_z = Array.Empty<double>();
    double p_o, p_s;
    double ds, dd, err;
    int iter;
    bool ready;
    // Параметры для NPT-ансамбля
    double Q_T = 100.0;   // Параметр термостата
    double P_ext = 25.0; // Внешнее давление
    double Q_B = 50.0;   // Параметр баростата
    double eta;           // Параметр масштабирования объема
    double eta_dot;       // Производная параметра eta
    bool isNPT = true;    // Флаг для активации NPT-режима
    double xi;      // Параметр баростата
    double S_f;     // Интеграл термостата


    // Физические величины
    double lambda = 0; // параметр термостата
    double virial = 0;
    double kinetic_energy = 0;
    double potential_energy = 0;
    double H_energy = 0;
    double temperature = 0;
    double pressure = 0;
    double volume = 0;
    double density = 0;

    // Массивы частиц (инициализированы пустыми массивами)
    double[] pos_x = Array.Empty<double>();
    double[] pos_y = Array.Empty<double>();
    double[] pos_z = Array.Empty<double>();
    double[] vel_x = Array.Empty<double>();
    double[] vel_y = Array.Empty<double>();
    double[] vel_z = Array.Empty<double>();
    double[] acc_x = Array.Empty<double>();
    double[] acc_y = Array.Empty<double>();
    double[] acc_z = Array.Empty<double>();
    double[] e_pot = Array.Empty<double>();

    // Список Верле
    int[] List = Array.Empty<int>();
    int[] Marker1 = Array.Empty<int>();
    int[] Marker2 = Array.Empty<int>();
    int[] Advance = Array.Empty<int>();
    double[] disp_list_x = Array.Empty<double>();
    double[] disp_list_y = Array.Empty<double>();
    double[] disp_list_z = Array.Empty<double>();
    bool update_test = true;
    int list_size = 0, max_list_size = 0;

    // Таблицы потенциала
    double[,] table_pot = new double[0, 0];
    int table_size = 0;
    double sigma = 1.0, epsilon = 1.0;

    // Вспомогательные переменные
    Random rand = new Random();

    static void Main(string[] args)
    {
        MD_Simulation md = new MD_Simulation();

        // ЭТАП I: ПЛАВЛЕНИЕ (NVT-ансамбль) 
        Console.WriteLine("ЭТАП I: ПЛАВЛЕНИЕ КРИСТАЛЛА (NVT)");

        // Параметры моделирования для этапа I
        md.N = 1372;          // Количество частиц
        md.M = 20000;         // Всего шагов
        md.dt = 0.005;        // Временной шаг
        md.T0 = 0.01;         // Начальная температура
        md.T_ext = 2.8;       // Целевая температура
        md.vel_dT = 0.2;      // Скорость нагрева
        md.isNPT = false;     // Режим NVT
        md.nx = 7; md.ny = 7; md.nz = 7;  // Размеры решетки
        md.alat = 1.75;       // Параметр решетки

        // Инициализация системы
        md.InitializeSystem();
        md.GenerateFCCLattice();
        md.InitializeVelocities();

        // Создаем файлы для записи (очищаем существующие)
        File.WriteAllText("system_config.lammpstrj", string.Empty);
        File.WriteAllText("thermodynamics.csv", "time,temperature,pressure,kinetic_energy,potential_energy,total_energy\n");

        // Главный цикл моделирования для этапа I

        for (int step = 1; step <= md.M; step++)
        {
            md.time = step * md.dt;

            // 1. Половина шага интегрирования (обновление координат)
            md.SystemDynamics(1);

            // 2. Применение периодических граничных условий
            md.ApplyBoundaryConditions();

            // 3. Расчет сил
            md.CalculateForces();

            // 4. Вторая половина шага интегрирования (обновление скоростей)
            md.SystemDynamics(2);

            // Сохранение данных каждые 100 шагов
            if (step % 100 == 0)
            {
                md.SaveLAMMPSTrajectory(step);
                md.SaveThermodynamicParameters();

                Console.WriteLine($"Этап I: Шаг {step}/{md.M}, " +
                                 $"Время: {md.time:F2}τ, " +
                                 $"T: {md.temperature:F4}, " +
                                 $"P: {md.pressure:F4}");
            }

            // Линейное увеличение температуры до шага 2800
            if (step <= 2800)
            {
                md.T0 = 0.01 + (2.8 - 0.01) * step / 2800.0; // Линейный нагрев
            }
            else
            {
                md.T0 = 2.8; // Фиксированная температура
            }
        }

        // Сохранение финального состояния системы
        md.SaveSystemState("system_pos_vel_acc.txt");
        Console.WriteLine("Этап I завершен. Состояние системы сохранено.");
        md.CalculateRDF("liquid", begin_step: 100, end_step: 20000); // ρ≈0.75

        // ЭТАП II: СЖАТИЕ (NPT-ансамбль)
        Console.WriteLine("\nЭТАП II: СЖАТИЕ ЖИДКОСТИ (NPT)");

        // Параметры моделирования для этапа II
        md.M = 20000;         // Всего шагов
        md.isNPT = true;      // Режим NPT
        md.T0 = 2.8;          // Начальная температура
        md.T_ext = 2.8;       // Целевая температура (постоянная)
        md.vel_dT = 0.0;      // Без изменения температуры
        md.P_ext = 25.0;      // Внешнее давление
        md.Q_B = 50.0;        // Параметр баростата

        // Загрузка состояния системы из этапа I
        md.LoadSystemFromFile("system_pos_vel_acc.txt");

        // Пересчет начальных параметров
        md.volume = md.L_x * md.L_y * md.L_z;
        md.density = md.N / md.volume;

        // Главный цикл моделирования для этапа II
        for (int step = 1; step <= md.M; step++)
        {
            md.time = (md.M + step) * md.dt;  // Продолжаем время

            // 1. Половина шага интегрирования (обновление координат)
            md.SystemDynamics(1);

            // 2. Применение периодических граничных условий
            md.ApplyBoundaryConditions();

            // 3. Расчет сил
            md.CalculateForces();

            // 4. Вторая половина шага интегрирования (обновление скоростей)
            md.SystemDynamics(2);

            // Сохранение данных каждые 100 шагов
            if (step % 100 == 0)
            {
                md.SaveLAMMPSTrajectory(md.M + step);
                md.SaveThermodynamicParameters();

                Console.WriteLine($"Этап II: Шаг {step}/{md.M}, " +
                                 $"Время: {md.time:F2}τ, " +
                                 $"T: {md.temperature:F4}, " +
                                 $"P: {md.pressure:F4}, " +
                                 $"ρ: {md.density:F4}, " +
                                 $"L: {md.L_x:F2}×{md.L_y:F2}×{md.L_z:F2}");
            }
        }

        // Финализация
        md.SaveSystemState("system_final.txt");
        Console.WriteLine("Этап II завершен. Начинаем расчет RDF...");

        // Расчет функции радиального распределения
        md.CalculateRDF("compressed", begin_step: 20001, end_step: 40000); // ρ≈1.05

        Console.WriteLine("\nМОДЕЛИРОВАНИЕ УСПЕШНО ЗАВЕРШЕНО!");
        Console.WriteLine("Результаты сохранены в файлы:");
        Console.WriteLine("- Траектория: system_config.lammpstrj");
        Console.WriteLine("- Термодинамика: thermodynamics.csv");
        Console.WriteLine("- Финальное состояние: system_final.txt");
        Console.WriteLine("- RDF: rdf.csv");
        Console.ReadLine();
    }

    public void LoadSystemFromFile(string filename)
    {
        try
        {
            using (StreamReader reader = new StreamReader(filename))
            {
                L_x = double.Parse(reader.ReadLine() ?? "0", CultureInfo.InvariantCulture);
                L_y = double.Parse(reader.ReadLine() ?? "0", CultureInfo.InvariantCulture);
                L_z = double.Parse(reader.ReadLine() ?? "0", CultureInfo.InvariantCulture);

                for (int i = 0; i < N; i++)
                {
                    string? line = reader.ReadLine();

                    // Проверяем на конец файла
                    if (line == null)
                    {
                        throw new Exception($"Файл содержит только {i} строк из требуемых {N}");
                    }

                    string[] parts = line.Split();

                    // Проверяем количество элементов в строке
                    if (parts.Length < 6)
                    {
                        throw new Exception($"Строка {i + 1} содержит недостаточно данных: {parts.Length} элементов вместо 6");
                    }

                    pos_x[i] = double.Parse(parts[0], CultureInfo.InvariantCulture);
                    pos_y[i] = double.Parse(parts[1], CultureInfo.InvariantCulture);
                    pos_z[i] = double.Parse(parts[2], CultureInfo.InvariantCulture);
                    vel_x[i] = double.Parse(parts[3], CultureInfo.InvariantCulture);
                    vel_y[i] = double.Parse(parts[4], CultureInfo.InvariantCulture);
                    vel_z[i] = double.Parse(parts[5], CultureInfo.InvariantCulture);
                }
            }
            Console.WriteLine("Система загружена из файла.");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Ошибка загрузки: {ex.Message}");
            // Инициализация по умолчанию при ошибке
            GenerateFCCLattice();
            InitializeVelocities();
        }
        for (int i = 0; i < N; i++)
        {
            acc_x[i] = 0;
            acc_y[i] = 0;
            acc_z[i] = 0;
        }
    }



    // Инициализация системы
    public void InitializeSystem(bool keepArrays = false)
    {

        if (!keepArrays)
        {
            // Инициализация массивов
            pos_x = new double[N];
            pos_y = new double[N];
            pos_z = new double[N];
            vel_x = new double[N];
            vel_y = new double[N];
            vel_z = new double[N];
            acc_x = new double[N];
            acc_y = new double[N];
            acc_z = new double[N];
            e_pot = new double[N];
            disp_list_x = new double[N];
            disp_list_y = new double[N];
            disp_list_z = new double[N];

            v_x = new double[N];
            v_y = new double[N];
            v_z = new double[N];
            v_xo = new double[N];
            v_yo = new double[N];
            v_zo = new double[N];
            b_x = new double[N];
            b_y = new double[N];
            b_z = new double[N];

            max_list_size = 100 * N;
            List = new int[max_list_size];
            Marker1 = new int[N];
            Marker2 = new int[N];
            Advance = new int[N];
        }

        L_x = nx * alat;
        L_y = ny * alat;
        L_z = nz * alat;
        volume = L_x * L_y * L_z;
        density = N / volume;

        GenerateLJPotentialTable();

        // Инициализация параметров баростата
        eta = 0.0;
        eta_dot = 0.0;
        // Инициализация параметроров Нозе-Гувера
        xi = 0.0;
        S_f = 0.0;
    }

    // Генерация ГЦК решетки
    public void GenerateFCCLattice()
    {
        double[,] fcc_basis = {
            {0, 0, 0}, {0.5, 0.5, 0},
            {0.5, 0, 0.5}, {0, 0.5, 0.5}
        };

        int atom_count = 0;
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                for (int z = 0; z < nz; z++)
                    for (int b = 0; b < 4; b++)
                    {
                        if (atom_count >= N) break;
                        pos_x[atom_count] = alat * (x + fcc_basis[b, 0]);
                        pos_y[atom_count] = alat * (y + fcc_basis[b, 1]);
                        pos_z[atom_count] = alat * (z + fcc_basis[b, 2]);
                        atom_count++;
                    }
        CenterSystem();
    }

    // Центрирование системы
    public void CenterSystem()
    {
        double cx = 0, cy = 0, cz = 0;
        for (int i = 0; i < N; i++)
        {
            cx += pos_x[i];
            cy += pos_y[i];
            cz += pos_z[i];
        }
        cx /= N; cy /= N; cz /= N;
        for (int i = 0; i < N; i++)
        {
            pos_x[i] -= cx;
            pos_y[i] -= cy;
            pos_z[i] -= cz;
        }
    }

    // Инициализация скоростей
    public void InitializeVelocities()
    {
        double vx_sum = 0, vy_sum = 0, vz_sum = 0;
        double v2_sum = 0;

        for (int i = 0; i < N; i++)
        {
            vel_x[i] = rand.NextDouble() - 0.5;
            vel_y[i] = rand.NextDouble() - 0.5;
            vel_z[i] = rand.NextDouble() - 0.5;
            vx_sum += vel_x[i];
            vy_sum += vel_y[i];
            vz_sum += vel_z[i];
            v2_sum += vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] + vel_z[i] * vel_z[i];
        }

        // Удаление общего импульса
        vx_sum /= N; vy_sum /= N; vz_sum /= N;
        for (int i = 0; i < N; i++)
        {
            vel_x[i] -= vx_sum;
            vel_y[i] -= vy_sum;
            vel_z[i] -= vz_sum;
            v2_sum += vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] + vel_z[i] * vel_z[i];
        }

        // Масштабирование до начальной температуры
        double temp = v2_sum / (3.0 * N);
        double scale = Math.Sqrt(T0 / temp);
        for (int i = 0; i < N; i++)
        {
            vel_x[i] *= scale;
            vel_y[i] *= scale;
            vel_z[i] *= scale;
        }
    }



    // В конце первого этапа (после цикла моделирования)
    public void SaveSystemState(string filename)
    {
        using (StreamWriter writer = new StreamWriter(filename))
        {
            for (int i = 0; i < N; i++)
            {
                writer.WriteLine(L_x.ToString(CultureInfo.InvariantCulture));
                writer.WriteLine(L_y.ToString(CultureInfo.InvariantCulture));
                writer.WriteLine(L_z.ToString(CultureInfo.InvariantCulture));

                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "{0:F6} {1:F6} {2:F6} {3:F6} {4:F6} {5:F6}",
                    pos_x[i], pos_y[i], pos_z[i],
                    vel_x[i], vel_y[i], vel_z[i]));
            }
        }
        Console.WriteLine($"Сохранено состояние системы в файл: {filename}");
    }


    // Интегрирование уравнений движения
    public void SystemDynamics(int stage_id)
    {
        if (!isNPT) // Режим NVT (Этап I - Плавление)
        {
            if (stage_id == 1)
            {
                // 1. Обновление координат (уравнение движения)
                for (int i = 0; i < N; i++)
                {
                    pos_x[i] += vel_x[i] * dt + acc_x[i] * dt * dt / 2.0;
                    pos_y[i] += vel_y[i] * dt + acc_y[i] * dt * dt / 2.0;
                    pos_z[i] += vel_z[i] * dt + acc_z[i] * dt * dt / 2.0;
                }

                // 2. Расчет кинетической энергии и температуры
                kinetic_energy = 0;
                for (int i = 0; i < N; i++)
                {
                    kinetic_energy += vel_x[i] * vel_x[i]
                                    + vel_y[i] * vel_y[i]
                                    + vel_z[i] * vel_z[i];
                }
                temperature = kinetic_energy / (3.0 * N);  // T = sum(v^2)/(3N)

                // 3. Линейное изменение температуры

                if (temperature < 1e-10) temperature = T0;

                // 4. Вычисление параметра термостата
                lambda = Math.Sqrt(T0 / temperature);

                // 5. Обновление скоростей (первая половина шага) с масштабированием
                for (int i = 0; i < N; i++)
                {
                    vel_x[i] = lambda * vel_x[i] + acc_x[i] * dt / 2.0;
                    vel_y[i] = lambda * vel_y[i] + acc_y[i] * dt / 2.0;
                    vel_z[i] = lambda * vel_z[i] + acc_z[i] * dt / 2.0;
                }
            }
            else if (stage_id == 2)
            {
                // 1. Расчет термодинамических величин
                volume = L_x * L_y * L_z;
                density = N / volume;
                pressure = density * temperature + virial / volume;

                // 2. Обновление скоростей (вторая половина шага)
                for (int i = 0; i < N; i++)
                {
                    vel_x[i] += acc_x[i] * dt / 2.0;
                    vel_y[i] += acc_y[i] * dt / 2.0;
                    vel_z[i] += acc_z[i] * dt / 2.0;
                }

                // 3. Расчет полной энергии
                kinetic_energy = 0;
                for (int i = 0; i < N; i++)
                {
                    kinetic_energy += 0.5 * (vel_x[i] * vel_x[i]
                                      + vel_y[i] * vel_y[i]
                                      + vel_z[i] * vel_z[i]);
                }
                H_energy = kinetic_energy + potential_energy;
            }
        }
        else // Режим NPT (Этап II - Сжатие)
        {
            if (stage_id == 1)
            {
                // 1. Обновление координат с поправкой Нозе-Гувера
                for (int i = 0; i < N; i++)
                {
                    pos_x[i] += vel_x[i] * dt
                              + (acc_x[i] - (lambda + xi) * vel_x[i]) * dt * dt / 2.0;
                    pos_y[i] += vel_y[i] * dt
                              + (acc_y[i] - (lambda + xi) * vel_y[i]) * dt * dt / 2.0;
                    pos_z[i] += vel_z[i] * dt
                              + (acc_z[i] - (lambda + xi) * vel_z[i]) * dt * dt / 2.0;
                }

                // 2. Расчет кинетической энергии и температуры
                kinetic_energy = 0;
                for (int i = 0; i < N; i++)
                {
                    kinetic_energy += vel_x[i] * vel_x[i]
                                    + vel_y[i] * vel_y[i]
                                    + vel_z[i] * vel_z[i];
                }
                temperature = kinetic_energy / (3.0 * N);

                // 3. Обновление параметров термостата
                S_f += lambda * dt + (kinetic_energy - 3.0 * N * T0) * dt * dt / Q_T;
                lambda += (kinetic_energy - 3.0 * N * T0) * dt / (2.0 * Q_T);

                // 4. Обновление скоростей с поправками
                for (int i = 0; i < N; i++)
                {
                    vel_x[i] += (acc_x[i] - (lambda + xi) * vel_x[i]) * dt / 2.0;
                    vel_y[i] += (acc_y[i] - (lambda + xi) * vel_y[i]) * dt / 2.0;
                    vel_z[i] += (acc_z[i] - (lambda + xi) * vel_z[i]) * dt / 2.0;
                }
            }
            else if (stage_id == 2)
            {
                // 1. Обновление баростата
                volume = L_x * L_y * L_z;
                density = N / volume;
                pressure = density * temperature + virial / volume;

                xi += (pressure - P_ext) * dt / (Q_B * density);

                // 2. Масштабирование объема
                L_x *= (1.0 + xi * dt);
                L_y *= (1.0 + xi * dt);
                L_z *= (1.0 + xi * dt);
                volume = L_x * L_y * L_z;
                density = N / volume;

                // Сохраняем старые размеры
                double old_L_x = L_x / (1.0 + xi * dt);
                double old_L_y = L_y / (1.0 + xi * dt);
                double old_L_z = L_z / (1.0 + xi * dt);

                // Масштабируем координаты
                for (int i = 0; i < N; i++)
                {
                    pos_x[i] *= L_x / old_L_x;
                    pos_y[i] *= L_y / old_L_y;
                    pos_z[i] *= L_z / old_L_z;
                }


                // 3. Итерационное уточнение скоростей (упрощенная версия)
                // фрагмент 1.7
                kinetic_energy = 0;
                // Перерасчет скоростей частиц, взаимодействующих с тепловой ванной
                for (int i = 0; i < N; i++)
                {
                    v_x[i] = vel_x[i];
                    v_y[i] = vel_y[i];
                    v_z[i] = vel_z[i];
                    kinetic_energy += v_x[i] * v_x[i] + v_y[i] * v_y[i] + v_z[i] * v_z[i];
                }

                p_s = lambda;
                ready = false;
                iter = 0;
                err = 1.0E-5; // Допустимая погрешность

                while ((!ready) && (iter < 50))
                {
                    iter++;
                    p_o = p_s; // Сохраняем предыдущее значение
                    ds = 0;

                    for (int i = 0; i < N; i++)
                    {
                        // Сохраняем скорости с предыдущей итерации
                        v_xo[i] = v_x[i];
                        v_yo[i] = v_y[i];
                        v_zo[i] = v_z[i];

                        // Вычисляем b_x, b_y, b_z (невязка)
                        b_x[i] = -(dt / 2.0) * (acc_x[i] - p_o * v_x[i]) - (vel_x[i] - v_xo[i]);
                        ds += v_xo[i] * dt * b_x[i] / Q_T;

                        b_y[i] = -(dt / 2.0) * (acc_y[i] - p_o * v_y[i]) - (vel_y[i] - v_yo[i]);
                        ds += v_yo[i] * dt * b_y[i] / Q_T;

                        b_z[i] = -(dt / 2.0) * (acc_z[i] - p_o * v_z[i]) - (vel_z[i] - v_zo[i]);
                        ds += v_zo[i] * dt * b_z[i] / Q_T;
                    }

                    dd = -(p_o * (dt / 2.0) + 1);
                    ds = ds - dd * ((-kinetic_energy + 3.0 * N * T0) * (dt / 2.0) / Q_T - (lambda - p_o));
                    ds = ds / (-dt * (dt / 2.0) * kinetic_energy / Q_T + dd);

                    kinetic_energy = 0;
                    for (int i = 0; i < N; i++)
                    {
                        v_x[i] = v_x[i] + (b_x[i] + (dt / 2.0) * v_xo[i] * ds) / dd;
                        v_y[i] = v_y[i] + (b_y[i] + (dt / 2.0) * v_yo[i] * ds) / dd;
                        v_z[i] = v_z[i] + (b_z[i] + (dt / 2.0) * v_zo[i] * ds) / dd;
                        kinetic_energy += v_x[i] * v_x[i] + v_y[i] * v_y[i] + v_z[i] * v_z[i];
                    }
                    p_s = p_o + ds;

                    // Проверяем сходимость
                    ready = true;
                    int j = 0;
                    while ((j < N) && ready)
                    {
                        if (Math.Abs(v_x[j]) > 1e-10 && Math.Abs((v_x[j] - v_xo[j]) / v_x[j]) > err)
                            ready = false;
                        if (Math.Abs(v_y[j]) > 1e-10 && Math.Abs((v_y[j] - v_yo[j]) / v_y[j]) > err)
                            ready = false;
                        if (Math.Abs(v_z[j]) > 1e-10 && Math.Abs((v_z[j] - v_zo[j]) / v_z[j]) > err)
                            ready = false;
                        j++;
                    }

                    // Также проверяем сходимость для p_s
                    if (Math.Abs(p_s) > 1e-10 && Math.Abs((p_s - p_o) / p_s) > err)
                        ready = false;
                }

                // Обновляем скорости и параметр термостата
                for (int i = 0; i < N; i++)
                {
                    vel_x[i] = v_x[i];
                    vel_y[i] = v_y[i];
                    vel_z[i] = v_z[i];
                }
                lambda = p_s;

                // Пересчет кинетической энергии (теперь с коэффициентом 0.5)
                double kinetic = 0;
                for (int i = 0; i < N; i++)
                {
                    kinetic += 0.5 * (vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] + vel_z[i] * vel_z[i]);
                }

                // Расчет полной энергии
                H_energy = kinetic + potential_energy + (lambda * lambda * Q_T) / 2.0 + 3.0 * N * T0 * S_f;

            }
        }
    }


    // Периодические граничные условия
    public void ApplyBoundaryConditions()
    {
        for (int i = 0; i < N; i++)
        {
            if (pos_x[i] > L_x / 2.0) pos_x[i] -= L_x;
            if (pos_x[i] < -L_x / 2.0) pos_x[i] += L_x;
            if (pos_y[i] > L_y / 2.0) pos_y[i] -= L_y;
            if (pos_y[i] < -L_y / 2.0) pos_y[i] += L_y;
            if (pos_z[i] > L_z / 2.0) pos_z[i] -= L_z;
            if (pos_z[i] < -L_z / 2.0) pos_z[i] += L_z;
        }
    }

    // Генерация таблицы потенциала Леннард-Джонса
    public void GenerateLJPotentialTable()
    {
        table_size = 25000;
        table_pot = new double[table_size, 2];
        double energy_rcut = 4.0 * epsilon * (Math.Pow(sigma / r_cut, 12) - Math.Pow(sigma / r_cut, 6));

        double rij = 0.5;
        double dr = 0.0001;

        for (int i = 0; i < table_size; i++)
        {
            rij += dr;
            double sr6 = Math.Pow(sigma / rij, 6);
            double sr12 = sr6 * sr6;
            double energy = 4.0 * epsilon * (sr12 - sr6) - energy_rcut;
            double force = 24.0 * epsilon * (2 * sr12 - sr6) / rij;
            table_pot[i, 0] = energy;
            table_pot[i, 1] = force;
        }
    }


    // Расчет сил с использованием списка Верле
    public void CalculateForces()
    {
        // Обнуление массивов
        for (int i = 0; i < N; i++)
        {
            acc_x[i] = 0;
            acc_y[i] = 0;
            acc_z[i] = 0;
            e_pot[i] = 0;
        }
        virial = 0;
        potential_energy = 0;

        // Обновление списка Верле при необходимости
        if (update_test)
        {
            UpdateVerletList();
            update_test = false;
        }

        // Расчет сил
        for (int i = 0; i < N - 1; i++)
        {
            for (int j_idx = Marker1[i]; j_idx <= Marker2[i]; j_idx++)
            {
                int j = List[j_idx];
                double dx = pos_x[i] - pos_x[j];
                double dy = pos_y[i] - pos_y[j];
                double dz = pos_z[i] - pos_z[j];

                // Учет периодических граничных условий
                dx -= Math.Round(dx / L_x) * L_x;
                dy -= Math.Round(dy / L_y) * L_y;
                dz -= Math.Round(dz / L_z) * L_z;

                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 < r_cut * r_cut)
                {
                    double r = Math.Sqrt(r2);
                    int idx = Math.Clamp((int)((r - 0.5) / 0.0001), 0, table_size - 1);

                    double energy = table_pot[idx, 0];
                    double force = table_pot[idx, 1];

                    e_pot[i] += energy / 2.0;
                    e_pot[j] += energy / 2.0;
                    virial += force * r;

                    double fx = force * dx / r;
                    double fy = force * dy / r;
                    double fz = force * dz / r;

                    acc_x[i] += fx;
                    acc_y[i] += fy;
                    acc_z[i] += fz;

                    acc_x[j] -= fx;
                    acc_y[j] -= fy;
                    acc_z[j] -= fz;
                }
            }
        }

        // Нормализация ускорений (масса = 1)
        // Ускорения уже рассчитаны правильно

        // Расчет потенциальной энергии
        for (int i = 0; i < N; i++)
            potential_energy += e_pot[i];

        virial /= 3.0;

        // Обновление перемещений для списка Верле
        for (int i = 0; i < N; i++)
        {
            disp_list_x[i] += vel_x[i] * dt + acc_x[i] * dt * dt / 2.0;
            disp_list_y[i] += vel_y[i] * dt + acc_y[i] * dt * dt / 2.0;
            disp_list_z[i] += vel_z[i] * dt + acc_z[i] * dt * dt / 2.0;
        }
        CheckVerletListUpdate();
    }

    // Обновление списка Верле
    public void UpdateVerletList()
    {
        double rng = (r_cut + dr_s) * (r_cut + dr_s);
        int L = 1;

        for (int i = 0; i < N - 1; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                double dx = pos_x[i] - pos_x[j];
                double dy = pos_y[i] - pos_y[j];
                double dz = pos_z[i] - pos_z[j];

                dx -= Math.Round(dx / L_x) * L_x;
                dy -= Math.Round(dy / L_y) * L_y;
                dz -= Math.Round(dz / L_z) * L_z;

                double r2 = dx * dx + dy * dy + dz * dz;
                Advance[j] = (r2 < rng) ? 1 : 0;

            }

            Marker1[i] = L;
            for (int j = i + 1; j < N; j++)
            {
                if (Advance[j] == 1)
                {
                    if (L >= max_list_size)
                    {
                        Console.WriteLine("Ошибка: список Верле переполнен");
                        break;
                    }
                    List[L] = j;
                    L++;
                }
            }
            Marker2[i] = L - 1;

        }
        list_size = L - 1;

        // Сброс перемещений
        for (int i = 0; i < N; i++)
        {
            disp_list_x[i] = 0;
            disp_list_y[i] = 0;
            disp_list_z[i] = 0;
        }
    }

    // Проверка необходимости обновления списка Верле
    public void CheckVerletListUpdate()
    {
        double max_disp1 = 0, max_disp2 = 0;
        for (int i = 0; i < N; i++)
        {
            double disp = Math.Sqrt(
                disp_list_x[i] * disp_list_x[i] +
                disp_list_y[i] * disp_list_y[i] +
                disp_list_z[i] * disp_list_z[i]);

            if (disp > max_disp1)
            {
                max_disp2 = max_disp1;
                max_disp1 = disp;
            }
            else if (disp > max_disp2)
            {
                max_disp2 = disp;
            }
        }
        update_test = (max_disp1 + max_disp2 > dr_s);
    }
    // Сохранение траектории в формате LAMMPS (с правильным форматом чисел)
    public void SaveLAMMPSTrajectory(int step)
    {
        using (StreamWriter writer = new StreamWriter("system_config.lammpstrj", true))
        {
            writer.WriteLine("ITEM: TIMESTEP");
            writer.WriteLine(step);
            writer.WriteLine("ITEM: NUMBER OF ATOMS");
            writer.WriteLine(N);
            writer.WriteLine("ITEM: BOX BOUNDS pp pp pp");

            // Форматирование чисел с точкой вместо запятой
            writer.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0} {1}", -L_x / 2.0, L_x / 2.0));
            writer.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0} {1}", -L_y / 2.0, L_y / 2.0));
            writer.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0} {1}", -L_z / 2.0, L_z / 2.0));

            writer.WriteLine("ITEM: ATOMS id type x y z");

            for (int i = 0; i < N; i++)
            {
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "{0} 1 {1:F6} {2:F6} {3:F6}",
                    i + 1, pos_x[i], pos_y[i], pos_z[i]));
            }
        }
    }

    // Сохранение термодинамических параметров
    public void SaveThermodynamicParameters()
    {
        // Раздельные файлы для каждого параметра как в методичке
        SaveParameter("system_temperature.txt", time, temperature);
        SaveParameter("system_pressure.txt", time, pressure);
        SaveParameter("system_kinetic_energy.txt", time, kinetic_energy);
        SaveParameter("system_potential_energy.txt", time, potential_energy);
        SaveParameter("system_total_energy.txt", time, H_energy);
    }
    private void SaveParameter(string filename, double time, double value)
    {
        try
        {
            using (StreamWriter writer = new StreamWriter(filename, true))
            {
                // Форматирование строки как в методичке: время + значение
                string line = $"{time.ToString(CultureInfo.InvariantCulture)} {value.ToString(CultureInfo.InvariantCulture)}";
                writer.WriteLine(line);
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Ошибка сохранения параметра {filename}: {ex.Message}");
        }
    }

    // Расчет функции радиального распределения
    public void CalculateRDF(string phaseName, int begin_step = 1, int end_step = 1000)
    {
        // Параметры из методички (фрагмент 4.3)
        double dr = 0.01;
        int num_bins = 500;

        // Инициализация массивов как в методичке
        double[] g_sum = new double[num_bins];
        int frame_count = 0;
        int max_bin = 0;
        double total_volume = 0; //Сумма объемов всех кадров

        // Чтение координат
        double[] x = new double[N];
        double[] y = new double[N];
        double[] z = new double[N];

        try
        {
            using (var reader = new StreamReader("system_config.lammpstrj")) // Имя из методички
            {
                while (!reader.EndOfStream)
                {
                    // Чтение шага
                    reader.ReadLine(); // "ITEM: TIMESTEP"
                    string? stepLine = reader.ReadLine();
                    if (stepLine == null) break;
                    int current_step = int.Parse(stepLine);

                    // Пропускаем кадры до begin_step
                    if (current_step < begin_step || current_step > end_step)
                    {
                        // Пропускаем 8 строк кадра
                        SkipFrame(reader);
                        continue;

                    }

                    // Чтение количества атомов
                    reader.ReadLine(); // "ITEM: NUMBER OF ATOMS"
                    string? atomsLine = reader.ReadLine();
                    if (atomsLine == null) break;
                    int num_atoms = int.Parse(atomsLine);

                    // Чтение размеров ячейки
                    reader.ReadLine(); // "ITEM: BOX BOUNDS pp pp pp"
                    string? boxX = reader.ReadLine();
                    string? boxY = reader.ReadLine();
                    string? boxZ = reader.ReadLine();

                    if (boxX == null || boxY == null || boxZ == null) break;

                    double Lx = double.Parse(boxX.Split()[1], CultureInfo.InvariantCulture)
                            - double.Parse(boxX.Split()[0], CultureInfo.InvariantCulture);

                    double Ly = double.Parse(boxY.Split()[1], CultureInfo.InvariantCulture)
                            - double.Parse(boxY.Split()[0], CultureInfo.InvariantCulture);

                    double Lz = double.Parse(boxZ.Split()[1], CultureInfo.InvariantCulture)
                            - double.Parse(boxZ.Split()[0], CultureInfo.InvariantCulture);

                    // Максимальный радиус - половина минимального размера
                    double max_r = Math.Min(Math.Min(Lx, Ly), Lz) * 0.5;

                    // Пропуск заголовка атомов
                    reader.ReadLine(); // "ITEM: ATOMS id type x y z"

                    for (int i = 0; i < num_atoms; i++)
                    {
                        string? atomLine = reader.ReadLine();
                        if (atomLine == null) break;

                        string[] parts = atomLine.Split();
                        if (parts.Length < 5) continue;

                        x[i] = double.Parse(parts[2], CultureInfo.InvariantCulture);
                        y[i] = double.Parse(parts[3], CultureInfo.InvariantCulture);
                        z[i] = double.Parse(parts[4], CultureInfo.InvariantCulture);
                    }

                    double volume_frame = Lx * Ly * Lz;
                    total_volume += volume_frame; // Важно: накапливаем объем

                    // Проверка количества атомов
                    if (num_atoms != N)
                    {
                        Console.WriteLine($"Ошибка: количество атомов изменилось ({num_atoms} vs {N})");
                        return;
                    }

                    // Расчет парных расстояний (как в методичке)
                    for (int i = 0; i < num_atoms - 1; i++)
                    {
                        for (int j = i + 1; j < num_atoms; j++)
                        {
                            double dx = x[i] - x[j];
                            double dy = y[i] - y[j];
                            double dz = z[i] - z[j];

                            // Учет PBC (методичка 1.2)
                            dx -= Math.Round(dx / Lx) * Lx;
                            dy -= Math.Round(dy / Ly) * Ly;
                            dz -= Math.Round(dz / Lz) * Lz;

                            double r = Math.Sqrt(dx * dx + dy * dy + dz * dz);

                            if (r < max_r)
                            {
                                int bin = (int)(r / dr);
                                if (bin < num_bins)
                                {
                                    g_sum[bin]++;
                                    if (bin > max_bin) max_bin = bin;
                                }
                            }
                        }
                    }

                    frame_count++;
                    if (frame_count % 100 == 0)
                        Console.WriteLine($"Обработано кадров: {frame_count}");
                }
            }

            // Проверка наличия кадров
            if (frame_count == 0)
            {
                Console.WriteLine($"Предупреждение: для фазы {phaseName} не найдено ни одного кадра в диапазоне {begin_step}-{end_step}");
                return;
            }

            // выбор нормаированного объема в зависимости от фазы
            double norm_volume;
            if (phaseName == "liquid")
            {
                // Для жидкости (NVT) - постоянный объем
                norm_volume = L_x * L_y * L_z;
            }
            else // "compressed"
            {
                // Для сжатой фазы (NPT) - средний объем
                norm_volume = total_volume / frame_count;
            }

            // Формируем имя файла на основе фазы
            string filename = $"rdf_{phaseName}.csv";

            using (StreamWriter writer = new StreamWriter(filename))
            {
                writer.WriteLine("r,g_r");
                for (int i = 0; i <= max_bin; i++)
                {
                    double r_inner = i * dr;
                    double r_outer = (i + 1) * dr;
                    double shell_volume = 4.0 / 3.0 * Math.PI *
                        (Math.Pow(r_outer, 3) - Math.Pow(r_inner, 3));

                    double r = (i + 0.5) * dr;
                    double g_r = g_sum[i] / (frame_count * 0.5 * N * (N - 1) * shell_volume / norm_volume);
                    writer.WriteLine($"{r.ToString(CultureInfo.InvariantCulture)},{g_r.ToString(CultureInfo.InvariantCulture)}");
                }
            }
            Console.WriteLine($"Расчет RDF завершен. Кадров: {frame_count}, Макс. bin: {max_bin}");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Ошибка при расчете RDF: {ex.Message}");
        }
    }
    private void SkipFrame(StreamReader reader)
    {
        reader.ReadLine(); // "ITEM: NUMBER OF ATOMS"
        int num_atoms = int.Parse(reader.ReadLine() ?? "0");
        for (int i = 0; i < 3; i++) reader.ReadLine(); // BOX BOUNDS
        reader.ReadLine(); // "ITEM: ATOMS"
        for (int i = 0; i < num_atoms; i++) reader.ReadLine();
    }
}   