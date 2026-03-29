# Проект по фурье-фильтрации

Этот проект решает задачу выделения полезного сигнала из зашумлённого бинарного файла измерений с помощью преобразования Фурье.

## Что делает программа

Программа:
- считывает бинарный файл с `100000` отсчётами;
- строит амплитудный спектр;
- оценивает частоты полезного сигнала;
- восстанавливает очищенный сигнал;
- сохраняет графики и численные результаты.

## Структура проекта

```text
fourier_filtering_project/
├── docs/
│   ├── PLAN.md
│   ├── SOLUTION.md
├── examples/
│   ├── demo_output/
│   └── demo_signal_float32.bin
├── scripts/
│   └── generate_demo_signal.py
├── src/
│   └── fourier_filter/
│       ├── __init__.py
│       └── pipeline.py
├── .gitignore
├── main.py
├── README.md
├── requirements.txt
├── run_demo.ps1
└── run_demo.sh
```

## Установка

### Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Linux / macOS

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Запуск на реальном файле варианта

Если у тебя уже есть бинарный файл варианта, например `variant_5.bin`, запусти:

```bash
python main.py --input variant_5.bin --output output
```

Если тип данных известен заранее, можно указать его явно:

```bash
python main.py --input variant_5.bin --output output --dtype float32
```

## Основные аргументы программы

- `--input` — путь к бинарному файлу сигнала;
- `--output` — папка для сохранения результатов;
- `--samples` — ожидаемое число отсчётов, по умолчанию `100000`;
- `--dt-ms` — шаг дискретизации в миллисекундах, по умолчанию `0.1`;
- `--dtype` — тип данных: `auto`, `float32`, `float64`, `int16`, `int32`, `int64`;
- `--endianness` — порядок байтов: `little` или `big`;
- `--min-peak-ratio` — порог по отношению к уровню шума, по умолчанию `2.0`;
- `--neighbor-bins` — число сохраняемых бинов вокруг каждого пика, по умолчанию `2`;
- `--top-k` — максимальное число выводимых пиков, по умолчанию `20`.

## Какие файлы создаются после запуска

После выполнения программы формируются:
- `detected_frequencies.txt`;
- `summary.json`;
- `filtered_signal.csv`;
- `spectrum.csv`;
- `signal_time_domain.png`;
- `filtered_signal.png`;
- `spectrum.png`.

## Демонстрационный пример

В проекте есть генератор тестового сигнала. Ниже собраны все основные команды, которые удобно использовать для демонстрации работы программы.

### 1. Сгенерировать демо-сигнал

```bash
python scripts/generate_demo_signal.py
```

После этого появится файл:

```text
examples/demo_signal_float32.bin
```

### 2. Базовый запуск демо-примера

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output --dtype float32
```

### 3. Запуск через готовые скрипты

#### Windows PowerShell

```powershell
./run_demo.ps1
```

Если PowerShell блокирует выполнение локальных скриптов:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_demo.ps1
```

#### Linux / macOS

```bash
chmod +x run_demo.sh
./run_demo.sh
```

### 4. Просмотр всех доступных аргументов программы

```bash
python main.py --help
```

### 5. Демо с другим именем выходной папки

```bash
python main.py --input examples/demo_signal_float32.bin --output demo_result --dtype float32
```

### 6. Демо без окна Ханна

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output_no_hann --dtype float32 --no-hann-window
```

### 7. Демо с ручной настройкой порога поиска пиков

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output_ratio_25 --dtype float32 --min-peak-ratio 2.5
```

### 8. Демо с изменением числа соседних бинов вокруг пика

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output_bins_4 --dtype float32 --neighbor-bins 4
```

### 9. Демо с ограничением числа выводимых пиков

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output_top3 --dtype float32 --top-k 3
```

### 10. Демо с явным указанием числа отсчётов и шага дискретизации

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output_manual --dtype float32 --samples 100000 --dt-ms 0.1
```

### 11. Демо с указанием порядка байтов

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output_little --dtype float32 --endianness little
```

Для big-endian:

```bash
python main.py --input examples/demo_signal_float32.bin --output examples/demo_output_big --dtype float32 --endianness big
```

### 12. Какие результаты показать после Демоа

После выполнения удобно открыть:
- `examples/demo_output/detected_frequencies.txt`;
- `examples/demo_output/spectrum.png`;
- `examples/demo_output/filtered_signal.png`;
- `examples/demo_output/signal_time_domain.png`.

Для демонстрационного сигнала ожидаются основные частоты примерно:
- `50 Hz`;
- `120 Hz`;
- `400 Hz`.