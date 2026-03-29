\
#!/usr/bin/env bash
set -euo pipefail

python3 scripts/generate_demo_signal.py
python3 main.py --input examples/demo_signal_float32.bin --output examples/demo_output --dtype float32
