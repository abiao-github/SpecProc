#!/usr/bin/env python3
"""
SpecProc Application Launcher

PyQt-based GUI or CLI for echelle spectrograph FITS data reduction.

Note: If you encounter import errors, please use the proper launch method:
- For GUI:  python3 run.py --mode gui
- For CLI:  python3 run.py --mode cli

Or install in development mode:
  pip install -e .

Make sure SpecProc is installed in your Python environment.
"""

import argparse
import os
from pathlib import Path
from typing import List

# Add src directory to path
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from src.config.config_manager import ConfigManager
from src.core.processing_pipeline import ProcessingPipeline
from src.utils.fits_io import read_fits_image


def parse_selection(selection: str, max_index: int) -> List[int]:
    """Parse selection string like '1,3,5-7' into zero-based indices."""
    indices = set()
    for part in selection.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            try:
                start, end = part.split('-', 1)
                i0 = int(start) - 1
                i1 = int(end) - 1
                for i in range(min(i0, i1), max(i0, i1) + 1):
                    if 0 <= i < max_index:
                        indices.add(i)
            except ValueError:
                continue
        else:
            try:
                i = int(part) - 1
                if 0 <= i < max_index:
                    indices.add(i)
            except ValueError:
                continue
    return sorted(indices)


def choose_files(files: List[Path], prompt: str) -> List[Path]:
    if not files:
        return []
    print(prompt)
    for idx, f in enumerate(files, start=1):
        print(f"  [{idx}] {f.name}")

    selection = input("输入序号/范围 (例如 1,3-5)：").strip()
    chosen = parse_selection(selection, len(files))
    return [files[i] for i in chosen]


def run_cli(config_path: str = None):
    config = ConfigManager(config_path)

    rawdata_path = config.get_rawdata_path()
    if not Path(rawdata_path).exists():
        rawdata_path = input(f"当前 rawdata_path 路径 '{rawdata_path}' 不存在，请输入 RawData 路径或回车使用默认：").strip() or rawdata_path
        # Store the path as entered (without expanding for config file)
        config.set('data', 'rawdata_path', rawdata_path)
        config.save()
        # Re-read with expansion for use
        rawdata_path = Path(rawdata_path).expanduser().as_posix()
    Path(rawdata_path).mkdir(parents=True, exist_ok=True)

    files = sorted(Path(rawdata_path).glob("*.fits"))
    if not files:
        print('没有发现 FITS 文件。请检查 rawdata 路径。')
        return

    bias_files = choose_files(files, '请选择 Bias 文件：')
    flat_files = choose_files(files, '请选择 Flat 文件：')
    science_files = choose_files(files, '请选择 Science 文件：')
    calib_files = choose_files(files, '请选择 Wavelength Calibration 文件（选 1 个即可）：')

    if not science_files:
        print('必须选择至少一个 Science 文件，程序退出。')
        return
    if not bias_files:
        print('建议选择 Bias 文件。')
    if not flat_files:
        print('建议选择 Flat 文件。')
    if not calib_files:
        print('建议选择 Calib 文件。')

    stages_input = input('选择处理步骤 (0-5，逗号/范围，回车=全步): ').strip() or '0-5'
    stage_indices = parse_selection(stages_input, 6)

    pipeline = ProcessingPipeline(config)

    for sci in science_files:
        try:
            print('\n=== 处理 Science 文件:', sci, '===')

            current_science = sci
            if 0 in stage_indices:
                all_inputs = [str(current_science)] + [str(p) for p in bias_files] + [str(p) for p in flat_files]
                if calib_files:
                    all_inputs.append(str(calib_files[0]))
                corrected = pipeline.stage_overscan_correction(all_inputs)
                current_science = Path(corrected[0])
                if bias_files:
                    bias_files = [Path(x) for x in corrected[1: 1 + len(bias_files)]]
                if flat_files:
                    flat_files = [Path(x) for x in corrected[1+len(bias_files): 1+len(bias_files)+len(flat_files)]]
                if calib_files:
                    calib_files = [Path(corrected[-1])]

            master_bias = None
            if 1 in stage_indices and bias_files:
                master_bias = pipeline.stage_bias_correction([str(p) for p in bias_files])

            apertures = None
            if 2 in stage_indices and flat_files:
                _, apertures = pipeline.stage_flat_fielding([str(p) for p in flat_files])

            wave_calib = None
            if 3 in stage_indices and calib_files:
                wave_calib = pipeline.stage_wavelength_calibration(str(calib_files[0]))

            science_image, header = read_fits_image(str(current_science))
            if master_bias is not None:
                science_image = science_image.astype(float) - master_bias

            if 4 in stage_indices:
                science_image = pipeline.stage_scattered_light_subtraction_science(science_image)

            if 5 in stage_indices:
                if apertures is None:
                    raise RuntimeError('Stage 5 需要先执行 Stage 2 以得到 apertures')
                if wave_calib is None:
                    raise RuntimeError('Stage 5 需要先执行 Stage 3 以得到 wavelength calibration')
                pipeline.state.apertures = apertures
                pipeline.state.wavelength_calib = wave_calib
                pipeline.stage_extraction(science_image)

            print(f'文件 {sci.name} 处理完成')

        except Exception as e:
            print(f'处理失败: {e}')

    print('\nCLI 处理结束。')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SpecProc: GUI and CLI spectral reduction')
    parser.add_argument('--mode', choices=['gui', 'cli'], default='gui', help='运行模式')
    parser.add_argument('--config', default=None, help='配置文件路径')
    args = parser.parse_args()

    if args.mode == 'cli':
        run_cli(args.config)
    else:
        from src.main import main
        main(args.config)
