# coding: utf-8

from __future__ import annotations

__author__ = "Giyeok Lee, Joonsuk Cheung"
__email__ = "lgy4230@yonsei.ac.kr"
__date__ = "26/Feb/2026"
__maintainer__ = "Giyeok Lee (@lgyEthan)"
__copyright__ = "Copyright (c) Materials Theory Group @ Yonsei University (2026)"

import argparse
import glob
import os
import sys
from collections.abc import Iterable
from pathlib import Path

from autobskan import __version__

try:
    from tqdm import tqdm
except Exception:
    # Keep CLI usable even if tqdm is unavailable in constrained environments.
    def tqdm(iterable, **_kwargs):
        return iterable


def _resolve_current_files(current_option):
    if isinstance(current_option, str):
        pattern = current_option.strip()
        if pattern.upper() == "ALL":
            pattern = "*"
        return sorted(glob.glob(pattern))

    if isinstance(current_option, Iterable):
        return [str(path) for path in current_option]

    raise TypeError(f"CURRENT option has unsupported type: {type(current_option)}")


def _run_image_task(bskan_input):
    from autobskan.image import AR, serial_run, stmplot

    input_source = bskan_input.option_dict.get("INPUT_SOURCE", "BSKAN")
    simulation = bskan_input.option_dict.get("SIMULATION", "STM")
    volume_option = (
        bskan_input.option_dict.get("VOLUME")
        if input_source == "VASP"
        else bskan_input.option_dict["CURRENT"]
    )
    volume_files = _resolve_current_files(volume_option)
    if not volume_files:
        raise FileNotFoundError(
            "No volumetric files were found from CURRENT/VOLUME option in bskan.in."
        )

    completed = 0
    skipped = 0
    failed = 0
    generated_images = 0
    for volume_file in tqdm(volume_files, desc=f"Total {len(volume_files)} files"):
        if not os.path.isfile(volume_file):
            print(f"[WARN] Skip non-file path: {volume_file}")
            skipped += 1
            continue
        if os.path.getsize(volume_file) == 0:
            print(f"[WARN] Skip empty volumetric file: {volume_file}")
            skipped += 1
            continue
        try:
            outputs = None
            if input_source == "VASP":
                volumetric = AR.Chgcar(volume_file)
                outputs = serial_run.single_vasp_volume(
                    volumetric=volumetric,
                    bskan_input=bskan_input,
                    volume_path=volume_file,
                    image_dir=f"{volume_file}_images",
                )
            else:
                current = stmplot.Current(volume_file)
                if simulation == "PHI_APP":
                    outputs = serial_run.single_apparent_barrier_current(
                        current=current,
                        bskan_input=bskan_input,
                        image_dir=f"{volume_file}_images",
                    )
                elif simulation == "LWF":
                    raise ValueError(
                        "Phi_loc requires INPUT_SOURCE = VASP and a LOCPOT file."
                    )
                else:
                    outputs = serial_run.single_current(
                        current=current,
                        bskan_input=bskan_input,
                        image_dir=f"{volume_file}_images",
                    )
            completed += 1
            if isinstance(outputs, (list, tuple)):
                generated_images += len(outputs)
        except Exception as exc:
            print(f"[WARN] Failed for {volume_file}: {exc}")
            failed += 1

    print(
        "Image task summary: "
        f"{completed} completed, {skipped} skipped, {failed} failed"
        + (f", {generated_images} STM image(s) generated" if generated_images else "")
        + "."
    )
    if completed == 0:
        raise RuntimeError("No volumetric file produced an image.")
    return {
        "completed": completed,
        "skipped": skipped,
        "failed": failed,
        "generated_images": generated_images,
    }


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="autobskan CLI: image generation and bSKAN calculation workflow."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "-i",
        "--input",
        default="bskan.in",
        help="Input file path. Default: bskan.in",
    )
    subparsers = parser.add_subparsers(dest="command")
    precheck = subparsers.add_parser(
        "pre-check",
        aliases=["precheck"],
        help="Check and optionally correct the structure origin before SCF.",
    )
    precheck.add_argument(
        "structure",
        nargs="?",
        default="POSCAR",
        help="VASP POSCAR/CONTCAR to check. Default: POSCAR",
    )
    precheck.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Show each pre-check stage (default).",
    )
    precheck.add_argument(
        "-q",
        "--quiet",
        dest="verbose",
        action="store_false",
        help="Show only the result and output-file decision.",
    )
    precheck.set_defaults(verbose=True)
    return parser.parse_args(argv)


def _vector_text(vector):
    return " ".join(f"{float(value): .9f}" for value in vector)


def _default_shifted_path(structure):
    source = Path(os.path.abspath(os.path.expanduser(structure)))
    stem = source.stem if source.suffix else source.name
    return source.with_name(f"{stem}-autobskan_shifted.vasp")


def _show_precheck_stage(args, number, message):
    if args.verbose:
        print(f"[{number}/4] {message}")


def _run_precheck(args):
    from autobskan.calculation.symmetry_precheck import (
        PrecheckStatus,
        precheck_bskan_structure,
        write_bskan_asample,
        write_prechecked_poscar,
    )

    source = Path(os.path.abspath(os.path.expanduser(args.structure)))
    output = _default_shifted_path(source)
    asample = source.with_name("ASAMPLE")
    _show_precheck_stage(args, 1, f"Reading VASP structure: {source}")
    result = precheck_bskan_structure(source)
    if result.status is PrecheckStatus.UNSUPPORTED_GEOMETRY:
        _show_precheck_stage(
            args,
            2,
            "Checked lattice geometry; the c-axis surface-normal prerequisite failed.",
        )
    else:
        _show_precheck_stage(
            args,
            2,
            "Completed the bSKAN-compatible lattice and atomic-basis symmetry search.",
        )
    if args.verbose:
        print("AutoBSKAN bSKAN symmetry pre-check")
        print(f"Structure: {source}")
    print(f"Status: {result.status.value}")
    if args.verbose and result.status is PrecheckStatus.UNSUPPORTED_GEOMETRY:
        print("bSKAN basis operations: not evaluated")
    elif args.verbose:
        print(
            "bSKAN basis operations: "
            f"{result.operation_count}; nonzero translations: "
            f"{result.nonzero_operation_count}"
        )

    if result.status is PrecheckStatus.UNSUPPORTED_GEOMETRY:
        _show_precheck_stage(
            args,
            3,
            "Cannot guarantee bSKAN symmetry safety for a tilted c lattice vector.",
        )
        _show_precheck_stage(args, 4, "Skipped all output-file generation.")
        print(
            "Pre-check error: c must be perpendicular to the ab plane "
            "(c parallel to a x b)."
        )
        print(f"Reason: {result.unsafe_reason}")
        print(
            "Construct an equivalent surface cell satisfying this geometry, "
            "then run pre-check again before SCF."
        )
        print(f"No shifted file was generated: {output}")
        print(f"No ASAMPLE was generated or modified: {asample}")
        return 2

    if result.status is PrecheckStatus.READY:
        _show_precheck_stage(args, 3, "No origin translation is required.")
        asample_existed = asample.exists()
        written_asample = write_bskan_asample(source, asample, force=True)
        _show_precheck_stage(
            args,
            4,
            "Wrote and validated VASP 4 ASAMPLE with Selective dynamics and Direct coordinates.",
        )
        print(
            "Shift not required: all detected symmetry translations are zero "
            f"within {result.tolerance:.1e}."
        )
        print(f"No shifted file was generated: {output}")
        if output.exists():
            print("A file already present at that path was left unchanged.")
        action = "Updated" if asample_existed else "Generated"
        print(f"{action} ASAMPLE: {written_asample}")
        return 0

    if result.status is PrecheckStatus.UNSAFE:
        _show_precheck_stage(args, 3, "No valid common origin was found.")
        _show_precheck_stage(args, 4, "Skipped shifted-structure generation.")
        print(
            "Simple shift cannot resolve this structure: no single origin makes "
            "every detected bSKAN translation zero."
        )
        print(f"Reason: {result.unsafe_reason}")
        if args.verbose:
            print("Nonzero operations (fractional coordinates):")
            for index in result.nonzero_operation_indices:
                operation = result.operations[index]
                print(
                    f"  #{index + 1}: t={_vector_text(operation.translation)}; "
                    f"R={operation.rotation.tolist()}"
                )
        print(
            "Do not start the bSKAN workflow with this structure. Use a "
            "symmetry-free/full-BZ workflow or correct bSKAN's symmetry expansion."
        )
        print(f"No shifted file was generated: {output}")
        if output.exists():
            print("A file already present at that path was left unchanged.")
        print(f"No ASAMPLE was generated or modified: {asample}")
        if asample.exists():
            print("The existing ASAMPLE was left unchanged and must not be used for this structure.")
        return 3

    if args.verbose:
        print(f"Origin o (fractional): {_vector_text(result.origin_fractional)}")
        print(
            "ASE displacement -o*cell (angstrom): "
            f"{_vector_text(result.ase_displacement_angstrom)}"
        )
    _show_precheck_stage(
        args,
        3,
        "Applying the proposed displacement with ASE Atoms.translate().",
    )
    output_existed = output.exists()
    written, _before, after = write_prechecked_poscar(
        source,
        output,
        force=True,
    )
    _show_precheck_stage(
        args,
        4,
        "Re-read the translated structure; post-check status is READY.",
    )
    action = "Updated" if output_existed else "Created"
    print(
        f"Shift applied and verified: {action.lower()} {written} with "
        f"{after.nonzero_operation_count} nonzero translations after post-check."
    )
    print(f"Generated file: {written}")
    print(f"No ASAMPLE was generated or modified: {asample}")
    if asample.exists():
        print(
            "The existing ASAMPLE was left unchanged; regenerate it only after "
            "the new SCF is complete."
        )
    print("Use this file as POSCAR and run SCF from the beginning.")
    return 0


def main(argv=None):
    args = _parse_args(argv)
    if args.command in {"pre-check", "precheck"}:
        return _run_precheck(args)

    from autobskan.io.input import Options

    if os.path.exists(args.input):
        bskan_input = Options(args.input)
    else:
        print(f"[WARN] Cannot find {args.input}. Falling back to default options.")
        bskan_input = Options(None)

    task = bskan_input.option_dict["TASK"]
    if task == "CALCULATION":
        from autobskan.calculation import connect
        from autobskan.calculation.symmetry_precheck import (
            ExistingScfRequiresRestartError,
        )

        if not os.path.exists(args.input):
            raise FileNotFoundError(
                "TASK=CALCULATION requires an explicit input file (e.g. bskan.in)."
            )
        try:
            connect.main(args.input)
        except ExistingScfRequiresRestartError as exc:
            print(f"[ERROR] {exc}", file=sys.stderr)
            return 4
        return 0

    if task == "IMAGE":
        try:
            result = _run_image_task(bskan_input)
        except Exception as exc:
            print(f"[ERROR] Image generation failed: {exc}", file=sys.stderr)
            return 5
        return 5 if result["failed"] else 0

    if task == "POST_PROCESSING_ONLY":
        from autobskan.image import post_processing

        candidate = bskan_input.option_dict.get("CURRENT", "raw_images")
        raw_image_dir = candidate if os.path.isdir(candidate) else "raw_images"
        if not os.path.isdir(raw_image_dir):
            raise FileNotFoundError(
                "Cannot find raw image directory for POST_PROCESSING_ONLY. "
                "Set CURRENT=<raw_image_dir> or prepare ./raw_images."
            )
        post_processing.main(raw_image_dir, bskan_input)
        return 0

    raise IOError(
        "TASK in autobskan should be one of CALCULATION / IMAGE / "
        f"POST_PROCESSING_ONLY, not {task}"
    )


if __name__ == "__main__":
    raise SystemExit(main())
