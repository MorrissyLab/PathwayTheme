"""Command-line interface: ``pathwaytheme run config.yaml``."""

from __future__ import annotations

import argparse
import sys

from .config import PipelineConfig


def _cmd_run(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_yaml(args.config)
    if args.output_dir:
        config.output_dir = args.output_dir
    if args.no_figures:
        config.viz.make_figures = False
    from .pipeline import run
    run(config, verbose=not args.quiet)
    return 0


def _cmd_init(args: argparse.Namespace) -> int:
    """Write a default config YAML to stdout or a file."""
    cfg = PipelineConfig()
    if args.output:
        cfg.to_yaml(args.output)
        print(f"wrote default config to {args.output}")
    else:
        import yaml
        print(yaml.safe_dump(cfg.to_dict(), sort_keys=False))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pathwaytheme",
        description="Omic enrichment (ssGSEA/EnrichR/GoSlim) + PCA pathway analysis.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="run the full pipeline from a YAML config")
    p_run.add_argument("config", help="path to a pipeline config YAML")
    p_run.add_argument("--output-dir", default=None, help="override output_dir")
    p_run.add_argument("--no-figures", action="store_true", help="tables only")
    p_run.add_argument("--quiet", action="store_true")
    p_run.set_defaults(func=_cmd_run)

    p_init = sub.add_parser("init", help="print / write a default config YAML")
    p_init.add_argument("-o", "--output", default=None, help="write to file instead of stdout")
    p_init.set_defaults(func=_cmd_init)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
