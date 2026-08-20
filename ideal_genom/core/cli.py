#!/usr/bin/env python3
"""
Command-line interface for IDEAL-GENOM-QC pipeline.

This module provides the main CLI entry point for running genomic quality control
pipelines using YAML configuration files.
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Optional

# Package imports
from .config import load_config, ConfigurationError
from .pipeline import PipelineExecutor

# Try to get version from package
try:
    from importlib.metadata import version
    __version__ = version("ideal-genom")
except Exception:
    __version__ = "1.3.0"  # fallback version


def setup_logging(level: str = "INFO") -> None:
    """
    Setup basic logging configuration.
    
    Parameters
    ----------
    level : str
        Logging level (DEBUG, INFO, WARNING, ERROR)
    """
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {level}')
    
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def validate_config_file(config_path: str) -> Path:
    """
    Validate that the configuration file exists and is readable.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file
        
    Returns
    -------
    Path
        Validated configuration file path
        
    Raises
    ------
    FileNotFoundError
        If configuration file doesn't exist
    """
    config_file = Path(config_path)
    
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    if not config_file.is_file():
        raise FileNotFoundError(f"Configuration path is not a file: {config_path}")
    
    return config_file


def cmd_run(args: argparse.Namespace) -> int:
    """
    Execute the run command.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments
        
    Returns
    -------
    int
        Exit code (0 for success, 1 for failure)
    """
    try:
        # Validate configuration file
        config_file = validate_config_file(args.config)
        
        # Load and validate configuration
        config = load_config(str(config_file))
        
        # Setup logging based on configuration or CLI args
        log_level = args.log_level or config.get('settings', {}).get('logging', {}).get('level', 'INFO')
        setup_logging(log_level)
        
        logger = logging.getLogger('ideal_genom.cli')
        logger.info(f"Loading configuration from: {config_file}")
        
        # Create and configure pipeline executor
        executor = PipelineExecutor(config, dry_run=args.dry_run)
        
        # Show pipeline summary if requested
        if args.dry_run:
            summary = executor.get_pipeline_summary()
            print("\n" + "="*60)
            print("PIPELINE SUMMARY (DRY RUN)")
            print("="*60)
            print(f"Pipeline Name: {summary['pipeline_name']}")
            print(f"Output Directory: {summary['base_output_dir']}")
            print(f"Total Steps: {summary['total_steps']}")
            print(f"Enabled Steps: {summary['enabled_steps']}")
            print("\nEnabled Steps:")
            for i, step in enumerate(summary['steps'], 1):
                print(f"  {i}. {step['name']} ({step['module']}.{step['class']})")
            print("="*60)
            return 0
        
        # Execute pipeline
        logger.info("Starting pipeline execution")
        executor.execute()
        logger.info("Pipeline execution completed successfully")
        
        return 0
        
    except ConfigurationError as e:
        print(f"Configuration Error: {e}", file=sys.stderr)
        return 1
    except FileNotFoundError as e:
        print(f"File Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Pipeline Error: {e}", file=sys.stderr)
        logging.getLogger('ideal_genom.cli').error(f"Unexpected error: {e}", exc_info=True)
        return 1


def cmd_validate(args: argparse.Namespace) -> int:
    """
    Execute the validate command.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments
        
    Returns
    -------
    int
        Exit code (0 for success, 1 for failure)
    """
    try:
        # Validate configuration file
        config_file = validate_config_file(args.config)
        
        # Load and validate configuration
        config = load_config(str(config_file))
        
        # Create executor to validate pipeline configuration
        executor = PipelineExecutor(config, dry_run=True)
        summary = executor.get_pipeline_summary()
        
        print("✓ Configuration file is valid")
        print(f"✓ Pipeline '{summary['pipeline_name']}' configured with {summary['enabled_steps']}/{summary['total_steps']} enabled steps")
        
        return 0
        
    except ConfigurationError as e:
        print(f"✗ Configuration Error: {e}", file=sys.stderr)
        return 1
    except FileNotFoundError as e:
        print(f"✗ File Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"✗ Validation Error: {e}", file=sys.stderr)
        return 1


def cmd_template(args: argparse.Namespace) -> int:
    """
    Execute the template command.
    
    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments
        
    Returns
    -------
    int
        Exit code (0 for success, 1 for failure)
    """
    try:
        # Find the template file in the package
        template_path = Path(__file__).parent.parent / "pipeline_config_template.yaml"
        
        if not template_path.exists():
            print("✗ Template file not found in package", file=sys.stderr)
            return 1
        
        # Copy template to specified location
        output_path = Path(args.output)
        
        if output_path.exists() and not args.force:
            print(f"✗ Output file already exists: {output_path}", file=sys.stderr)
            print("Use --force to overwrite", file=sys.stderr)
            return 1
        
        # Read template and write to output
        with open(template_path, 'r') as src, open(output_path, 'w') as dst:
            dst.write(src.read())
        
        print(f"✓ Pipeline template created: {output_path}")
        print("Edit the configuration file to customize for your data and requirements")
        
        return 0
        
    except Exception as e:
        print(f"✗ Template Error: {e}", file=sys.stderr)
        return 1


def create_parser() -> argparse.ArgumentParser:
    """
    Create the command line argument parser.
    
    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        prog='ideal-genom',
        description='IDEAL-GENOM-QC: Comprehensive genomic quality control pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run pipeline with configuration
  ideal-genom run --config my_pipeline.yaml
  
  # Dry run to see what would be executed
  ideal-genom run --config my_pipeline.yaml --dry-run
  
  # Validate configuration without running
  ideal-genom validate --config my_pipeline.yaml
  
  # Generate configuration template
  ideal-genom template --output my_pipeline.yaml

For more information, see the documentation at:
https://ideal-genom-qc.readthedocs.io/
        """
    )
    
    parser.add_argument(
        '--version', 
        action='version', 
        version=f'%(prog)s {__version__}'
    )
    
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Set logging level (overrides config file)'
    )
    
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        metavar='COMMAND'
    )
    
    # Run command
    run_parser = subparsers.add_parser(
        'run',
        help='Execute the genomic QC pipeline'
    )
    run_parser.add_argument(
        '--config', '-c',
        required=True,
        help='Path to pipeline configuration YAML file'
    )
    run_parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be executed without running'
    )
    run_parser.set_defaults(func=cmd_run)
    
    # Validate command
    validate_parser = subparsers.add_parser(
        'validate',
        help='Validate pipeline configuration'
    )
    validate_parser.add_argument(
        '--config', '-c',
        required=True,
        help='Path to pipeline configuration YAML file'
    )
    validate_parser.set_defaults(func=cmd_validate)
    
    # Template command
    template_parser = subparsers.add_parser(
        'template',
        help='Generate pipeline configuration template'
    )
    template_parser.add_argument(
        '--output', '-o',
        default='pipeline_config.yaml',
        help='Output file for configuration template (default: pipeline_config.yaml)'
    )
    template_parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing output file'
    )
    template_parser.set_defaults(func=cmd_template)
    
    return parser


def main(argv: Optional[list] = None) -> int:
    """
    Main CLI entry point.
    
    Parameters
    ----------
    argv : list, optional
        Command line arguments (uses sys.argv if None)
        
    Returns
    -------
    int
        Exit code
    """
    parser = create_parser()
    args = parser.parse_args(argv)
    
    if not hasattr(args, 'func'):
        parser.print_help()
        return 1
    
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())