"""Pipeline orchestration engine for ideal_genom."""

import os
import re
import importlib
import logging
from pathlib import Path

from typing import Dict, Any, List

logger = logging.getLogger(__name__)


class PipelineExecutor:
    """
    Orchestrates execution of sub-pipeline classes based on configuration.
    
    Attributes
    ----------
    config : dict
        Pipeline configuration dictionary
    steps : dict
        Dictionary of instantiated sub-pipeline objects
    base_output_dir : str
        Base output directory for all pipeline steps
    """
    
    def __init__(self, config: Dict[str, Any], dry_run: bool = False):
        """
        Initialize pipeline executor.
        
        Parameters
        ----------
        config : dict
            Pipeline configuration dictionary (from config.load_config)
        dry_run : bool
            If True, skip directory creation and actual execution
        """
        self.config = config
        self.steps = {}  # Store instantiated sub-pipeline objects
        self.base_output_dir = config['pipeline']['base_output_dir']
        self.pipeline_name = config['pipeline']['name']
        self.dry_run = dry_run
        
        # Create base output directory (skip if dry run)
        if not dry_run:
            os.makedirs(self.base_output_dir, exist_ok=True)
        
        # Setup logging
        self.logger = self._setup_logging()
        
    def execute(self) -> None:
        """Execute all pipeline steps sequentially."""
        pipeline_steps = self.config['pipeline']['steps']
        
        # Filter enabled steps and validate dependencies
        enabled_steps = self._filter_enabled_steps(pipeline_steps)
        self._validate_step_dependencies(enabled_steps)
        
        self.logger.info(f"Starting pipeline: {self.pipeline_name}")
        self.logger.info(f"Output directory: {self.base_output_dir}")
        self.logger.info(f"Total steps: {len(pipeline_steps)}")
        self.logger.info(f"Enabled steps: {len(enabled_steps)}")
        
        # First pass: Instantiate ALL classes (enabled and disabled) for reference resolution
        self.logger.info(f"\nInstantiating all pipeline classes for reference resolution...")
        for step_config in pipeline_steps:
            step_name = step_config['name']
            is_enabled = step_config.get('enabled', True)
            
            if not is_enabled:
                self.logger.info(f"Instantiating disabled step for references: {step_name}")
                try:
                    self._instantiate_step(step_config)
                except Exception as e:
                    self.logger.warning(f"⚠️  Failed to instantiate disabled step '{step_name}': {e}")
                    # Continue anyway - this step won't be available for references
        
        # Second pass: Execute only enabled steps
        for i, step_config in enumerate(enabled_steps, 1):
            step_name = step_config['name']
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"Step {i}/{len(enabled_steps)}: {step_name}")
            self.logger.info(f"{'='*60}")
            
            try:
                # If already instantiated (from first pass), just execute
                if step_name in self.steps:
                    self._execute_existing_step(step_config)
                else:
                    # Instantiate and execute
                    self._execute_step(step_config)
                self.logger.info(f"✓ Completed step: {step_name}")
            except Exception as e:
                self.logger.error(f"✗ Failed step: {step_name}")
                self.logger.error(f"Error: {str(e)}", exc_info=True)
                raise RuntimeError(f"Pipeline failed at step '{step_name}': {str(e)}")
        
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"Pipeline completed successfully: {self.pipeline_name}")
        self.logger.info(f"{'='*60}")
    
    def _execute_step(self, step_config: Dict[str, Any]) -> None:
        """
        Execute a single pipeline step.
        
        Parameters
        ----------
        step_config : dict
            Configuration for the step to execute
        """
        step_name = step_config['name']
        
        # Import the class dynamically
        module_path = step_config['module']
        class_name = step_config['class']
        
        self.logger.info(f"Loading {class_name} from {module_path}")
        
        try:
            module = importlib.import_module(module_path)
            pipeline_class = getattr(module, class_name)
        except (ImportError, AttributeError) as e:
            raise ImportError(
                f"Failed to import {class_name} from {module_path}: {str(e)}"
            )
        
        # Resolve parameters (handle ${...} references)
        init_params = self._resolve_params(step_config['init_params'])
        execute_params = self._resolve_params(step_config.get('execute_params', {}))
        
        # Convert string paths to Path objects for parameters ending with '_path'
        init_params = self._convert_paths_to_path_objects(init_params)
        execute_params = self._convert_paths_to_path_objects(execute_params)
        
        self.logger.info(f"Initializing {class_name}")
        self.logger.debug(f"Init params: {init_params}")
        self.logger.debug(f"Execute params: {execute_params}")
        
        # Instantiate the sub-pipeline class
        pipeline_instance = pipeline_class(**init_params)
        
        # Determine the execute method name
        # Convention: execute_<step_name>_pipeline
        execute_method_name = f'execute_{step_name}_pipeline'
        
        if not hasattr(pipeline_instance, execute_method_name):
            raise AttributeError(
                f"{class_name} does not have method '{execute_method_name}'"
            )
        
        execute_method = getattr(pipeline_instance, execute_method_name)
        
        self.logger.info(f"Executing {execute_method_name}")
        execute_method(execute_params)
        
        # Store the instance for reference by subsequent steps
        self.steps[step_name] = pipeline_instance
        
        # Generate reports if configured
        self._generate_reports(step_name, pipeline_instance)
        
        # Perform cleanup if configured
        self._perform_cleanup(step_name, pipeline_instance)
        
        self.logger.info(f"Step output stored as: steps.{step_name}")
    
    def _instantiate_step(self, step_config: Dict[str, Any]) -> None:
        """
        Instantiate a pipeline step class without executing it.
        Used for disabled steps to enable reference resolution.
        
        Parameters
        ----------
        step_config : dict
            Configuration for the step to instantiate
        """
        step_name = step_config['name']
        
        # Import the class dynamically
        module_path = step_config['module']
        class_name = step_config['class']
        
        self.logger.debug(f"Loading {class_name} from {module_path} (instantiate only)")
        
        try:
            module = importlib.import_module(module_path)
            pipeline_class = getattr(module, class_name)
        except (ImportError, AttributeError) as e:
            raise ImportError(
                f"Failed to import {class_name} from {module_path}: {str(e)}"
            )
        
        # Resolve parameters (handle ${...} references)
        init_params = self._resolve_params(step_config['init_params'])
        
        # Convert string paths to Path objects for parameters ending with '_path'
        init_params = self._convert_paths_to_path_objects(init_params)
        
        self.logger.debug(f"Instantiating {class_name} (disabled step)")
        self.logger.debug(f"Init params: {init_params}")
        
        # Instantiate the sub-pipeline class
        pipeline_instance = pipeline_class(**init_params)
        
        # Store the instance for reference by subsequent steps
        self.steps[step_name] = pipeline_instance
        
        self.logger.debug(f"Disabled step instantiated and stored as: steps.{step_name}")
    
    def _execute_existing_step(self, step_config: Dict[str, Any]) -> None:
        """
        Execute a step that has already been instantiated.
        
        Parameters
        ----------
        step_config : dict
            Configuration for the step to execute
        """
        step_name = step_config['name']
        
        # Get the already instantiated pipeline object
        pipeline_instance = self.steps[step_name]
        
        # Resolve execute parameters
        execute_params = self._resolve_params(step_config.get('execute_params', {}))
        execute_params = self._convert_paths_to_path_objects(execute_params)
        
        # Determine the execute method name
        execute_method_name = f'execute_{step_name}_pipeline'
        
        if not hasattr(pipeline_instance, execute_method_name):
            raise AttributeError(
                f"{pipeline_instance.__class__.__name__} does not have method '{execute_method_name}'"
            )
        
        execute_method = getattr(pipeline_instance, execute_method_name)
        
        self.logger.info(f"Executing {execute_method_name}")
        self.logger.debug(f"Execute params: {execute_params}")
        execute_method(execute_params)

        # Generate reports if configured
        self._generate_reports(step_name, pipeline_instance)
        
        # Perform cleanup if configured
        self._perform_cleanup(step_name, pipeline_instance)

        return
    
    def _perform_cleanup(self, step_name: str, pipeline_instance: Any) -> None:
        """
        Perform cleanup of intermediate files for Sample QC and Variant QC steps.
        
        Parameters
        ----------
        step_name : str
            Name of the pipeline step
        pipeline_instance : Any
            Instance of the pipeline class that was executed
        """
        # Check if cleanup is disabled globally
        keep_intermediate = self.config.get('settings', {}).get('files', {}).get('keep_intermediate', True)
        
        self.logger.info(f"🔍 Cleanup check for {step_name}:")
        self.logger.info(f"   - keep_intermediate setting: {keep_intermediate}")
        self.logger.info(f"   - Config path: settings.files.keep_intermediate")
        
        if keep_intermediate:
            self.logger.info(f"⏭️  Skipping cleanup for {step_name} - keep_intermediate is True")
            return
        
        # Only handle sample_qc and variant_qc steps
        if step_name not in ['sample_qc', 'variant_qc']:
            self.logger.info(f"⏭️  No cleanup configured for step: {step_name}")
            return
        
        self.logger.info(f"🧹 Initiating cleanup for {step_name}...")
        
        try:
            if step_name == 'sample_qc':
                self._cleanup_sample_qc(pipeline_instance)
            elif step_name == 'variant_qc':
                self._cleanup_variant_qc(pipeline_instance)
                
            self.logger.info(f"✅ Cleanup completed for {step_name}")
            
        except Exception as e:
            self.logger.warning(f"⚠️  Cleanup failed for {step_name}: {e}")
            # Don't fail the pipeline for cleanup issues
    
    def _cleanup_sample_qc(self, pipeline_instance: Any) -> None:
        """Cleanup intermediate files for Sample QC step."""
        from ideal_genom.qc.sample_qc import SampleQCCleanUp
        
        # Get required paths from pipeline instance
        if not hasattr(pipeline_instance, 'output_path') or not hasattr(pipeline_instance, 'input_path'):
            self.logger.warning("❌ SampleQC instance missing required paths for cleanup")
            self.logger.warning(f"   Available attributes: {[attr for attr in dir(pipeline_instance) if not attr.startswith('_')]}")
            return
        
        output_path = pipeline_instance.results_dir
        input_path = pipeline_instance.input_path
        
        self.logger.info(f"🧹 Running Sample QC cleanup...")
        self.logger.info(f"   - Output path: {output_path}")
        self.logger.info(f"   - Input path: {input_path}")
        
        cleanup = SampleQCCleanUp(
            output_path=output_path,
            input_path=input_path
        )
        cleanup.clean_all()
    
    def _cleanup_variant_qc(self, pipeline_instance: Any) -> None:
        """Cleanup intermediate files for Variant QC step.""" 
        from ideal_genom.qc.variant_qc import VariantQCCleanUp
        
        # Get required paths from pipeline instance
        if not hasattr(pipeline_instance, 'output_path'):
            self.logger.warning("❌ VariantQC instance missing output_path for cleanup")
            self.logger.warning(f"   Available attributes: {[attr for attr in dir(pipeline_instance) if not attr.startswith('_')]}")
            return
        
        output_path = pipeline_instance.results_dir
        
        self.logger.info(f"🧹 Running Variant QC cleanup...")
        self.logger.info(f"   - Output path: {output_path}")
        
        cleanup = VariantQCCleanUp(output_path=output_path)
        cleanup.clean_all()
    
    def _generate_reports(self, step_name: str, pipeline_instance: Any) -> None:
        """
        Generate visualization reports for completed QC steps.
        
        Parameters
        ----------
        step_name : str
            Name of the pipeline step
        pipeline_instance : Any
            Instance of the pipeline class that was executed
        """
        # Check if report generation is enabled globally
        generate_reports = self.config.get('settings', {}).get('reports', {}).get('generate_reports', True)
        
        self.logger.info(f"📊 Report check for {step_name}:")
        self.logger.info(f"   - generate_reports setting: {generate_reports}")
        
        if not generate_reports:
            self.logger.info(f"⏭️  Skipping report generation for {step_name} - generate_reports is False")
            return
        
        # Only handle sample_qc, ancestry_qc, and variant_qc steps
        if step_name not in ['sample_qc', 'ancestry_qc', 'variant_qc']:
            self.logger.info(f"⏭️  No report generation configured for step: {step_name}")
            return
        
        self.logger.info(f"📊 Generating reports for {step_name}...")
        
        try:
            if step_name == 'sample_qc':
                self._generate_sample_qc_report(pipeline_instance)
            elif step_name == 'ancestry_qc':
                self._generate_ancestry_qc_report(pipeline_instance)
            elif step_name == 'variant_qc':
                self._generate_variant_qc_report(pipeline_instance)
                
            self.logger.info(f"✅ Report generation completed for {step_name}")
            
        except Exception as e:
            self.logger.warning(f"⚠️  Report generation failed for {step_name}: {e}")
            # Don't fail the pipeline for report generation issues
    
    def _generate_sample_qc_report(self, pipeline_instance: Any) -> None:
        """Generate visualization reports for Sample QC step."""
        from ideal_genom.qc.sample_qc import SampleQCReport
        
        # Check for required attributes
        required_attrs = ['results_dir']
        if not all(hasattr(pipeline_instance, attr) for attr in required_attrs):
            self.logger.warning("❌ SampleQC instance missing required attributes for report generation")
            self.logger.warning(f"   Available attributes: {[attr for attr in dir(pipeline_instance) if not attr.startswith('_')]}")
            return
        
        plots_dir = pipeline_instance.results_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        self.logger.info(f"📊 Generating Sample QC reports...")
        self.logger.info(f"   - Plots directory: {plots_dir}")
        
        # Create report instance
        report = SampleQCReport(output_path=plots_dir)
        
        # Get file paths from sample_qc instance
        call_rate_smiss = getattr(pipeline_instance, 'call_rate_miss', None)
        sexcheck_miss = getattr(pipeline_instance, 'sexcheck_miss', None)
        xchr_miss = getattr(pipeline_instance, 'xchr_miss', None)
        maf_greater_het = getattr(pipeline_instance, 'maf_greater_het', None)
        maf_less_het = getattr(pipeline_instance, 'maf_less_het', None)
        maf_greater_smiss = getattr(pipeline_instance, 'maf_greater_smiss', None)
        maf_less_smiss = getattr(pipeline_instance, 'maf_less_smiss', None)
        genome = getattr(pipeline_instance, 'genome', None)
        
        # Check if required files exist
        required_files = {
            'call_rate_smiss': call_rate_smiss,
            'sexcheck_miss': sexcheck_miss, 
            'xchr_miss': xchr_miss,
            'maf_greater_het': maf_greater_het,
            'maf_less_het': maf_less_het,
            'maf_greater_smiss': maf_greater_smiss,
            'maf_less_smiss': maf_less_smiss
        }
        
        missing_files = []
        for name, file_path in required_files.items():
            if file_path is None:
                missing_files.append(f"{name} (not set)")
            elif not file_path.exists():
                missing_files.append(f"{name} ({file_path})")
        
        if missing_files:
            self.logger.warning(f"❌ Sample QC report: Missing required files: {missing_files}")
            return
        
        # Get step configuration for parameters
        step_config = self.config.get('sample_qc', {})
        call_rate_thres = step_config.get('call_rate_thres', 0.02)
        std_deviation_het = step_config.get('std_deviation_het', 3)
        f_coeff_threshold = step_config.get('sex_check', [0.2, 0.8])
        maf_het = step_config.get('maf_het', 0.01)
        ibd_threshold = step_config.get('ibd_threshold', None)
        generate_ibd_report = genome is not None and genome.exists() if genome else False
        
        try:
            # Type assertion since we've already validated these files exist
            assert call_rate_smiss is not None and call_rate_smiss.exists()
            assert sexcheck_miss is not None and sexcheck_miss.exists()
            assert xchr_miss is not None and xchr_miss.exists()
            assert maf_greater_het is not None and maf_greater_het.exists()
            assert maf_less_het is not None and maf_less_het.exists()
            assert maf_greater_smiss is not None and maf_greater_smiss.exists()
            assert maf_less_smiss is not None and maf_less_smiss.exists()
            
            # Generate the report
            report.report_sample_qc(
                call_rate_smiss=call_rate_smiss,
                sexcheck_miss=sexcheck_miss,
                xchr_miss=xchr_miss,
                maf_greater_het=maf_greater_het,
                maf_less_het=maf_less_het,
                maf_greater_smiss=maf_greater_smiss,
                maf_less_smiss=maf_less_smiss,
                genome=genome,
                generate_ibd_report=generate_ibd_report,
                f_coeff_thresholds=f_coeff_threshold,
                call_rate_thres=call_rate_thres,
                std_deviation_het=std_deviation_het,
                maf_het=maf_het,
                ibd_threshold=ibd_threshold
            )
            
            self.logger.info(f"✅ Sample QC reports generated successfully in {plots_dir}")
            
        except Exception as e:
            self.logger.error(f"❌ Failed to generate Sample QC report: {e}")
            self.logger.debug(f"Sample QC report error details", exc_info=True)
    
    def _generate_ancestry_qc_report(self, pipeline_instance: Any) -> None:
        """Generate visualization reports for Ancestry QC step."""
        from ideal_genom.qc.ancestry_qc import AncestryQCReport
        
        # Check for required attributes  
        required_attrs = ['results_dir', 'eigenvectors', 'eigenvalues', 'ancestry_fails', 'population_tags']
        missing_attrs = [attr for attr in required_attrs if not hasattr(pipeline_instance, attr)]
        
        if missing_attrs:
            self.logger.warning(f"❌ AncestryQC instance missing required attributes: {missing_attrs}")
            return
        
        # Ensure all file paths exist
        file_attrs = ['eigenvectors', 'eigenvalues', 'ancestry_fails', 'population_tags']
        missing_files = []
        for attr in file_attrs:
            file_path = getattr(pipeline_instance, attr)
            if file_path is None or not file_path.exists():
                missing_files.append(attr)
        
        if missing_files:
            self.logger.warning(f"❌ AncestryQC missing or invalid file paths: {missing_files}")
            return
        
        plots_dir = pipeline_instance.results_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        self.logger.info(f"📊 Generating Ancestry QC reports...")
        self.logger.info(f"   - Plots directory: {plots_dir}")
        
        # Create report instance
        report = AncestryQCReport(
            output_path=plots_dir,
            einvectors=pipeline_instance.eigenvectors,
            eigenvalues=pipeline_instance.eigenvalues,
            ancestry_fails=pipeline_instance.ancestry_fails,
            population_tags=pipeline_instance.population_tags
        )
        
        # Generate comprehensive ancestry QC report
        # Get reference population from config or use default
        ref_pop = self.config.get('pipeline', {}).get('steps', [])
        ancestry_step = next((step for step in ref_pop if step.get('name') == 'ancestry_qc'), {})
        reference_pop = ancestry_step.get('execute_params', {}).get('reference_pop', 'EUR')
        
        # Get plot format from config
        plot_format = self.config.get('settings', {}).get('reports', {}).get('plot_format', 'png')
        
        report.report_ancestry_qc(
            reference_pop=reference_pop,
            aspect_ratio='equal',
            format=plot_format
        )
    
    def _generate_variant_qc_report(self, pipeline_instance: Any) -> None:
        """Generate visualization reports for Variant QC step.""" 
        from ideal_genom.qc.variant_qc import VariantQCReport
        
        # Check for required attributes
        required_attrs = ['results_dir']
        if not all(hasattr(pipeline_instance, attr) for attr in required_attrs):
            self.logger.warning("❌ VariantQC instance missing required attributes for report generation")
            return
        
        plots_dir = pipeline_instance.results_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        self.logger.info(f"📊 Generating Variant QC reports...")
        self.logger.info(f"   - Plots directory: {plots_dir}")
        
        # Create report instance  
        report = VariantQCReport(output_path=plots_dir)
        
        # Get file paths from variant_qc instance  
        males_missing_data = getattr(pipeline_instance, 'males_missing_data', None)
        females_missing_data = getattr(pipeline_instance, 'females_missing_data', None)
        
        # Check if required files exist
        required_files = {
            'males_missing_data': males_missing_data,
            'females_missing_data': females_missing_data
        }
        
        missing_files = []
        for name, file_path in required_files.items():
            if file_path is None:
                missing_files.append(f"{name} (not set)")
            elif not file_path.exists():
                missing_files.append(f"{name} ({file_path})")
        
        if missing_files:
            self.logger.warning(f"❌ Variant QC report: Missing required files: {missing_files}")
            return
        
        # Get step configuration for parameters
        step_config = self.config.get('variant_qc', {})
        y_axis_cap = step_config.get('y_axis_cap', 100)
        missing_data_threshold = step_config.get('missing_data_threshold', 0.1)
        
        try:
            # Type assertion since we've already validated these files exist
            assert males_missing_data is not None and males_missing_data.exists()
            assert females_missing_data is not None and females_missing_data.exists()
            
            # Generate the report
            report.report_variant_qc(
                missing_data_rate_male=males_missing_data,
                missing_data_rate_female=females_missing_data,
                y_axis_cap=y_axis_cap,
                missing_data_threshold=missing_data_threshold
            )
            
            self.logger.info(f"✅ Variant QC reports generated successfully in {plots_dir}")
            
        except Exception as e:
            self.logger.error(f"❌ Failed to generate Variant QC report: {e}")
            self.logger.debug(f"Variant QC report error details", exc_info=True)
        self.logger.info("📊 Variant QC report generation requires specific file detection")
        self.logger.info("📊 Report framework ready - implement file detection for full automation")
    
    
    def _resolve_params(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Resolve parameter values, including references to previous steps.
        
        Parameters
        ----------
        params : dict
            Parameters dictionary potentially containing references
            
        Returns
        -------
        dict
            Parameters with all references resolved
        """
        resolved = {}
        
        for key, value in params.items():
            if isinstance(value, str):
                resolved[key] = self._resolve_string_value(value)
            elif isinstance(value, dict):
                resolved[key] = self._resolve_params(value)
            elif isinstance(value, list):
                resolved[key] = [
                    self._resolve_string_value(v) if isinstance(v, str) else v
                    for v in value
                ]
            else:
                resolved[key] = value
        
        return resolved
    
    def _resolve_string_value(self, value: str) -> Any:
        """
        Resolve a string value that may contain references.
        
        Supports:
        - ${base_output_dir} - pipeline base output directory
        - ${steps.step_name.attribute} - attribute from previous step
        
        Parameters
        ----------
        value : str
            String value potentially containing references
            
        Returns
        -------
        Any
            Resolved value
        """
        # Find all ${...} patterns
        pattern = r'\$\{([^}]+)\}'
        matches = re.findall(pattern, value)
        
        if not matches:
            return value
        
        resolved_value = value
        
        for match in matches:
            ref_value = self._resolve_reference(match)
            resolved_value = resolved_value.replace(f'${{{match}}}', str(ref_value))
        
        return resolved_value
    
    def _resolve_reference(self, reference: str) -> Any:
        """
        Resolve a reference string to its actual value.
        
        Parameters
        ----------
        reference : str
            Reference string (without ${ })
            
        Returns
        -------
        Any
            Resolved value
            
        Raises
        ------
        ValueError
            If reference is invalid or step not found
        """
        # Handle base_output_dir
        if reference == 'base_output_dir':
            return self.base_output_dir
        
        # Handle steps.step_name.attribute
        if reference.startswith('steps.'):
            parts = reference.split('.')
            
            if len(parts) < 3:
                raise ValueError(
                    f"Invalid step reference: {reference}. "
                    f"Expected format: steps.step_name.attribute"
                )
            
            step_name = parts[1]
            attr_path = parts[2:]
            
            if step_name not in self.steps:
                # Check if step exists in config but hasn't been processed yet
                all_step_names = [s['name'] for s in self.config['pipeline']['steps']]
                if step_name in all_step_names:
                    raise ValueError(
                        f"Step '{step_name}' referenced but not yet instantiated. "
                        f"This may indicate the step is disabled or comes later in the pipeline. "
                        f"Available steps: {list(self.steps.keys())}"
                    )
                else:
                    raise ValueError(
                        f"Step '{step_name}' not found in pipeline configuration. "
                        f"The current step references '${{{reference}}}' but '{step_name}' "
                        f"is not defined in your YAML configuration file. "
                        f"Available steps in config: {all_step_names}. "
                        f"Please either add the '{step_name}' step to your configuration "
                        f"or update the reference to use an available step."
                    )
            
            # Navigate through nested attributes
            obj = self.steps[step_name]
            for attr in attr_path:
                if not hasattr(obj, attr):
                    raise ValueError(
                        f"Step '{step_name}' does not have attribute '{attr}'"
                    )
                obj = getattr(obj, attr)
            
            return obj
        
        raise ValueError(f"Unknown reference: {reference}")
    
    def _convert_paths_to_path_objects(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Convert string parameters ending with '_path' or '_file' to Path objects.
        
        Parameters
        ----------
        params : dict
            Parameters dictionary potentially containing path strings
            
        Returns
        -------
        dict
            Parameters with path strings converted to Path objects
        """
        converted = {}
        
        for key, value in params.items():
            if (key.endswith('_path') or key.endswith('_file')) and isinstance(value, str):
                converted[key] = Path(value)
            else:
                converted[key] = value
        
        return converted
    
    def _filter_enabled_steps(self, pipeline_steps: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Filter pipeline steps to include only enabled ones.
        
        Parameters
        ----------
        pipeline_steps : list
            List of all pipeline step configurations
            
        Returns
        -------
        list
            List of enabled pipeline steps
        """
        enabled_steps = []
        
        for step in pipeline_steps:
            # Default to enabled if not specified
            if step.get('enabled', True):
                enabled_steps.append(step)
            else:
                self.logger.info(f"Skipping disabled step: {step['name']}")
        
        return enabled_steps
    
    def _validate_step_dependencies(self, enabled_steps: List[Dict[str, Any]]) -> None:
        """
        Validate step dependencies and issue warnings for potential issues.
        
        Parameters
        ----------
        enabled_steps : list
            List of enabled pipeline steps
        """
        enabled_step_names = {step['name'] for step in enabled_steps}
        all_step_names = {step['name'] for step in self.config['pipeline']['steps']}
        
        # Define step dependencies for genomic QC workflow
        dependencies = {
            'variant_qc': ['sample_qc'],
            'ancestry_qc': ['sample_qc'],  
            'population_analysis': []  # Can run standalone or after any QC step
        }
        
        # Check for missing dependencies and issue warnings
        for step in enabled_steps:
            step_name = step['name']
            
            if step_name in dependencies:
                missing_deps = []
                unavailable_deps = []
                
                for dep in dependencies[step_name]:
                    if dep not in enabled_step_names:
                        missing_deps.append(dep)
                        if dep not in all_step_names:
                            unavailable_deps.append(dep)
                
                if missing_deps:
                    if unavailable_deps:
                        self.logger.warning(
                            f"⚠️  Step '{step_name}' is enabled but dependency "
                            f"step(s) {unavailable_deps} are not defined in configuration."
                        )
                    else:
                        self.logger.info(
                            f"ℹ️  Step '{step_name}' is enabled but dependency "
                            f"step(s) {missing_deps} are disabled. "
                            f"Disabled steps will be instantiated for reference resolution."
                        )
        
        # Check for proper step ordering
        self._validate_step_ordering(enabled_steps)
    
    def _validate_step_ordering(self, enabled_steps: List[Dict[str, Any]]) -> None:
        """
        Validate that steps are in proper execution order.
        
        Parameters
        ----------
        enabled_steps : list
            List of enabled pipeline steps
        """
        # Define preferred step order
        preferred_order = ['sample_qc', 'ancestry_qc', 'variant_qc', 'population_analysis']
        
        # Get positions of enabled steps in preferred order
        step_positions = {}
        for step in enabled_steps:
            step_name = step['name']
            if step_name in preferred_order:
                step_positions[step_name] = preferred_order.index(step_name)
        
        # Check if steps are in correct order
        if len(step_positions) > 1:
            current_order = [step['name'] for step in enabled_steps if step['name'] in preferred_order]
            sorted_order = sorted(current_order, key=lambda x: preferred_order.index(x))
            
            if current_order != sorted_order:
                self.logger.warning(
                    f"⚠️  Steps may be out of optimal order. "
                    f"Current: {current_order}, Recommended: {sorted_order}. "
                    f"This may cause issues with step dependencies."
                )

    def _setup_logging(self) -> logging.Logger:
        """
        Setup logging for the pipeline.
        
        Returns
        -------
        logging.Logger
            Configured logger instance
        """
        # Create logger
        pipeline_logger = logging.getLogger(f'ideal_genom.pipeline.{self.pipeline_name}')
        pipeline_logger.setLevel(logging.DEBUG)
        
        # Avoid duplicate handlers
        if pipeline_logger.handlers:
            return pipeline_logger
        
        # Console handler (always available)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # File handler (only if not dry run)
        if not self.dry_run:
            log_dir = os.path.join(self.base_output_dir, 'pipeline_logs')
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(log_dir, f'{self.pipeline_name}.log')
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            pipeline_logger.addHandler(file_handler)

        console_handler.setFormatter(formatter)
        pipeline_logger.addHandler(console_handler)
        
        return pipeline_logger
    
    def get_step_output(self, step_name: str, attribute: str = 'output_path') -> Any:
        """
        Get output from a completed step.
        
        Parameters
        ----------
        step_name : str
            Name of the step
        attribute : str, optional
            Attribute to retrieve (default: 'output_path')
            
        Returns
        -------
        Any
            Value of the requested attribute
            
        Raises
        ------
        ValueError
            If step not found or attribute doesn't exist
        """
        if step_name not in self.steps:
            raise ValueError(f"Step '{step_name}' not found")
        
        step_instance = self.steps[step_name]
        
        if not hasattr(step_instance, attribute):
            raise ValueError(
                f"Step '{step_name}' does not have attribute '{attribute}'"
            )
        
        return getattr(step_instance, attribute)
    
    def get_pipeline_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the pipeline configuration and status.
        
        Returns
        -------
        dict
            Pipeline summary including enabled steps, dependencies, and configuration
        """
        pipeline_steps = self.config['pipeline']['steps']
        enabled_steps = self._filter_enabled_steps(pipeline_steps)
        
        summary = {
            'pipeline_name': self.pipeline_name,
            'base_output_dir': self.base_output_dir,
            'total_steps': len(pipeline_steps),
            'enabled_steps': len(enabled_steps),
            'steps': []
        }
        
        for step in enabled_steps:
            step_info = {
                'name': step['name'],
                'module': step['module'],
                'class': step['class'],
                'enabled': step.get('enabled', True)
            }
            summary['steps'].append(step_info)
        
        return summary
