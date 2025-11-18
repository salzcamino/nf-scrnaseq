"""
Test suite for nf-scrnaseq pipeline modules

Tests validate:
- Python syntax in all Nextflow modules
- Nextflow variable substitution produces valid code
- Workflow structure and dependencies
- Configuration file validity
"""

import pytest
import re
from pathlib import Path
import textwrap


class TestModuleSyntax:
    """Test Python syntax in Nextflow modules"""

    @pytest.fixture
    def module_dir(self):
        """Get modules directory"""
        return Path(__file__).parent.parent / 'modules' / 'local'

    @pytest.fixture
    def all_modules(self, module_dir):
        """Get all .nf module files"""
        return list(module_dir.glob('*.nf'))

    def extract_python_code(self, nf_file):
        """Extract Python code from Nextflow shell block"""
        with open(nf_file, 'r') as f:
            content = f.read()

        # Extract shell block
        shell_match = re.search(r"shell:\s*'''(.*?)'''", content, re.DOTALL)
        if not shell_match:
            return None

        python_code = shell_match.group(1)

        # Remove Nextflow variable substitutions
        python_code = re.sub(r'!\{[^}]+\}', '"TEST_VALUE"', python_code)

        # Remove leading whitespace consistently
        python_code = textwrap.dedent(python_code)

        return python_code

    def test_all_modules_have_shell_blocks(self, all_modules):
        """Verify all modules have shell blocks"""
        modules_without_shell = []
        for module in all_modules:
            code = self.extract_python_code(module)
            if code is None:
                modules_without_shell.append(module.name)

        assert len(modules_without_shell) == 0, \
            f"Modules without shell blocks: {modules_without_shell}"

    def test_python_syntax_valid(self, all_modules):
        """Test that all module Python code is syntactically valid"""
        syntax_errors = []

        for module in all_modules:
            python_code = self.extract_python_code(module)
            if python_code:
                try:
                    compile(python_code, str(module), 'exec')
                except SyntaxError as e:
                    syntax_errors.append({
                        'module': module.name,
                        'line': e.lineno,
                        'message': e.msg
                    })

        assert len(syntax_errors) == 0, \
            f"Syntax errors found: {syntax_errors}"

    def test_nextflow_substitution_produces_valid_code(self, module_dir):
        """Test that Nextflow variable substitution produces valid Python"""
        import_data = module_dir / 'import_data.nf'

        if import_data.exists():
            code = self.extract_python_code(import_data)
            # Should be valid after substitution
            compile(code, str(import_data), 'exec')


class TestWorkflowStructure:
    """Test workflow structure and dependencies"""

    @pytest.fixture
    def main_workflow(self):
        """Load main workflow file"""
        main_nf = Path(__file__).parent.parent / 'main.nf'
        with open(main_nf, 'r') as f:
            return f.read()

    @pytest.fixture
    def module_dir(self):
        """Get modules directory"""
        return Path(__file__).parent.parent / 'modules' / 'local'

    def test_all_includes_have_modules(self, main_workflow, module_dir):
        """Verify all included processes have corresponding module files"""
        # Expected process to filename mapping
        process_to_file = {
            'IMPORT_DATA': 'import_data.nf',
            'QC_FILTER': 'qc_filter.nf',
            'DOUBLET_DECONTAM': 'doublet_decontam.nf',
            'NORMALIZE': 'normalize.nf',
            'HIGHLY_VARIABLE_GENES': 'highly_variable_genes.nf',
            'REDUCE_DIMS': 'reduce_dims.nf',
            'CLUSTERING': 'clustering.nf',
            'DIFF_EXPRESSION': 'diff_expression.nf',
            'CELL_TYPE_ANNOTATION': 'cell_type_annotation.nf',
            'GENE_SET_ENRICHMENT': 'gsea.nf',
            'BATCH_CORRECTION': 'batch_correction.nf',
            'CELL_CYCLE_SCORING': 'cell_cycle.nf',
            'TRAJECTORY_ANALYSIS': 'trajectory.nf',
            'CELL_COMMUNICATION': 'cell_communication.nf',
            'HTML_REPORT': 'html_report.nf',
            'SAMPLE_INTEGRATION': 'sample_integration.nf',
            'DATA_EXPORT': 'data_export.nf'
        }

        # Find all includes
        includes = re.findall(r'include\s+\{\s+(\w+)\s+\}', main_workflow)

        missing_modules = []
        for process_name in includes:
            filename = process_to_file.get(process_name)
            if filename:
                filepath = module_dir / filename
                if not filepath.exists():
                    missing_modules.append(filename)

        assert len(missing_modules) == 0, \
            f"Missing module files: {missing_modules}"

    def test_required_parameters_defined(self, main_workflow):
        """Verify required parameters are referenced in workflow"""
        required_params = [
            'input', 'outdir', 'min_genes', 'min_cells',
            'n_pcs', 'n_neighbors', 'leiden_resolution'
        ]

        missing_params = []
        for param in required_params:
            if f"params.{param}" not in main_workflow:
                missing_params.append(param)

        assert len(missing_params) == 0, \
            f"Required parameters not found in workflow: {missing_params}"

    def test_conditional_blocks_use_boolean_params(self, main_workflow):
        """Verify conditional blocks use valid boolean parameters"""
        conditionals = re.findall(r'if\s+\(params\.(\w+)\)', main_workflow)

        # These should be boolean flags
        expected_bool_params = [
            'run_doublet_detection', 'run_integration', 'run_cell_cycle',
            'run_batch_correction', 'run_diff_expression', 'run_annotation',
            'run_gsea', 'run_trajectory', 'run_communication', 'generate_report'
        ]

        # All conditional params should be in expected list
        unexpected = set(conditionals) - set(expected_bool_params) - {'help'}

        assert len(unexpected) == 0, \
            f"Unexpected conditional parameters: {unexpected}"


class TestConfiguration:
    """Test Nextflow configuration files"""

    @pytest.fixture
    def nextflow_config(self):
        """Load main Nextflow config"""
        config_file = Path(__file__).parent.parent / 'nextflow.config'
        with open(config_file, 'r') as f:
            return f.read()

    def test_config_has_required_profiles(self, nextflow_config):
        """Verify all required execution profiles are defined"""
        required_profiles = ['docker', 'conda', 'singularity', 'test']

        missing_profiles = []
        for profile in required_profiles:
            if f"{profile} {{" not in nextflow_config:
                missing_profiles.append(profile)

        assert len(missing_profiles) == 0, \
            f"Missing profiles: {missing_profiles}"

    def test_config_has_resource_labels(self, nextflow_config):
        """Verify resource labels are defined"""
        required_labels = ['process_low', 'process_medium', 'process_high']

        missing_labels = []
        for label in required_labels:
            if label not in nextflow_config:
                missing_labels.append(label)

        assert len(missing_labels) == 0, \
            f"Missing resource labels: {missing_labels}"

    def test_cloud_configs_exist(self):
        """Verify cloud execution config files exist"""
        conf_dir = Path(__file__).parent.parent / 'conf'
        required_configs = ['aws.config', 'gcp.config', 'slurm.config', 'pbs.config']

        missing_configs = []
        for config_file in required_configs:
            if not (conf_dir / config_file).exists():
                missing_configs.append(config_file)

        assert len(missing_configs) == 0, \
            f"Missing config files: {missing_configs}"


class TestTestData:
    """Test the test dataset validity"""

    @pytest.fixture
    def test_data_dir(self):
        """Get test data directory"""
        return Path(__file__).parent.parent / 'test_data' / '10x_sample'

    def test_required_files_exist(self, test_data_dir):
        """Verify all required test data files exist"""
        required_files = [
            'barcodes.tsv', 'genes.tsv', 'matrix.mtx',
            'cell_types_ground_truth.csv'
        ]

        missing_files = []
        for filename in required_files:
            if not (test_data_dir / filename).exists():
                missing_files.append(filename)

        assert len(missing_files) == 0, \
            f"Missing test data files: {missing_files}"

    def test_barcode_count(self, test_data_dir):
        """Verify test dataset has expected number of cells"""
        barcodes_file = test_data_dir / 'barcodes.tsv'
        if barcodes_file.exists():
            with open(barcodes_file) as f:
                barcodes = f.readlines()

            assert len(barcodes) == 200, \
                f"Expected 200 barcodes, found {len(barcodes)}"

    def test_gene_count(self, test_data_dir):
        """Verify test dataset has expected number of genes"""
        genes_file = test_data_dir / 'genes.tsv'
        if genes_file.exists():
            with open(genes_file) as f:
                genes = f.readlines()

            assert len(genes) == 105, \
                f"Expected 105 genes, found {len(genes)}"

    def test_key_markers_present(self, test_data_dir):
        """Verify key marker genes are in test dataset"""
        genes_file = test_data_dir / 'genes.tsv'
        if genes_file.exists():
            with open(genes_file) as f:
                gene_symbols = [line.strip().split('\t')[1] for line in f]

            key_markers = ['MT-CO1', 'CD3D', 'CD19', 'CD14', 'NKG7']
            missing_markers = [m for m in key_markers if m not in gene_symbols]

            assert len(missing_markers) == 0, \
                f"Missing key markers: {missing_markers}"


class TestDocumentation:
    """Test documentation completeness"""

    def test_readme_exists(self):
        """Verify README.md exists"""
        readme = Path(__file__).parent.parent / 'README.md'
        assert readme.exists(), "README.md not found"

    def test_readme_has_key_sections(self):
        """Verify README has essential sections"""
        readme = Path(__file__).parent.parent / 'README.md'
        with open(readme) as f:
            content = f.read()

        required_sections = [
            'Quick Start', 'Parameters', 'Output', 'Installation'
        ]

        missing_sections = []
        for section in required_sections:
            if section not in content:
                missing_sections.append(section)

        assert len(missing_sections) == 0, \
            f"Missing README sections: {missing_sections}"

    def test_dockerfile_exists(self):
        """Verify Dockerfile exists"""
        dockerfile = Path(__file__).parent.parent / 'Dockerfile'
        assert dockerfile.exists(), "Dockerfile not found"

    def test_environment_yml_exists(self):
        """Verify environment.yml exists"""
        env_file = Path(__file__).parent.parent / 'environment.yml'
        assert env_file.exists(), "environment.yml not found"
