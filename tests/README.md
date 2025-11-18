# nf-scrnaseq Test Suite

Automated tests for the nf-scrnaseq pipeline.

## Test Coverage

### 1. Module Syntax Tests (`test_modules.py::TestModuleSyntax`)
- ✅ Verify all modules have shell blocks
- ✅ Validate Python syntax in all 17 modules
- ✅ Test Nextflow variable substitution

### 2. Workflow Structure Tests (`test_modules.py::TestWorkflowStructure`)
- ✅ Verify all included processes have module files
- ✅ Check required parameters are defined
- ✅ Validate conditional blocks use boolean parameters

### 3. Configuration Tests (`test_modules.py::TestConfiguration`)
- ✅ Verify required execution profiles (docker, conda, singularity, test)
- ✅ Check resource labels are defined
- ✅ Validate cloud config files exist

### 4. Test Data Tests (`test_modules.py::TestTestData`)
- ✅ Verify required test data files exist
- ✅ Validate cell count (200 cells)
- ✅ Validate gene count (105 genes)
- ✅ Check key marker genes are present

### 5. Documentation Tests (`test_modules.py::TestDocumentation`)
- ✅ Verify README.md exists
- ✅ Check essential README sections
- ✅ Validate Dockerfile exists
- ✅ Validate environment.yml exists

## Running Tests

### Install Dependencies

```bash
pip install -r requirements-dev.txt
```

### Run All Tests

```bash
pytest
```

### Run Specific Test Classes

```bash
# Test module syntax only
pytest tests/test_modules.py::TestModuleSyntax -v

# Test workflow structure
pytest tests/test_modules.py::TestWorkflowStructure -v

# Test configuration
pytest tests/test_modules.py::TestConfiguration -v
```

### Run with Coverage

```bash
pytest --cov=modules --cov-report=html --cov-report=term
```

### Run Fast Tests Only

```bash
pytest -m "not slow"
```

## Test Markers

- `@pytest.mark.unit` - Unit tests
- `@pytest.mark.integration` - Integration tests
- `@pytest.mark.slow` - Slow-running tests

## Continuous Integration

Tests run automatically on:
- Every push to `main` or `dev` branches
- Every pull request

See `.github/workflows/ci.yml` for CI configuration.

## Adding New Tests

1. Add test functions to appropriate test class
2. Use descriptive test names: `test_<what_is_being_tested>`
3. Add docstrings explaining what the test validates
4. Run locally before committing: `pytest -v`

### Example

```python
def test_new_feature(self):
    """Verify new feature works correctly"""
    result = some_function()
    assert result == expected_value, "Feature should return expected value"
```

## Test Reports

After running tests, view coverage report:

```bash
# Generate HTML coverage report
pytest --cov=modules --cov-report=html

# Open in browser
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

## Troubleshooting

### Import Errors

If you get import errors, make sure you're running from the project root:

```bash
cd /path/to/nf-scrnaseq
pytest
```

### Missing Dependencies

Install all development dependencies:

```bash
pip install -r requirements-dev.txt
```

### Tests Failing

1. Check if you're on the correct branch
2. Ensure all module files are present
3. Run with verbose output: `pytest -vv`
4. Check individual test: `pytest tests/test_modules.py::TestClass::test_function -v`
