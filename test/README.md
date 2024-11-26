# Test Directory Structure

This directory is meant for the following purposes:
- Integration testing of individual functions/class objects within Troute
- Hosting example data to run full-length test runs

## Quick Start

```bash
# Run all tests
pytest test/

# Run tests from a specific category
pytest test/troute-config/
pytest test/troute-network/
```

## Test Categories

### Config Tests (`test/troute-config/`)

**Usage:**
- Validating existing Troute example configs
- Testing the Pydantic Config Validations

### Network Tests (`test/troute-network/`)

**Usage:**
- Using YAML files from `NHD` and `HYFeatures` in order to test network creation

### NWM Tests (`test/troute-nwm/`)

**Usage:**
- end to end tests of individual functions within the following functions:
  - `main_v03()`
  - `main_v04()`

### Routing Tests (`test/troute-config/`)

**Usage:**
- Integration tests for NHD and HYFeature util functions

### End to End Examples:
##### `test/LowerColorado_TX/`
##### `test/LowerColorado_TX_v4/`
##### `test/unit_test_hyfeature/`

**Usage:**
- Full datasets and input files for verifying troute is up and running

## Test Fixtures

This section describes the global test fixtures available in `conftest.py` for testing the Troute.

## Data Fixtures

### `hyfeature_qlat_data`
Returns a list of dictionaries containing lateral flow data for HYFeatures testing:
```python
{
    "qlat_files": ["202304020000.CHRTOUT_DOMAIN1.csv", ...],  # 24 hourly files
    "nts": 288,
    "final_timestamp": datetime(2023, 4, 2, 23, 0)
}
```
Used for testing lateral flow input handling in HYFeatures format.

### `nhd_qlat_data`
Returns a dictionary containing NHD lateral flow test data:
```python
{
    "qlat_files": ["202108231400.CHRTOUT_DOMAIN1", ...],  # 24 hourly files
    "nts": 288,
    "final_timestamp": datetime(2021, 8, 24, 13, 0)
}
```
Used for testing lateral flow input handling in NHD format.

### `nhd_validation_files`
Returns a dictionary containing a list of validation files used for testing:
```python
{
    "validation_files": ["202108231400.CHRTOUT_DOMAIN1", ...]  # 24 hourly files
}
```

## Network Configuration Fixtures

### `nhd_built_test_network`
Provides a fully built NHD test network with preprocessed data including:
- Network connections
- Parameter dataframes
- Waterbody connections and types
- Network topology
- Gage crosswalks
- Diffusive network data
- Topobathy data

### `warmstart_nhd_test`
Returns warm start configuration for NHD testing including:
- Waterbodies dataframe
- Initial conditions (q0)
- Time settings (t0)
- Last observation dataframe
- Data assimilation parameters

## Path and Configuration Fixtures

### `nhd_test_files`
Returns a tuple of paths needed for NHD testing:
```python
(cwd, path, config)  # Current working dir, test path, and config file path
```
Default test directory: `test/LowerColorado_TX/`

### `hyfeatures_test_data`
Returns a tuple of paths needed for HYFeatures testing:
```python
(cwd, path, config)  # Current working dir, test path, and config file path
```
Default test directory: `test/LowerColorado_TX_v4/`

### `validated_config`
Provides a validation function for configuration files that:
- Validates configs against schema
- Fixes relative paths
- Supports strict and non-strict validation modes
- Returns a validated Config object

## Network Object Fixtures

### `nhd_test_network`
Returns a complete NHD test network configuration dictionary containing:
- Logging parameters
- Preprocessing settings
- Network topology settings
- Waterbody configuration
- Computation parameters
- Forcing configuration
- Restart settings
- Hybrid routing parameters
- Output configuration
- Parity check settings
- Data assimilation parameters

### `hyfeatures_test_network`
Similar to `nhd_test_network` but configured for HYFeatures testing.

### `hyfeatures_network_object`
Returns an initialized `HYFeaturesNetwork` object with all parameters configured for testing.

## Usage Examples

```python
def test_qlat_processing(hyfeature_qlat_data):
    # Test lateral flow data processing
    assert len(hyfeature_qlat_data[0]['qlat_files']) == 24
    assert hyfeature_qlat_data[0]['nts'] == 288

def test_network_build(nhd_built_test_network):
    # Test network construction
    connections = nhd_built_test_network['connections']
    param_df = nhd_built_test_network['param_df']
    assert not param_df.empty
    assert len(connections) > 0

def test_config_validation(validated_config):
    # Test configuration validation
    config_path = Path("test/config.yaml")
    config_data = {"key": "value"}
    validated = validated_config(config_path, config_data)
    assert validated is not None
```

## Best Practices

1. Always use the appropriate fixture for your test case to minimize setup code
2. Use `temporarily_change_dir` context manager when working with paths to ensure you're able to run integration tests successfully
3. Check fixture contents before using them in tests
4. Use fixture combinations when testing complex interactions

## Adding New Fixtures

When adding new fixtures:
1. Follow the existing naming convention
2. Document the fixture structure and purpose
3. Include example usage
4. Add any necessary dependencies
5. Update this README

## Adding New Tests

1. Create new test files in the appropriate category directory
2. Use the naming convention `test_*.py` for test files
3. Document any new fixtures in the relevant conftest.py
