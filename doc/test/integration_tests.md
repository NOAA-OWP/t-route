# Integration Tests

This directory details the integration tests within troute

## Quick Start

```bash
# Run all integration tests
pytest test/

# Run tests from a specific category
pytest test/troute-nwm
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

### Routing Tests (`test/troute-routing/`)

**Usage:**
- Integration tests for NHD and HYFeature util functions

### End to End Examples:
##### `test/LowerColorado_TX/`
##### `test/LowerColorado_TX_v4/`
##### `test/LowerColorado_TX_HYFeatures_v22/`

## Adding New Tests

1. Create new test files in the appropriate category directory
2. Use the naming convention `test_*.py` for test files
3. Document any new fixtures in the relevant conftest.py
