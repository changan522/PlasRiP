# Anaconda.org upload commands

## 1. Build the conda package
```bash
conda install -n base -y -c conda-forge conda-build anaconda-client
conda build conda.recipe
```

## 2. Inspect the generated artifact path
```bash
conda build conda.recipe --output
```

## 3. Login to Anaconda.org
```bash
anaconda login
```

## 4. Upload the package
```bash
anaconda upload $(conda build conda.recipe --output)
```

## Optional channel / label examples
```bash
anaconda upload --label main $(conda build conda.recipe --output)
```

## Optional token-based upload for CI
```bash
anaconda -t "$ANACONDA_API_TOKEN" upload $(conda build conda.recipe --output)
```
