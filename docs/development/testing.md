# Testing

## Running tests

Test cases are located under `test/`. Each subdirectory is a self-contained case with input files and an expected output for comparison.

```bash
cd test/mono-disperse
../../bin/ICE
```

## Adding a new test case

1. Create a new directory under `test/`.
2. Add the required input files (`input.ini`, `ic/`, `bc/`, `grid/`).
3. Document the expected solution and the analytical or reference result in a `README.md` inside the case directory.

---
