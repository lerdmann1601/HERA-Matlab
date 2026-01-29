# Advanced Usage (Developer Mode)

Developers can call `HERA.run_ranking` directly with data matrices, bypassing
file I/O. This is useful for integration into other pipelines or simulation
studies.

## Syntax

```matlab
results = HERA.run_ranking(userInput);
```

## Input Structure (`userInput`)

Instead of `folderPath`, provide `custom_data`:

```matlab
% 1. Prepare Data
% Cell array of (N_Subjects x N_Methods) matrices
data_m1 = randn(50, 5); 
data_m2 = randn(50, 5);
custom_data = {data_m1, data_m2};

% 2. Configure User Input
userInput = struct();
userInput.custom_data = custom_data;
userInput.metric_names = {'Accuracy', 'Runtime'}; % Must match data count
userInput.dataset_names = {'Method A', 'Method B', 'Method C', 'Method D',
    'Method E'};
userInput.ranking_mode = 'M1_M2';
userInput.output_dir = pwd;

% 3. Run
results = HERA.run_ranking(userInput);
```

## Output (`results` struct)

* `.final_rank`: Final ranking vector.
* `.d_vals_all`: Effect sizes (Cliff's Delta).
* `.p_vals_all`: Raw p-values.
* `.ci_d_all`: Confidence intervals.

For the full reference of the `results` struct, see [Results Structure Reference](https://lerdmann1601.github.io/HERA-Matlab/Results_Structure_Reference).
