# Input Data Specification

* **Format**: CSV or Excel (`.xlsx`).
* **Organization**: One file per metric.
* **Filename**: The filename (excluding extension) must strictly match the
    corresponding entry in `metric_names` (e.g., `Accuracy.csv` for metric
    `Accuracy`).
* **Dimensions**: Rows = Observations (*n*, e.g., Subjects), Columns =
    Datasets (*N*, e.g., Methods).
* **Consistency**: All files must have identical dimensions. Uneven sample sizes
    (missing data) are handled by automatic `NaN` padding (empty cells) to
    ensure a uniform matrix size, and pairwise deletion is applied during
    analysis.
