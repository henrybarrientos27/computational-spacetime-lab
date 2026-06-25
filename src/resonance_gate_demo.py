from pathlib import Path
import pandas as pd

DATA_PATH = Path("data/sample_scan.csv")

def load_scan(path: Path = DATA_PATH) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing data file: {path}")
    return pd.read_csv(path)

def add_correction_term(df: pd.DataFrame) -> pd.DataFrame:
    required = {"frequency_hz", "amplitude", "baseline_response"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df.copy()

    threshold = df["amplitude"].median()
    gated = df["amplitude"] > threshold

    df["gate_active"] = gated
    df["correction"] = 0.0
    df.loc[gated, "correction"] = (
        0.05 * df.loc[gated, "baseline_response"] * df.loc[gated, "amplitude"]
    )

    df["corrected_response"] = df["baseline_response"] + df["correction"]
    return df

def main() -> None:
    df = load_scan()
    corrected = add_correction_term(df)
    print(corrected.to_string(index=False))

if __name__ == "__main__":
    main()
