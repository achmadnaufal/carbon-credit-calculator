# Carbon Credit Calculator

Python tool for estimating carbon sequestration and credits from tree planting projects

## Features
- Data ingestion from CSV/Excel input files
- Automated analysis and KPI calculation
- Summary statistics and trend reporting
- Sample data generator for testing and development

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

```python
from src.main import CarbonCreditCalculator

analyzer = CarbonCreditCalculator()
df = analyzer.load_data("data/sample.csv")
result = analyzer.analyze(df)
print(result)
```

## Data Format

Expected CSV columns: `tree_id, species, dbh_cm, height_m, age_years, plot_id`

## Project Structure

```
carbon-credit-calculator/
├── src/
│   ├── main.py          # Core analysis logic
│   └── data_generator.py # Sample data generator
├── data/                # Data directory (gitignored for real data)
├── examples/            # Usage examples
├── requirements.txt
└── README.md
```

## License

MIT License — free to use, modify, and distribute.
