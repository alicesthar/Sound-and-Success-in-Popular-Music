# Sound and Success in Popular Music

This project explores how the acoustic characteristics of popular music have evolved between 2000 and 2018, and how these features relate to commercial success.  
The analysis uses the MusicOSet dataset, combining chart performance data with acoustic features.

## Research Questions
- **EDA**: How have average acoustic characteristics changed over time (2000–2018)?
- **RQ1**: Which acoustic features changed most between 2000-2018?
- **RQ2**: Which acoustic features are associated with being a hit (top 25%)?
- **RQ3**: How do acoustic feature distributions differ between early (2000–2003) and late (2015–2018) eras?

## Dataset
- MusicOSet (2000–2018)  
- Files used:
  - `song_chart.csv`
  - `song_pop.csv`
  - `acoustic_features.csv`

## Repository Structure
```text
├── data_raw/            # Original MusicOSet TSV files
├── outputs/             # CSV summaries (e.g. regression outputs)
├── plots_final/         # Final figures used in the report
├── analysis.R           # Main analysis script
└── README.md            # Project overview and instructions

## How to run

	1.	Place the MusicOSet files in the data_raw/ folder.
	2.	Run Sound-and-Success-in-Popular-Music-25196971.R from start to finish in RStudio.

## Outputs
	•	Figures are saved in plots_final/
	•	Summary tables are saved in outputs/
	•	Key results are also printed to the R console for transparency.