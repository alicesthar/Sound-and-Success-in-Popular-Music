# Sound and Success in Popular Music

This project examines how the acoustic characteristics of popular music have evolved between 2000 and 2018, and how these characteristics are associated with commercial success. Using the MusicOSet dataset, the analysis combines exploratory data analysis with statistical modelling to investigate long-term trends in musical sound and differences between historical chart eras.

The project has been developed as part of the IJC437 Introduction to Data Science module in the MSc Data Science programme at the University of Sheffield.

## Research Questions
- **EDA**: How have average acoustic characteristics changed over time (2000–2018)?
- **RQ1**: Which acoustic features changed most between 2000-2018?
- **RQ2**: Which acoustic features are associated with being a hit (top 25%)?
- **RQ3**: How do acoustic feature distributions differ between early (2000–2003) and late (2015–2018) eras?

## Key Findings
-	Several acoustic features exhibit clear long-term trends between 2000 and 2018, indicating gradual changes in the sound of popular chart music over time. However, these changes are not uniform across features, suggesting that musical evolution has occurred along specific dimensions rather than as a broad shift in all acoustic properties.
-	Linear trend analysis shows that valence and energy experience the largest overall declines across the study period, while tempo and acousticness show modest increases. In contrast, danceability and loudness remain relatively stable, indicating that not all aspects of musical structure change at the same rate.
-	Logistic regression analysis identifies several acoustic features that are statistically associated with chart success. Higher loudness and danceability are linked to increased odds of a song being classified as a hit (top 25%), whereas higher energy and acousticness are associated with lower hit odds. Tempo shows no strong association with hit status in the multivariate model.
-	Comparisons between early (2000–2003) and late (2015–2018) eras reveal clear differences in the distributions of multiple acoustic features. Songs from the late era tend to exhibit lower valence and energy, alongside slightly higher tempo and acousticness, highlighting structural differences in the acoustic profiles of popular music across the two periods.
	
## Dataset
This project uses the MusicOSet dataset (2000–2018), an open dataset designed for music data mining research.

- MusicOSet (2000–2018)  
- Files used:
  - `song_chart.csv`
  - `song_pop.csv`
  - `acoustic_features.csv`
All datasets are used in their original form and are placed in the data_raw/ directory.

## R Code

All analysis is conducted in R and implemented in a single, fully commented script (Sound-and-Success-in-Popular-Music-25196971.R).
The code follows a consistent style and naming convention and is structured to reflect the analytical workflow:
	1.	Data loading and auditing
	2.	Data cleaning and feature scaling
	3.	Exploratory data analysis
	4.	Statistical modelling for each research question
	5.	Figure and table generation

Intermediate results and audit information are printed to the R console to support transparency and reproducibility.

## Instructions for Running the Code

To reproduce the analysis:
	1.	Clone or download this repository.
	2.	Place the required MusicOSet data files in the data_raw/ directory.
	3.	Open Sound-and-Success-in-Popular-Music-2519697.R in RStudio.
	4.	Ensure the required R packages are installed (the script automatically checks and installs missing packages).
	5.	Run the script sequentially from top to bottom.

All figures are saved automatically to the plots_final/ directory, and summary tables are saved to the outputs/ directory.

## Repository Structure
```text
├── data_raw/            # Original MusicOSet TSV files
├── outputs/             # CSV summaries (e.g. regression outputs)
├── plots_final/         # Final figures used in the report
├── analysis.R           # Main analysis script
└── README.md            # Project overview and instructions
