# Urine_predict

This is an Rshiny app allowing the user to preprocess fullscan mzML files (one per sample) in order to run urinary tract infection prediction on the data. The preprocessing function takes precalculated target data and reduced target data (both csvs are provided) which provides the necessary information to extract target peptides from the fullscan data. Prediction will run automatically after preprocessing. 

The following figure illustrates the different steps performed during preprocessing: 

<img width="2000" height="1500" alt="ms_pipeline_clean" src="https://github.com/user-attachments/assets/601353c7-3861-41da-a381-6ab974963eeb" />


Should you wish to run prediction on a previously preprocessed file, the model ready csv can be directly uploaded and prediction will be performed using the pretrained random forest models contained in models.Rdata. 
