Main codes for the paper 'Apnoea suppresses brain activity in infants' (bioRxiv DOI: https://doi.org/10.1101/2024.02.16.580547)

Repository contains codes:
- compute_Hilbert_and_stats.m: gets EEG amplitudes for apnoeas and breathing pauses. Also performs time-frequency statistics to check if these amplitudes are different compared with control periods.
- linear_models_eeg_physiology.m: linear models to test whether the EEG amplitude changes are associated with (1) heart rate, (2) oxygen saturation, (3) pause duration, (4) age, and (5) sleep state.

Note that the apnoeas/breathing pauses should have been computed at this point. The algorithm to determine the inter-breath intervals from the impedance pneumographic signal is available at: https://gitlab.com/paediatric_neuroimaging/identify_ibi_from_ip.git. Usage instructions can be found at: https://doi.org/10.1136/bmjresp-2021-001042.
