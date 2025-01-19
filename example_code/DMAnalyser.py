from ednabp import SampleData
from ednabp import DMAnalyser

def main():
    # sample_data = SampleData()
    # sample_data.unpickle_data("CJ_test/eDNA_samples_2025-01-08.pkl")
    # dm_analyser = DMAnalyser(sampledata=sample_data, no_verbose=False)
    # dm_analyser.write_abundance_table(
    #     save_dir="CJ_test/",
    #     taxa_level="species",
    #     process=None,
    #     sample_id_list=None
    # )
    # dm_analyser.write_richness_table(save_dir="CJ_test/",
    #                                  taxa_level="kingdom",
    #                                  unit_level="species",
    #                                  sample_id_list=None)

    # dm_analyser = DMAnalyser()
    # dm_analyser.plot_heatmap(
    #     csv_path="CJ_test/family_species_richness.csv",
    #     taxa_column="family",
    #     metric_column="richness",
    #     x_categories=["Year", "Month", "Site"],
    # )

    # dm_analyser.plot_barchart(csv_path="TI_test/family_norm_abundance.csv",
    #                           taxa_column="family",
    #                           metric_column="abundance",
    #                           save_dir=None,
    #                           overwrite=True)
    
    # dm_analyser.plot_contour(csv_path="TI_test/species_log_abundance.csv",
    #                          shp_path="TI_test/gis/kueishan.shp",
    #                          metric_column="abundance",
    #                          save_dir="TI_test",
    #                          overwrite=True)
    return

if __name__ == "__main__":
    main()