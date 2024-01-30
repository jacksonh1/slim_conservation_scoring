from local_conservation_analysis_pipeline import conservation_pipeline

STEPS_TO_RUN = [
    "s1setup_folder",
    "multiprocess_steps",
    "s8calculate_annotations",
    "s9add_annotations2table",
]


if __name__ == '__main__':
    config_file = './params.yaml'
    n_cores = 6
    conservation_pipeline.main(config_file, n_cores, STEPS_TO_RUN)






