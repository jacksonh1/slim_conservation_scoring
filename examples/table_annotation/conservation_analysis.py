from local_conservation_analysis_pipeline import conservation_pipeline



if __name__ == '__main__':
    config_file = './params.yaml'
    n_cores = 6
    conservation_pipeline.main(config_file, n_cores)






