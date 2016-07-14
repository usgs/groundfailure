#!/usr/bin/env python


def correct_config_filepaths(config):
    # Pull out the input filepath
    input_path = config['input']['folder']
    # Pull all other filepaths that need editing
    paths_to_correct = []
    paths_to_correct.append(config['mapdata']['dem']['file'])
    paths_to_correct.append(config['mapdata']['roads']['folder'])
    paths_to_correct.append(config['mapdata']['cities']['file'])
    paths_to_correct.append(config['mapdata']['ocean']['file'])
    paths_to_correct.append(config['mechanistic_models']['godt_2008']['layers']['cohesion']['file'])
    paths_to_correct.append(config['mechanistic_models']['godt_2008']['layers']['friction']['file'])
    paths_to_correct.append(config['mechanistic_models']['godt_2008']['layers']['slope']['filepath'])
    paths_to_correct.append(config['mechanistic_models']['classic_newmark']['layers']['cohesion']['file'])
    paths_to_correct.append(config['mechanistic_models']['classic_newmark']['layers']['friction']['file'])
    paths_to_correct.append(config['mechanistic_models']['classic_newmark']['layers']['slope']['file'])
    paths_to_correct.append(config['mechanistic_models']['classic_newmark']['layers']['watertable']['file'])
    paths_to_correct.append(config['mechanistic_models']['hazus']['layers']['susceptibility']['file'])
    paths_to_correct.append(config['logistic_models']['nowicki_2015']['layers']['slope'])
    paths_to_correct.append(config['logistic_models']['nowicki_2015']['layers']['rock'])
    paths_to_correct.append(config['logistic_models']['nowicki_2015']['layers']['landcover'])
    paths_to_correct.append(config['logistic_models']['nowicki_2015']['layers']['precip'])
    paths_to_correct.append(config['logistic_models']['nowicki_2015']['layers']['cti'])
    paths_to_correct.append(config['logistic_models']['nowicki_2015']['layers']['elev'])
    paths_to_correct.append(config['logistic_models']['zhu_2014']['layers']['vs30'])

    x = len(paths_to_correct)

    # for loop to add together the strings
    corrected_paths = []
    paths_to_correct.append(config['logistic_models']['zhu_2014']['layers']['cti'])
    for key in range(0, x+1):
        corrected_paths.append(str(input_path+paths_to_correct[key]))

    # debugging - delete later
    print(corrected_paths)

    # Replace values in config from corrected paths
    config['mapdata']['dem']['file'] = corrected_paths[0]
    config['mapdata']['roads']['folder'] = corrected_paths[1]
    config['mapdata']['cities']['file'] = corrected_paths[2]
    config['mapdata']['ocean']['file'] = corrected_paths[3]
    config['mechanistic_models']['godt_2008']['layers']['cohesion']['file'] = corrected_paths[4]
    config['mechanistic_models']['godt_2008']['layers']['friction']['file'] = corrected_paths[5]
    config['mechanistic_models']['godt_2008']['layers']['slope']['filepath'] = corrected_paths[6]
    config['mechanistic_models']['classic_newmark']['layers']['cohesion']['file'] = corrected_paths[7]
    config['mechanistic_models']['classic_newmark']['layers']['friction']['file'] = corrected_paths[8]
    config['mechanistic_models']['classic_newmark']['layers']['slope']['file'] = corrected_paths[9]
    config['mechanistic_models']['classic_newmark']['layers']['watertable']['file'] = corrected_paths[10]
    config['mechanistic_models']['hazus']['layers']['susceptibility']['file'] = corrected_paths[11]
    config['logistic_models']['nowicki_2015']['layers']['slope'] = corrected_paths[12]
    config['logistic_models']['nowicki_2015']['layers']['rock'] = corrected_paths[13]
    config['logistic_models']['nowicki_2015']['layers']['landcover'] = corrected_paths[14]
    config['logistic_models']['nowicki_2015']['layers']['precip'] = corrected_paths[15]
    config['logistic_models']['nowicki_2015']['layers']['cti'] = corrected_paths[16]
    config['logistic_models']['nowicki_2015']['layers']['elev'] = corrected_paths[17]
    config['logistic_models']['zhu_2014']['layers']['vs30'] = corrected_paths[18]
    config['logistic_models']['zhu_2014']['layers']['cti'] = corrected_paths[19]

    return config
