#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <limits>

typedef std::vector<unsigned> container;

int main (void)
{
    // nan
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    // filename of labels file
    std::string label_file =
        "../data/Myconnectome/labels_hammer_myconnectome_new2.csv";

    // open labels file
    std::ifstream stream_label (label_file);

    // check if file opening works correctly
    if (!stream_label)
    {
        std::cerr << "No label file => I cowardly refuse to go on"
                  << std::endl;
        return 1; //EXIT_FAILURE
    }

    // number of regions
    constexpr unsigned n_regions = 83;

    // number of time samples for each node
    constexpr unsigned n_time_samples = 518;

    // the map region_to_index stores, for every region, the indices of nodes
    // that belong to it
    std::map<unsigned, container> region_to_index;

    // Build the map
    unsigned counter = 0;
    unsigned idx = 0;
    unsigned n_nodes = 0;
    while (stream_label >> idx)
    {
        if (idx >= 1 && idx <= n_regions)
        {
            region_to_index[idx].push_back (counter++);
        }
        n_nodes++;
    }
    // folder names
    std::string data_folder = "../data/Myconnectome/time_series/";
    std::string result_folder = "../data/Myconnectome/time_means/";

    for (unsigned ses = 11; ses <= 104; ses++)
    {
        if (ses != 42 && ses != 52 && ses != 90 && ses != 80)
        {
            // Build file name
            std::ostringstream data_file_names;
            std::ostringstream result_file_names;

            data_file_names << data_folder << "ses";
            result_file_names << result_folder << "ses";
            if (ses < 100)
            {
                data_file_names << "0";
                result_file_names << "0";
            }
            data_file_names << ses << "_time_series.csv";
            result_file_names << ses << "_time_means.csv";

            std::cout << "processing session  ses" << ses << std::endl;

            std::string series_file = data_file_names.str();
            std::string result_file = result_file_names.str();

            // open files
            std::ifstream stream_series (series_file);
            std::ofstream stream_result (result_file);

            // check if file opening works correctly
            if (!stream_series)
            {
                std::cerr << "No series file => I cowardly refuse to go on"
                          << std::endl;
                return 1; //EXIT_FAILURE
            }

            if (!stream_result)
            {
                std::cerr << "No result file => I cowardly refuse to go on"
                          << std::endl;
                return 1; //EXIT_FAILURE
            }

            // matrix coming from the time series observations
            // number of rows: the number of nodes in the mesh
            // number of columns: the number of time samples
            std::vector<std::vector<double>> matrix (n_nodes,
                                                     std::vector<double> (n_time_samples, 0.) );

            // Build the matrix
            std::string row;
            unsigned row_idx = 0;
            while ( getline (stream_series, row) )
            {
                std::istringstream row_stream (row);
                unsigned j = 0;
                std::string elem;
                while (getline (row_stream, elem, ',') && (j < n_time_samples) )
                {
                    matrix[row_idx][j] = std::stod (elem);
                    j++;
                }
                row_idx++;
            }

            // mean_matrix: n_regions x n_time_samples matrix, storing in each row
            // the means of observations in matrix belonging to the relative region
            std::vector<std::vector<double>> mean_matrix (n_regions,
                                                          std::vector<double> (n_time_samples, 0.) );

            // Build mean_matrix
            for (unsigned reg_idx = 1; reg_idx <= n_regions; reg_idx++)
            {
                std::map<unsigned, container>::const_iterator it =
                    region_to_index.find (reg_idx);
                if (it != region_to_index.cend() )
                {
                    const container& nodes_in_region = it->second;
                    unsigned region_size = nodes_in_region.size();
                    for (unsigned k = 0; k < n_time_samples; k++)
                    {
                        for (unsigned j = 0; j < region_size; j++)
                        {
                            mean_matrix[reg_idx - 1][k] += matrix[nodes_in_region[j]][k];
                        }
                        mean_matrix[reg_idx - 1][k] /= region_size;
                    }
                }
                else
                {
                    for (unsigned k = 0; k < n_time_samples; k++)
                    {
                        mean_matrix[reg_idx - 1][k] = nan;
                    }
                }
            }

            // Print mean_matrix on csv
            for (unsigned j = 0; j < n_regions; j++)
            {
                for (unsigned k = 0; k < n_time_samples; k++)
                {
                    char del = (k < n_time_samples - 1) ? ',' : '\n';
                    stream_result << mean_matrix[j][k] << del;
                }
            }
        }
    }
    std::cout << "Work finished" << std::endl;
    return 0; //EXIT_SUCCESS
}