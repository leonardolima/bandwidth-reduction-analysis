#include <iostream>
#include <filesystem>
#include <algorithm>
#include "io.h"
#include "test.h"
#include "basic.h"

int main (void)
{
    std::vector<std::string> exc_files = {"Pipe_to_Table_yxmc.txt", "Legend_Builder_yxmc.txt"};
    std::vector<std::pair<int, std::string>> pairs;
    std::string folder = "data_sets/";

    for (const auto& file : std::filesystem::directory_iterator(folder))
    {
        std::string file_name = file.path().filename().string();

        if (!(std::find(exc_files.begin(), exc_files.end(), file_name) != exc_files.end()))
        {
            int N = read_N_from_file(file_name);
            std::pair<int, std::string> pair = std::make_pair(N, file_name);
            pairs.push_back(pair);
        }
    }

    std::sort(pairs.begin(), pairs.end());

    for(std::vector<std::pair<std::string, int>>::size_type i = 0; i < pairs.size(); ++i)
    {
        std::cout << "File: " << pairs[i].second << std::endl;
        test_adapted_cuthill_mckee("data_sets/" + pairs[i].second);
        std::cout << std::endl;
    }

    std::cout << "=============================================" << std::endl;

    for(std::vector<std::pair<std::string, int>>::size_type i = 0; i < pairs.size(); ++i)
    {
        std::cout << "File: " << pairs[i].second << std::endl;
        test_levels("data_sets/" + pairs[i].second);
        std::cout << std::endl;
    }

    // for(std::vector<std::pair<std::string, int>>::size_type i = 0; i < pairs.size(); ++i)
    // {
    //     std::cout << "File: " << pairs[i].second << std::endl;
    //     test_topological("data_sets/" + file_name);
    //     std::cout << std::endl;
    // }

    // test_all("data_sets/Base64_Encoder_yxmc.txt");
    //test_all("data_sets/example.txt");
}
