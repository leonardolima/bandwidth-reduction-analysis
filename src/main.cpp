#include <iostream>
#include <filesystem>
#include <algorithm>
#include "peak_mem.h"

int main (void)
{
    std::vector<std::string> exc_files = {"Pipe_to_Table_yxmc.txt", "Legend_Builder_yxmc.txt"};

    std::string folder = "data_sets/";
    for (const auto& file : std::filesystem::directory_iterator(folder))
    {
        std::string file_name = file.path().filename().string();

        if (!(std::find(exc_files.begin(), exc_files.end(), file_name) != exc_files.end()))
        {
            std::cout << "File: " << file_name << std::endl;
            peak_mem("data_sets/" + file_name);
            std::cout << std::endl;
        }
    }

    // peak_mem("data_sets/Base64_Encoder_yxmc.txt");
}
