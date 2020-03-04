#include <iostream>
#include <filesystem>
#include "peak_mem.h"

int main (void)
{
    std::string folder = "data_sets/";
    for (const auto& file : std::filesystem::directory_iterator(folder))
    {
        std::string file_name = file.path().filename().string();

        std::cout << "File: " << file_name << std::endl;
        peak_mem("data_sets/" + file_name);
        std::cout << std::endl;
    }
}
