#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "matrix.hpp"

void readCSV(Matrix &A, const std::string &filename) {
    std::ifstream file;
    file.open(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            row.push_back(std::stod(cell));
        }
        A.cols = row.size();
        for (size_t i = 0; i < A.cols; i++) {
            A(A.rows, i) = row[i];
        }
        A.rows++;
    }
    file.close();
}