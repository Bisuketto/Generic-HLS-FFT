/**
 * @file DataVector.cpp
 * @author Bordeaux INP, Bertrand LE GAL (bertrand.legal@ims-bordeaux.fr) [http://legal.vvv.enseirb-matmeca.fr]
 * @brief This file is part of a LDPC library for realtime LDPC decoding on processor cores.
 * @version 0.0.0
 * @date 2020
 * 
 * @license This source is released under the GNU GENERAL PUBLIC LICENSE Version 3
 * 
 */

#include "DataVector.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


DataVector::DataVector(const DataVector &in)
{
    I.resize(in.I.size());
    Q.resize(in.Q.size());
}

DataVector::DataVector(const std::string filename)
{
    std::ifstream ifs(filename);
    if (ifs.is_open() == false)
    {
        std::cout << "(EE) Erreur critique à l'ouverture du fichier (" << filename << ")" << std::endl;
        std::cout << "(EE) L'execution du programme doit s'interrompre..." << std::endl;
        exit( EXIT_FAILURE );
    }

    std::string line;
    while (!ifs.eof()) {
        getline(ifs, line);
        if (line.size() == 0) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);
        std::string tmp;
        getline(iss, tmp, ' '); // ech number
        assert(I.size() == std::stoi(tmp));
        getline(iss, tmp, ' '); // ech number
        I.push_back(std::stof(tmp));
        getline(iss, tmp, ' '); // ech number
        Q.push_back(std::stof(tmp));
    }
}

bool DataVector::save(const std::string filename)
{
    std::string dir = filename.substr(0, filename.find_last_of("/"));
    int result = mkdir(dir.c_str(), 0777);
    // if(result != 0){
    //     // printf("Unable to create directory : %s\n", dir.c_str());
    //     perror("Unable to create result directory");
    //     return false;
    // }    

    FILE *fo = fopen(filename.c_str(), "w+");
    if (fo == NULL)
    {
        printf("(EE) Erreur à l'ouverture du fichier de sortie (%s).\n", filename.c_str());
        return false;
    }
    fprintf(fo, "indice real imag mag\n");
    std::complex<double> vabs;
    for (int i = 0; i < I.size(); i += 1) {
        vabs.real(I[i]); 
        vabs.imag(Q[i]);
        fprintf(fo, "%d %+23.23f %+23.23f %+23.23f\n", i, I[i], Q[i], 10*log10(std::abs(vabs)));
    }
    fclose(fo);

    return true;
}
