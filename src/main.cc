#include "sys/types.h"
#include "sys/sysinfo.h"
#include "stdlib.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "myUtils.hh"
#include "Settings.hh"
#include "cosmic_ray_generator_PARMA.hh"

int _mkdir(const char *path)
{
    return ::mkdir(path, 0755); // not sure if this works on mac
}

int main(int argc, char **argv)
{
    const double pi = 3.14159265359;

    srand((unsigned)time(NULL));

    Cosmic_Ray_Generator_PARMA *CR_GEN = new Cosmic_Ray_Generator_PARMA();

    cosmic_ray CR = CR_GEN->generate_One_Cosmic_ray();

    std::cout << "First particle: " << std::endl;
    std::cout << "  Type name: " << CR.type_name << std::endl;
    std::cout << "  Momentum direction (vx,vy,vz; z is vertical): " << CR.momentum_ini.x << " " << CR.momentum_ini.y << " " << CR.momentum_ini.z << std::endl;
    std::cout << "  Energy: " << CR.energy << " MeV" << std::endl;

    const double nadir_angle = std::acos(-CR.momentum_ini.z);

    std::cout << "  Nadir angle: " << nadir_angle * 180.0 / pi << " deg" << std::endl;

    std::ofstream myfile;
    _mkdir("./output");
    myfile.open("output/CR_list.txt");

    const uint nb = 100000;

    std::cout << "\nGenerating particles..." << std::endl;

    myfile << "NAME   ENERGY(MeV)    MOM_DIR_X    MOM_DIR_Y    MOM_DIR_Z    NADIR_ANG(deg)   WEIGHT(%)" << std::endl;

    for (uint ii = 0; ii < nb; ii++)
    {
        cosmic_ray CR = CR_GEN->generate_One_Cosmic_ray();

        myfile << std::scientific;

        myfile << CR.type_name
               << " " << CR.energy
               << " " << CR.momentum_ini.x << " " << CR.momentum_ini.y << " " << CR.momentum_ini.z
               << " " << std::acos(-CR.momentum_ini.z) * 180.0 / pi
               << " " << CR.weight
               << std::endl;
        // output energy in MeV
    }

    myfile.close();
}