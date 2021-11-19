#pragma once

#include <vector>
#include <string>
// INITIAL PARTICLES CAN BE: PHOTON, ELECTRON, POSITRON, MUONN, MUONP, NEUTRON, PROTON

typedef unsigned int uint;

class Settings
{
private:
    Settings() = default; // Private so that it can not be called

    Settings(Settings const &) {}

    // copy constructor is private
    // assignment operator is private
    static Settings *instance;

public:
    static Settings *getInstance();

public:
    const uint i_p = 0;
    const uint i_e = 1;
    const uint i_ep = 2;
    const uint i_mn = 3;
    const uint i_mp = 4;
    const uint i_ne = 5;
    const uint i_pr = 6;
    const uint i_al = 7;

    const std::string NAMES[8]{"photon", "electron", "positron", "muonN", "muonP", "neutron", "proton","alpha"};

    const std::string PART_NAME_TO_SAMPLE = "proton";

    enum PDG_nb : int
    {
        pdg_phot = 22,
        pdg_elec = 11,
        pdg_posi = -11,
        pdg_muN = 13,
        pdg_muP = -13,
        pdg_neut = 2112,
        pdg_prot = 2212,
        pdg_alpha = 9999
    };

    // list of PDG number of particles that can be generated

    const std::vector<int> PDG_LIST_ALL = {pdg_phot, pdg_elec, pdg_posi, pdg_muN, pdg_muP, pdg_neut, pdg_prot, pdg_alpha};

    std::vector<double> WEIGHTS;

    const double CR_GENERATION_ENER_MIN = 1000;   // MeV
    const double CR_GENERATION_ENER_MAX = 100000.0; // MeV

    const double altitude = 50.0; // km
    const double CR_SAMPLING_XZ_HALF_SIZE = 1;

    const double latitude = 78.0;
    const double longitude = 0.0;
};