#pragma once

#include <locale.h>
#include <numeric>
#include <stdlib.h>
#include <random>
#include "myUtils.hh"
#include "Settings.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct cosmic_ray_parma_output
{
    double energy;
    double u;
    double v;
    double w;
};

struct ThreeVector
{
    double x;
    double y;
    double z;
};

struct cosmic_ray
{
    ThreeVector momentum_ini;
    ThreeVector position_ini;
    double time;
    double energy;
    double weight;
    std::string type_name;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Cosmic_Ray_Generator_PARMA
{
public:
    Cosmic_Ray_Generator_PARMA();

    ~Cosmic_Ray_Generator_PARMA() = default;

    cosmic_ray generate_One_Cosmic_ray();

    std::vector<double> Calculate_Weights_from_PARMA();

    const double MeV = 1.0;
    const double pi = 3.14159265359;
    const double km = 1.0e-6;

private:
    double rand_double();

    const int i_phot = 0;
    const int i_elec = 1;
    const int i_posi = 2;
    const int i_muN = 3;
    const int i_muP = 4;
    const int i_neut = 5;
    const int i_prot = 6;
    const int i_alpha = 7;

    const int pdg_phot = 22;
    const int pdg_elec = 11;
    const int pdg_posi = -11;
    const int pdg_muN = 13;
    const int pdg_muP = -13;
    const int pdg_neut = 2112;
    const int pdg_prot = 2212;
    const int pdg_alpha = 9999;

    Settings *settings = Settings::getInstance();

    enum Parma_ID : int
    {
        parma_phot = 33,
        parma_elec = 31,
        parma_posi = 32,
        parma_muN = 30,
        parma_muP = 29,
        parma_neut = 0,
        parma_prot = 1,
        parma_alpha = 2
    };

    std::vector<int> parmaID_list_ALL{parma_phot,
                                      parma_elec,
                                      parma_posi,
                                      parma_muN,
                                      parma_muP,
                                      parma_neut,
                                      parma_prot,
                                      parma_alpha};

    ////////////////////////////////////////////////

    double min_cr_ener = settings->CR_GENERATION_ENER_MIN; // MeV
    double max_cr_ener = settings->CR_GENERATION_ENER_MAX; // MeV
    double min_cosAng = -1.0;
    double max_cosAng = 1.0;

    double glat = settings->latitude;
    double glong = settings->longitude;
    double alt = settings->altitude;

    // the bigger the numbers, more precise will be the sampling, but the more memory it will take

    const static int nb_to_generate = 10000;
    int counter = 0;           // indexing of the sampled cosmic rays
    int seed_cr_smpl = 123456; // dummy value, random seed of CR generator

    double output_energies[nb_to_generate]{};
    double output_u[nb_to_generate]{};
    double output_v[nb_to_generate]{};
    double output_w[nb_to_generate]{};

    double selected_weight = 0.0;

    const int n_types_for_weights = 8;

    // Functions

    ThreeVector CR_direction_rand_sample(const double &u, const double &v, const double &w);

    ThreeVector sample_CR_secondary_position(const double &altitude);

    cosmic_ray_parma_output sample_One_CR_from_PARMA();

    void generate_output_for_test();

    void generate_output_for_test_momentum(const ThreeVector &mom);

    void Generate_CR_samples_list_from_PARMA(const std::string &part_name);

    int find_parma_ID_from_PDG(const int &PDG_in);

    ///
    std::string name_outFile_mom = "./tests/cr_sampl_test_mom.txt";
    std::string filename_cr_sampling_test = "./tests/cr_sampl_test_";
};
