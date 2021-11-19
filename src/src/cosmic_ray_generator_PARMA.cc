#include "cosmic_ray_generator_PARMA.hh"

extern "C"
{
    void gen_parma_cr_(const int *,    // rng seed
                       const int *,    // Particle ID (Particle ID, 31:e-, 32:e+, 33:photon)
                       const double *, // emin MeV
                       const double *, // emax MeV
                       const int *,    // nb to generate
                       const double *, // glat deg -90 =< glat =< 90
                       const double *, // glong deg -180 =< glat =< 180
                       const double *, // alt in km
                       double[],       // MeV
                       double[],       // normalized
                       double[],       // normalized
                       double[]        // normalized
    );                                 // number of wanted Particle ID list
                                       //  output and input energies in MeV
    void get_parma_particles_weights_(int[],
                                      const int *,
                                      const double *,
                                      const double *,
                                      const double *,
                                      const double *,
                                      const double *,
                                      double[]);
}

Cosmic_Ray_Generator_PARMA::Cosmic_Ray_Generator_PARMA()
{
    setlocale(LC_ALL, "C"); // just in case

    for (uint i = 0; i < 10; i++)
    {
        std::cout << rand_double() << std::endl;
    }
    std::cout << std::endl;

    // first call to PARMA to generate the list of cosmic rays
    Generate_CR_samples_list_from_PARMA(settings->PART_NAME_TO_SAMPLE);

    std::cout << " " << std::endl;

    std::cout << "CR_GENERATION_ALT:  " << settings->altitude << " km" << std::endl;

    std::cout << "CR_GENERATION_ENER_MIN:  " << settings->CR_GENERATION_ENER_MIN << " MeV" << std::endl;
    std::cout << "CR_GENERATION_ENER_MAX:  " << settings->CR_GENERATION_ENER_MAX << " MeV" << std::endl;

    std::cout << " " << std::endl;

    settings->WEIGHTS.clear();
    settings->WEIGHTS = Calculate_Weights_from_PARMA();

    std::cout << "\nWeights in the selected energy range: " << std::endl;

    int ii = 0;
    for (uint ii = 0; ii < settings->PDG_LIST_ALL.size(); ii++)
    {
        std::cout << "  " << settings->NAMES[ii] << ": " << settings->WEIGHTS[ii] << std::endl;

        if (settings->NAMES[ii] == settings->PART_NAME_TO_SAMPLE)
        {
            selected_weight = settings->WEIGHTS[ii] * 100.0;
        }
    }

    std::cout << std::endl;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Cosmic_Ray_Generator_PARMA::find_parma_ID_from_PDG(const int &PDG_in)
{

    const int pdg_phot = 22;
    const int pdg_elec = 11;
    const int pdg_posi = -11;
    const int pdg_muN = 13;
    const int pdg_muP = -13;
    const int pdg_neut = 2112;
    const int pdg_prot = 2212;
    const int pdg_alpha = 9999;

    const int parma_phot = 33;
    const int parma_elec = 31;
    const int parma_posi = 32;
    const int parma_muN = 30;
    const int parma_muP = 29;
    const int parma_neut = 0;
    const int parma_prot = 1;
    const int parma_alpha = 2;

    switch (PDG_in)
    {
    case pdg_phot:
        return parma_phot;

    case pdg_elec:
        return parma_elec;

    case pdg_posi:
        return parma_posi;

    case pdg_muN:
        return parma_muN;

    case pdg_muP:
        return parma_muP;

    case pdg_neut:
        return parma_neut;

    case pdg_prot:
        return parma_prot;

    case pdg_alpha:
        return parma_alpha;

    default:
        std::cout << "Error: not a valid PDG number in find_parma_ID_from_PDG in cosmic_ray_generator_PARMA.cc" << std::endl;
        std::abort();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

cosmic_ray Cosmic_Ray_Generator_PARMA::generate_One_Cosmic_ray()
{

    cosmic_ray cosmic_part = {};

    //            sampled_set = read_particles[index_sampling_part];
    cosmic_ray_parma_output sampled_set = sample_One_CR_from_PARMA();

    ThreeVector position_ini = sample_CR_secondary_position(settings->altitude);

    // from PARMA OUTPUT:
    //   cos(theta)=1,  indicates  the  vertical  downward  direction,
    //   while  90  degree,  i.e.  cos(theta)=0, indicates the horizontal direction.
    ThreeVector momentum_ini = CR_direction_rand_sample(sampled_set.u, sampled_set.v, sampled_set.w);
    // multiplication by -1 is important, to make sure that when sampled_set.cos_zenith_angle is 1, the particle is sampled vertical downward

    const double time = 0.0;
    const double energy = sampled_set.energy; // MeV

    cosmic_part.energy = energy;
    cosmic_part.time = time;
    cosmic_part.momentum_ini = momentum_ini;
    cosmic_part.position_ini = position_ini;
    cosmic_part.type_name = settings->PART_NAME_TO_SAMPLE;

    cosmic_part.weight = selected_weight;

    return cosmic_part;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Following angular distributions compuer from PARMA code
ThreeVector Cosmic_Ray_Generator_PARMA::CR_direction_rand_sample(const double &u, const double &v, const double &w)
{
    ThreeVector momentum_ini{};

    // if cos_Sampled == 1 (zenith) then direction should be (0,0,1)

    momentum_ini.x = u;
    momentum_ini.y = v;
    momentum_ini.z = w;

    const double norm = std::sqrt(momentum_ini.x * momentum_ini.x + momentum_ini.y * momentum_ini.y + momentum_ini.z * momentum_ini.z);

    momentum_ini.x = momentum_ini.x / norm;
    momentum_ini.y = momentum_ini.y / norm;
    momentum_ini.z = momentum_ini.z / norm;

    return momentum_ini;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ThreeVector Cosmic_Ray_Generator_PARMA::sample_CR_secondary_position(const double &altitude)
{
    double r1 = rand_double();
    double r2 = rand_double();

    ThreeVector position = {(r1 - 0.5) * settings->CR_SAMPLING_XZ_HALF_SIZE * 2.0 * km,
                            (r2 - 0.5) * settings->CR_SAMPLING_XZ_HALF_SIZE * 2.0 * km,
                            altitude * km};

    return position;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// uses cumulative distribution sampling
cosmic_ray_parma_output Cosmic_Ray_Generator_PARMA::sample_One_CR_from_PARMA()
{

    if (counter >= nb_to_generate - 3) // if all the particles generated from Parma have been already used, generate new ones
    {
        Generate_CR_samples_list_from_PARMA(settings->PART_NAME_TO_SAMPLE);

        std::cout << "Generated " << nb_to_generate << " new random cosmic ray particles." << std::endl;
    }

    const double eRand = output_energies[counter];
    const double u_rand = output_u[counter];
    const double v_rand = output_v[counter];
    const double w_rand = output_w[counter];
    counter++;

#ifndef NDEBUG // if debug mode, some sanity checks

    if ((eRand < min_cr_ener) || (eRand > max_cr_ener))
    {
        std::cout << "Energy is out of range. Aborting." << std::endl;
        std::abort();
    }

#endif // ifndef NDEBUG

    cosmic_ray_parma_output spld_set{eRand, u_rand, v_rand, w_rand};

    return spld_set;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Cosmic_Ray_Generator_PARMA::Generate_CR_samples_list_from_PARMA(const std::string &part_name)
{

    const int nb = nb_to_generate;

    int type_index;

    if (part_name == "proton")
    {
        type_index = find_parma_ID_from_PDG(pdg_prot);
    }
    else if (part_name == "photon")
    {
        type_index = find_parma_ID_from_PDG(pdg_phot);
    }
    else if (part_name == "electron")
    {
        type_index = find_parma_ID_from_PDG(pdg_elec);
    }
    else if (part_name == "positron")
    {
        type_index = find_parma_ID_from_PDG(pdg_posi);
    }
    else if (part_name == "muonN")
    {
        type_index = find_parma_ID_from_PDG(pdg_muN);
    }
    else if (part_name == "muonP")
    {
        type_index = find_parma_ID_from_PDG(pdg_muP);
    }
    else if (part_name == "neutron")
    {
        type_index = find_parma_ID_from_PDG(pdg_neut);
    }
    else
    {
        std::cout << "Wrong Sample Particle name. It should be photon, electron, positron, muonN, muonP, neutron or proton" << std::endl;
        std::abort();
    }

    const int the_seed = myUtils::generate_a_unique_ID_int32();
    // seed, type, emin, emax, glat, glong, Alti, ener_out, u_out, v_out, w_out
    gen_parma_cr_(&the_seed,
                  &type_index,
                  &min_cr_ener,
                  &max_cr_ener,
                  &nb,
                  &glat,
                  &glong,
                  &alt,
                  output_energies,
                  output_u,
                  output_v,
                  output_w);

    counter = 0;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> Cosmic_Ray_Generator_PARMA::Calculate_Weights_from_PARMA()
{
    double output_w[n_types_for_weights];

    int *parmaID_list_wanted2 = &parmaID_list_ALL[0]; // trick to transform a C++ vector into a C array

    get_parma_particles_weights_(parmaID_list_wanted2,
                                 &n_types_for_weights,
                                 &min_cr_ener,
                                 &max_cr_ener,
                                 &alt,
                                 &glat,
                                 &glong,
                                 output_w);

    //get_parma_particles_weights(type_list, ntypes, emin, emax, alt, glat, glong, weight_list)

    std::vector<double> output;
    output.clear();
    for (int ii = 0; ii < n_types_for_weights; ++ii)
    {
        output.push_back(output_w[ii]);
    }

    return output;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Cosmic_Ray_Generator_PARMA::rand_double()
{
    return ((double)rand() / (double)RAND_MAX);
}