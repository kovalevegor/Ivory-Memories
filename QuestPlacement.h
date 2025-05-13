#pragma once
#include "Knot.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <random>
#include <unordered_map>

struct Chromosome {
    std::vector<int> node_ids; // node_ids[i] - ID узла для тега i
};

class QuestPlacement {
private:
    const std::vector<Knot>& knots;
    const std::unordered_map<Knot, std::vector<Knot>, std::hash<Knot>>& adjacency_list;
    std::vector<std::string> tag_list; // Список тегов, например, ["spawn", "quest_giver_heron"]
    std::map<std::string, int> desired_distances; // Желаемые дистанции от spawn_knot
    std::vector<std::tuple<std::string, std::string, int>> min_distance_constraints; // Минимальные расстояния между тегами
    int population_size;
    int generations;
    float mutation_rate;

    std::vector<Chromosome> initialize_population(int num_nodes, int num_tags, std::mt19937& rng);
    float compute_fitness(const Chromosome& chrom, const std::vector<int>& distances, const std::vector<std::vector<int>>& dist_matrix);
    Chromosome crossover(const Chromosome& parent1, const Chromosome& parent2, int num_nodes, std::mt19937& rng);
    void mutate(Chromosome& chrom, int num_nodes, std::mt19937& rng);
    Chromosome tournament_selection(const std::vector<Chromosome>& population, const std::vector<float>& fitnesses, int tournament_size, std::mt19937& rng);

public:
    QuestPlacement(
        const std::vector<Knot>& knots,
        const std::unordered_map<Knot, std::vector<Knot>, std::hash<Knot>>& adj_list,
        const std::vector<std::string>& tags,
        const std::map<std::string, int>& desired_dists,
        const std::vector<std::tuple<std::string, std::string, int>>& min_dist_constraints,
        int pop_size = 100,
        int gens = 50,
        float mut_rate = 0.01f
    );
    std::map<std::string, int> run_genetic_algorithm(int spawn_node_id, std::mt19937& rng);
};