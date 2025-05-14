#include "QuestPlacement.h"
#include <queue>
#include <algorithm>
#include <limits>
#include <iostream>
#include <unordered_map>

QuestPlacement::QuestPlacement(
    const std::vector<Knot>& knots,
    const std::unordered_map<Knot, std::vector<Knot>, std::hash<Knot>>& adj_list,
    const std::vector<std::string>& tags,
    const std::map<std::string, int>& desired_dists,
    const std::vector<std::tuple<std::string, std::string, int>>& min_dist_constraints,
    int spawn_node_id,
    int pop_size,
    int gens,
    float mut_rate
) : knots(knots), adjacency_list(adj_list), tag_list(tags), desired_distances(desired_dists),
min_distance_constraints(min_dist_constraints), spawn_node_id(spawn_node_id),
population_size(pop_size), generations(gens), mutation_rate(mut_rate) {
}

std::vector<Chromosome> QuestPlacement::initialize_population(int num_nodes, int num_tags, std::mt19937& rng) {
    std::vector<Chromosome> population;
    int spawn_tag_idx = -1;
    for (size_t i = 0; i < tag_list.size(); ++i) {
        if (tag_list[i] == "spawn") {
            spawn_tag_idx = i;
            break;
        }
    }

    for (int i = 0; i < population_size; ++i) {
        Chromosome chrom;
        chrom.node_ids.resize(num_tags);
        std::vector<int> available_nodes;
        for (int j = 0; j < num_nodes; ++j) if (j != spawn_node_id) available_nodes.push_back(j);
        std::shuffle(available_nodes.begin(), available_nodes.end(), rng);

        // Фиксируем spawn на spawn_node_id
        if (spawn_tag_idx != -1) {
            chrom.node_ids[spawn_tag_idx] = spawn_node_id;
        }

        // Назначаем остальные теги
        int avail_idx = 0;
        for (int j = 0; j < num_tags; ++j) {
            if (j != spawn_tag_idx) {
                chrom.node_ids[j] = available_nodes[avail_idx++];
            }
        }
        population.push_back(chrom);
    }
    return population;
}

float QuestPlacement::compute_fitness(const Chromosome& chrom, const std::vector<int>& distances,
    const std::vector<std::vector<int>>& dist_matrix) {
    // Проверка уникальности назначений
    std::set<int> assigned_nodes;
    for (int node_id : chrom.node_ids) {
        if (assigned_nodes.count(node_id)) return -1e9;
        assigned_nodes.insert(node_id);
    }

    // Проверка, что spawn назначен на spawn_node_id
    int spawn_tag_idx = -1;
    for (size_t i = 0; i < tag_list.size(); ++i) {
        if (tag_list[i] == "spawn") {
            spawn_tag_idx = i;
            break;
        }
    }
    if (spawn_tag_idx != -1 && chrom.node_ids[spawn_tag_idx] != spawn_node_id) {
        return -1e9; // Штраф за неправильное назначение spawn
    }

    float fitness = 0.0f;
    // Оценка соответствия желаемым дистанциям от spawn_knot
    for (size_t i = 0; i < tag_list.size(); ++i) {
        std::string tag = tag_list[i];
        int node_id = chrom.node_ids[i];
        int actual_distance = distances[node_id];
        int desired = desired_distances.at(tag);
        fitness -= std::abs(actual_distance - desired);
    }

    // Штраф за нарушение минимальных расстояний между тегами
    for (const auto& [tag1, tag2, min_dist] : min_distance_constraints) {
        auto it1 = std::find(tag_list.begin(), tag_list.end(), tag1);
        auto it2 = std::find(tag_list.begin(), tag_list.end(), tag2);
        if (it1 == tag_list.end() || it2 == tag_list.end()) continue;
        int idx1 = it1 - tag_list.begin();
        int idx2 = it2 - tag_list.begin();
        int node1 = chrom.node_ids[idx1];
        int node2 = chrom.node_ids[idx2];
        int dist = dist_matrix[node1][node2];
        if (dist < min_dist) fitness -= (min_dist - dist) * 10.0f;
    }

    return fitness;
}

Chromosome QuestPlacement::crossover(const Chromosome& parent1, const Chromosome& parent2, int num_nodes, std::mt19937& rng) {
    Chromosome child;
    child.node_ids.resize(parent1.node_ids.size());
    std::set<int> assigned;
    std::uniform_int_distribution<> dist(0, 1);
    int spawn_tag_idx = -1;
    for (size_t i = 0; i < tag_list.size(); ++i) {
        if (tag_list[i] == "spawn") {
            spawn_tag_idx = i;
            break;
        }
    }

    // Фиксируем spawn на spawn_node_id
    if (spawn_tag_idx != -1) {
        child.node_ids[spawn_tag_idx] = spawn_node_id;
        assigned.insert(spawn_node_id);
    }

    // Скрещивание остальных тегов
    for (size_t i = 0; i < parent1.node_ids.size(); ++i) {
        if (i == spawn_tag_idx) continue;
        int node = (dist(rng) == 0) ? parent1.node_ids[i] : parent2.node_ids[i];
        if (assigned.count(node)) {
            std::vector<int> available;
            for (int j = 0; j < num_nodes; ++j) if (!assigned.count(j)) available.push_back(j);
            if (!available.empty()) node = available[rng() % available.size()];
            else node = parent1.node_ids[i]; // Fallback
        }
        child.node_ids[i] = node;
        assigned.insert(node);
    }
    return child;
}

void QuestPlacement::mutate(Chromosome& chrom, int num_nodes, std::mt19937& rng) {
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    int spawn_tag_idx = -1;
    for (size_t i = 0; i < tag_list.size(); ++i) {
        if (tag_list[i] == "spawn") {
            spawn_tag_idx = i;
            break;
        }
    }

    if (dist(rng) < mutation_rate * (chrom.node_ids.size() - 1)) {
        int i = rng() % chrom.node_ids.size();
        while (i == spawn_tag_idx) i = rng() % chrom.node_ids.size();
        int j = rng() % chrom.node_ids.size();
        while (j == spawn_tag_idx || j == i) j = rng() % chrom.node_ids.size();
        std::swap(chrom.node_ids[i], chrom.node_ids[j]);
    }
}

Chromosome QuestPlacement::tournament_selection(const std::vector<Chromosome>& population,
    const std::vector<float>& fitnesses,
    int tournament_size, std::mt19937& rng) {
    std::uniform_int_distribution<> dist(0, population.size() - 1);
    float best_fitness = -std::numeric_limits<float>::max();
    Chromosome best;
    for (int i = 0; i < tournament_size; ++i) {
        int idx = dist(rng);
        if (fitnesses[idx] > best_fitness) {
            best_fitness = fitnesses[idx];
            best = population[idx];
        }
    }
    return best;
}

std::vector<int> compute_distances(const std::unordered_map<Knot, std::vector<Knot>, std::hash<Knot>>& adj_list,
    const std::vector<Knot>& knots, int start_id) {
    std::vector<int> distances(knots.size(), -1);
    std::queue<int> q;
    distances[start_id] = 0;
    q.push(start_id);
    while (!q.empty()) {
        int curr = q.front();
        q.pop();
        const Knot& current_knot = knots[curr];
        auto it = adj_list.find(current_knot);
        if (it == adj_list.end()) continue;
        for (const auto& nb : it->second) {
            int nb_id = -1;
            for (size_t i = 0; i < knots.size(); ++i) {
                if (nb.x == knots[i].x && nb.y == knots[i].y) {
                    nb_id = i;
                    break;
                }
            }
            if (nb_id != -1 && distances[nb_id] == -1) {
                distances[nb_id] = distances[curr] + 1;
                q.push(nb_id);
            }
        }
    }
    return distances;
}

std::vector<std::vector<int>> floyd_warshall(const std::unordered_map<Knot, std::vector<Knot>, std::hash<Knot>>& adj_list,
    const std::vector<Knot>& knots) {
    int n = knots.size();
    std::vector<std::vector<int>> dist(n, std::vector<int>(n, std::numeric_limits<int>::max()));
    for (int i = 0; i < n; ++i) dist[i][i] = 0;

    for (const auto& [knot, neighbors] : adj_list) {
        int u = -1;
        for (size_t i = 0; i < knots.size(); ++i) {
            if (knot.x == knots[i].x && knot.y == knots[i].y) {
                u = i;
                break;
            }
        }
        if (u == -1) continue;
        for (const auto& nb : neighbors) {
            int v = -1;
            for (size_t i = 0; i < knots.size(); ++i) {
                if (nb.x == knots[i].x && nb.y == knots[i].y) {
                    v = i;
                    break;
                }
            }
            if (v != -1) dist[u][v] = dist[v][u] = 1;
        }
    }

    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist[i][k] != std::numeric_limits<int>::max() &&
                    dist[k][j] != std::numeric_limits<int>::max() &&
                    dist[i][j] > dist[i][k] + dist[k][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
    return dist;
}

std::map<std::string, int> QuestPlacement::run_genetic_algorithm(std::mt19937& rng) {
    std::cout << "Running GA with spawn_node_id: " << spawn_node_id << " (" << knots[spawn_node_id].x << ", " << knots[spawn_node_id].y << ")\n";
    int num_nodes = knots.size();
    int num_tags = tag_list.size();

    // Вычисление расстояний от spawn_knot
    std::vector<int> distances = compute_distances(adjacency_list, knots, spawn_node_id);

    // Вычисление матрицы расстояний
    std::vector<std::vector<int>> dist_matrix = floyd_warshall(adjacency_list, knots);

    // Инициализация популяции
    auto population = initialize_population(num_nodes, num_tags, rng);
    std::vector<float> fitnesses(population_size);

    // Основной цикл генетического алгоритма
    for (int gen = 0; gen < generations; ++gen) {
        // Вычисление фитнеса
        for (int i = 0; i < population_size; ++i) {
            fitnesses[i] = compute_fitness(population[i], distances, dist_matrix);
        }

        // Создание нового поколения
        std::vector<Chromosome> new_population;
        while (new_population.size() < population_size) {
            auto parent1 = tournament_selection(population, fitnesses, 3, rng);
            auto parent2 = tournament_selection(population, fitnesses, 3, rng);
            auto child = crossover(parent1, parent2, num_nodes, rng);
            mutate(child, num_nodes, rng);
            new_population.push_back(child);
        }
        population = new_population;
    }

    // Выбор лучшей хромосомы
    float best_fitness = -std::numeric_limits<float>::max();
    Chromosome best_chrom;
    for (const auto& chrom : population) {
        float fit = compute_fitness(chrom, distances, dist_matrix);
        if (fit > best_fitness) {
            best_fitness = fit;
            best_chrom = chrom;
        }
    }

    // Формирование результата
    std::map<std::string, int> tag_assignments;
    for (size_t i = 0; i < tag_list.size(); ++i) {
        tag_assignments[tag_list[i]] = best_chrom.node_ids[i];
    }
    return tag_assignments;
}