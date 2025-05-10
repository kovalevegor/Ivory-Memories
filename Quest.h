#pragma once

#include "Knot.h"
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <string>
#include <functional>
#include <random>
#include <unordered_map>
#include <optional>

// Специальный узел для игрока
const Knot WITH_PLAYER(-1.0f, -1.0f);
const int PLAYER_KNOT_ID = -1;

// Типы источников предметов
enum SourceType { NODE_SOURCE, CHARACTER_SOURCE, ENEMY_SOURCE };

struct Source {
    SourceType type;
    int id; // knot ID для NODE_SOURCE, char ID для CHARACTER_SOURCE, enemy ID для ENEMY_SOURCE

    bool operator==(const Source& other) const {
        return type == other.type && id == other.id;
    }
};

// Структура игрового мира
struct World {
    std::vector<Knot> knots;
    std::map<int, std::vector<int>> adjacency_list; // knot ID к списку соседних knot IDs
    std::map<int, int> character_locations; // char ID к knot ID
    std::map<int, int> enemy_locations; // enemy ID к knot ID
    std::map<int, std::vector<Source>> item_sources; // item ID к списку источников
    std::map<int, int> character_requirements; // char ID к требуемому item ID или -1
};

// Типы действий
enum ActionType { MOVE, PICK_UP, TALK_TO, DEFEAT };

// Структура действия
struct Action {
    ActionType type;
    int target; // knot ID для MOVE
    int id; // для PICK_UP, TALK_TO, DEFEAT
};

// Структура состояния
struct State {
    int player_location; // knot ID
    std::map<int, int> item_locations; // item ID к knot ID или -1 для игрока
    std::set<int> talked_to_characters;
    std::set<int> defeated_enemies;

    bool operator==(const State& other) const {
        return player_location == other.player_location &&
            item_locations == other.item_locations &&
            talked_to_characters == other.talked_to_characters &&
            defeated_enemies == other.defeated_enemies;
    }
};

// Хэш для State
namespace std {
    template<>
    struct hash<State> {
        size_t operator()(const State& s) const {
            size_t h = hash<int>()(s.player_location);
            for (const auto& [item, loc] : s.item_locations) {
                h ^= hash<int>()(item) ^ hash<int>()(loc);
            }
            for (int char_id : s.talked_to_characters) {
                h ^= hash<int>()(char_id);
            }
            for (int enemy_id : s.defeated_enemies) {
                h ^= hash<int>()(enemy_id);
            }
            return h;
        }
    };
}

// Структура цели
struct Goal {
    std::optional<int> player_location; // knot ID
    std::set<int> required_items; // предметы у игрока
    std::set<int> required_talked_to;
    std::set<int> required_defeated;

    bool operator==(const Goal& other) const {
        return player_location == other.player_location &&
            required_items == other.required_items &&
            required_talked_to == other.required_talked_to &&
            required_defeated == other.required_defeated;
    }
};

// Хэш для Goal
namespace std {
    template<>
    struct hash<Goal> {
        size_t operator()(const Goal& g) const {
            size_t h = g.player_location ? hash<int>()(*g.player_location) : 0;
            for (int item : g.required_items) {
                h ^= hash<int>()(item);
            }
            for (int char_id : g.required_talked_to) {
                h ^= hash<int>()(char_id);
            }
            for (int enemy_id : g.required_defeated) {
                h ^= hash<int>()(enemy_id);
            }
            return h;
        }
    };
}

// Хэш для std::pair<State, Goal>
namespace std {
    template<>
    struct hash<std::pair<State, Goal>> {
        size_t operator()(const std::pair<State, Goal>& p) const {
            return hash<State>()(p.first) ^ hash<Goal>()(p.second);
        }
    };
}

// Структура узла квестового дерева
struct QuestNode {
    std::string goal_description;
    std::vector<Action> plan;
    std::vector<QuestNode*> children;
    ~QuestNode() {
        for (auto child : children) {
            delete child;
        }
    }
};

// Индивидуум для генетического алгоритма
struct Individual {
    std::map<int, std::vector<Source>> selected_sources; // item ID к подмножеству источников
    double fitness;
};

// Объявления функций
State apply_action(const State& s, const Action& a);
bool is_applicable(const Action& a, const State& s, const World& w);
std::vector<Action> get_applicable_actions(const State& s, const World& w);
bool satisfies_goal(const State& s, const Goal& g);
std::vector<Action> find_plan(const State& initial, const Goal& goal, const World& w, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache);
QuestNode* build_obtain_item_tree(int item, const std::vector<Source>& available_sources, const World& w, const State& initial, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache);
QuestNode* build_fetch_quest_tree(int item, int quest_giver_knot, const Individual& ind, const World& w, const State& initial, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache);
void save_quest_tree(const QuestNode* node, const std::string& filename);
std::vector<Individual> initialize_population(int population_size, const World& w, int num_items);
double evaluate_fitness(Individual& ind, const World& w, int item, int quest_giver_knot, const State& initial, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache);
void crossover(const Individual& parent1, const Individual& parent2, Individual& child, std::mt19937& rng);
void mutate(Individual& ind, const World& w, std::mt19937& rng);
Individual genetic_algorithm(int population_size, int generations, const World& w, int item, int quest_giver_knot, const State& initial, std::mt19937& rng);

