#include "Quest.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>
#include <unordered_set>

State apply_action(const State& s, const Action& a) {
    State new_s = s;
    switch (a.type) {
    case MOVE:
        new_s.player_location = a.target;
        break;
    case PICK_UP:
        new_s.item_locations[a.id] = PLAYER_KNOT_ID;
        break;
    case TALK_TO:
        new_s.talked_to_characters.insert(a.id);
        break;
    case DEFEAT:
        new_s.defeated_enemies.insert(a.id);
        break;
    }
    return new_s;
}

bool is_applicable(const Action& a, const State& s, const World& w) {
    switch (a.type) {
    case MOVE: {
        auto it = w.adjacency_list.find(s.player_location);
        if (it == w.adjacency_list.end()) return false;
        return std::find(it->second.begin(), it->second.end(), a.target) != it->second.end();
    }
    case PICK_UP: {
        auto it = s.item_locations.find(a.id);
        if (it == s.item_locations.end()) return false;
        return it->second == s.player_location;
    }
    case TALK_TO: {
        auto char_it = w.character_locations.find(a.id);
        if (char_it == w.character_locations.end()) return false;
        return char_it->second == s.player_location && s.talked_to_characters.find(a.id) == s.talked_to_characters.end();
    }
    case DEFEAT: {
        auto enemy_it = w.enemy_locations.find(a.id);
        if (enemy_it == w.enemy_locations.end()) return false;
        return enemy_it->second == s.player_location && s.defeated_enemies.find(a.id) == s.defeated_enemies.end();
    }
    }
    return false;
}

std::vector<Action> get_applicable_actions(const State& s, const World& w) {
    std::vector<Action> actions;
    auto neighbors = w.adjacency_list.at(s.player_location);
    for (int neighbor : neighbors) {
        actions.push_back({ MOVE, neighbor, -1 });
    }
    for (const auto& [item, loc] : s.item_locations) {
        if (loc == s.player_location) {
            actions.push_back({ PICK_UP, -1, item });
        }
    }
    for (const auto& [char_id, char_loc] : w.character_locations) {
        if (char_loc == s.player_location && s.talked_to_characters.find(char_id) == s.talked_to_characters.end()) {
            actions.push_back({ TALK_TO, -1, char_id });
        }
    }
    for (const auto& [enemy_id, enemy_loc] : w.enemy_locations) {
        if (enemy_loc == s.player_location && s.defeated_enemies.find(enemy_id) == s.defeated_enemies.end()) {
            actions.push_back({ DEFEAT, -1, enemy_id });
        }
    }
    return actions;
}

bool satisfies_goal(const State& s, const Goal& g) {
    if (g.player_location && *g.player_location != s.player_location) return false;
    for (int item : g.required_items) {
        auto it = s.item_locations.find(item);
        if (it == s.item_locations.end() || it->second != PLAYER_KNOT_ID) return false;
    }
    for (int char_id : g.required_talked_to) {
        if (s.talked_to_characters.find(char_id) == s.talked_to_characters.end()) return false;
    }
    for (int enemy_id : g.required_defeated) {
        if (s.defeated_enemies.find(enemy_id) == s.defeated_enemies.end()) return false;
    }
    return true;
}

std::vector<Action> find_plan(const State& initial, const Goal& goal, const World& w, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache) {
    // Проверяем кэш
    auto cache_key = std::make_pair(initial, goal);
    auto cache_it = plan_cache.find(cache_key);
    if (cache_it != plan_cache.end()) {
        return cache_it->second;
    }

    std::queue<std::pair<State, std::vector<Action>>> q;
    std::unordered_set<State> visited;
    q.push({ initial, {} });
    visited.insert(initial);
    while (!q.empty()) {
        auto [current, path] = q.front();
        q.pop();
        if (satisfies_goal(current, goal)) {
            plan_cache[cache_key] = path; // Сохраняем в кэш
            return path;
        }
        for (const auto& action : get_applicable_actions(current, w)) {
            State new_s = apply_action(current, action);
            if (visited.find(new_s) == visited.end()) {
                visited.insert(new_s);
                auto new_path = path;
                new_path.push_back(action);
                q.push({ new_s, new_path });
            }
        }
    }
    plan_cache[cache_key] = {}; // Кэшируем пустой план
    return {};
}

QuestNode* build_obtain_item_tree(int item, const std::vector<Source>& available_sources, const World& w, const State& initial, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache) {
    QuestNode* node = new QuestNode{ "Obtain item " + std::to_string(item), {}, {} };
    for (const auto& source : available_sources) {
        QuestNode* child = new QuestNode;
        if (source.type == NODE_SOURCE) {
            int knot_id = source.id;
            child->goal_description = "Pick up item " + std::to_string(item) + " from knot " + std::to_string(knot_id);
            Goal g;
            g.required_items = { item };
            State init = initial;
            init.item_locations[item] = knot_id;
            child->plan = find_plan(init, g, w, plan_cache);
        }
        else if (source.type == CHARACTER_SOURCE) {
            int char_id = source.id;
            int req_item = w.character_requirements.at(char_id);
            if (req_item == -1) {
                child->goal_description = "Talk to character " + std::to_string(char_id) + " to get item " + std::to_string(item);
                Goal g;
                g.required_talked_to = { char_id };
                child->plan = find_plan(initial, g, w, plan_cache);
            }
            else {
                child->goal_description = "Obtain item " + std::to_string(req_item) + " then talk to character " + std::to_string(char_id);
                QuestNode* obtain_req = build_obtain_item_tree(req_item, w.item_sources.at(req_item), w, initial, plan_cache);
                child->children.push_back(obtain_req);
                State after_req = initial;
                for (const auto& action : obtain_req->plan) {
                    after_req = apply_action(after_req, action);
                }
                Goal g;
                g.required_talked_to = { char_id };
                child->plan = find_plan(after_req, g, w, plan_cache);
            }
        }
        else if (source.type == ENEMY_SOURCE) {
            int enemy_id = source.id;
            child->goal_description = "Defeat enemy " + std::to_string(enemy_id) + " to get item " + std::to_string(item);
            Goal g;
            g.required_defeated = { enemy_id };
            child->plan = find_plan(initial, g, w, plan_cache);
        }
        if (!child->plan.empty() || !child->children.empty()) {
            node->children.push_back(child);
        }
        else {
            delete child;
        }
    }
    return node;
}

QuestNode* build_fetch_quest_tree(int item, int quest_giver_knot, const Individual& ind, const World& w, const State& initial, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache) {
    QuestNode* root = new QuestNode{ "Fetch item " + std::to_string(item) + " and bring to knot " + std::to_string(quest_giver_knot), {}, {} };
    QuestNode* obtain_item_node = build_obtain_item_tree(item, ind.selected_sources.at(item), w, initial, plan_cache);
    auto attach_go_to_Q = [&](QuestNode* node, auto& self) -> void {
        if (node->children.empty() && !node->plan.empty()) {
            QuestNode* go_to_Q = new QuestNode;
            go_to_Q->goal_description = "Go to knot " + std::to_string(quest_giver_knot);
            State after_obtain = initial;
            for (const auto& action : node->plan) {
                after_obtain = apply_action(after_obtain, action);
            }
            Goal g;
            g.player_location = quest_giver_knot;
            g.required_items = { item };
            go_to_Q->plan = find_plan(after_obtain, g, w, plan_cache);
            if (!go_to_Q->plan.empty()) {
                node->children.push_back(go_to_Q);
            }
            else {
                delete go_to_Q;
            }
        }
        else {
            for (auto child : node->children) {
                self(child, self);
            }
        }
        };
    attach_go_to_Q(obtain_item_node, attach_go_to_Q);
    if (!obtain_item_node->children.empty()) {
        root->children.push_back(obtain_item_node);
    }
    else {
        delete obtain_item_node;
    }
    return root;
}

void save_quest_tree(const QuestNode* node, const std::string& filename) {
    std::ofstream out(filename);
    std::function<void(const QuestNode*, int, int&)> traverse = [&](const QuestNode* n, int depth, int& node_id) {
        int current_id = node_id++;
        out << std::string(depth * 2, ' ') << "Node " << current_id << ": " << n->goal_description << "\n";
        out << std::string(depth * 2 + 2, ' ') << "Plan: ";
        for (const auto& action : n->plan) {
            if (action.type == MOVE) {
                out << "Move to knot " << action.target << "; ";
            }
            else if (action.type == PICK_UP) {
                out << "Pick up item " << action.id << "; ";
            }
            else if (action.type == TALK_TO) {
                out << "Talk to character " << action.id << "; ";
            }
            else {
                out << "Defeat enemy " << action.id << "; ";
            }
        }
        out << "\n";
        for (auto child : n->children) {
            traverse(child, depth + 1, node_id);
        }
        };
    int node_id = 0;
    traverse(node, 0, node_id);
    out.close();
}

std::vector<Individual> initialize_population(int population_size, const World& w, int num_items) {
    std::vector<Individual> population;
    std::random_device rd;
    std::mt19937 rng(rd());
    for (int i = 0; i < population_size; ++i) {
        Individual ind;
        for (int item = 0; item < num_items; ++item) {
            auto& sources = w.item_sources.at(item);
            std::vector<Source> selected;
            int num_sources = std::uniform_int_distribution<>(1, sources.size())(rng);
            std::vector<Source> shuffled = sources;
            std::shuffle(shuffled.begin(), shuffled.end(), rng);
            for (int j = 0; j < num_sources && j < shuffled.size(); ++j) {
                selected.push_back(shuffled[j]);
            }
            ind.selected_sources[item] = selected;
        }
        population.push_back(ind);
    }
    return population;
}

double evaluate_fitness(Individual& ind, const World& w, int item, int quest_giver_knot, const State& initial, std::unordered_map<std::pair<State, Goal>, std::vector<Action>>& plan_cache) {
    QuestNode* tree = build_fetch_quest_tree(item, quest_giver_knot, ind, w, initial, plan_cache);
    int num_branches = 0;
    double avg_plan_length = 0.0;
    std::vector<int> plan_lengths;
    std::function<void(QuestNode*)> count_branches = [&](QuestNode* node) {
        if (node->children.empty()) {
            int length = node->plan.size();
            for (auto child : node->children) {
                length += child->plan.size();
            }
            plan_lengths.push_back(length);
            num_branches++;
        }
        else {
            for (auto child : node->children) {
                count_branches(child);
            }
        }
        };
    count_branches(tree);
    if (num_branches > 0) {
        for (int len : plan_lengths) {
            avg_plan_length += len;
        }
        avg_plan_length /= num_branches;
    }
    delete tree;
    double fitness = 0.0;
    if (num_branches >= 2 && num_branches <= 4 && avg_plan_length >= 5 && avg_plan_length <= 15) {
        fitness = 1.0;
    }
    else {
        fitness += (num_branches >= 2 && num_branches <= 4) ? 0.4 : 0.0;
        fitness += (avg_plan_length >= 5 && avg_plan_length <= 15) ? 0.4 : 0.0;
        fitness += std::min(0.2, 0.2 * num_branches / 4.0);
    }
    ind.fitness = fitness;
    return fitness;
}

void crossover(const Individual& parent1, const Individual& parent2, Individual& child, std::mt19937& rng) {
    child.selected_sources.clear();
    for (const auto& [item, sources1] : parent1.selected_sources) {
        const auto& sources2 = parent2.selected_sources.at(item);
        std::vector<Source> child_sources;
        int size = std::uniform_int_distribution<>(1, std::max(sources1.size(), sources2.size()))(rng);
        std::vector<Source> combined = sources1;
        combined.insert(combined.end(), sources2.begin(), sources2.end());
        std::shuffle(combined.begin(), combined.end(), rng);
        for (int i = 0; i < size && i < combined.size(); ++i) {
            if (std::find(child_sources.begin(), child_sources.end(), combined[i]) == child_sources.end()) {
                child_sources.push_back(combined[i]);
            }
        }
        child.selected_sources[item] = child_sources;
    }
}

void mutate(Individual& ind, const World& w, std::mt19937& rng) {
    for (auto& [item, sources] : ind.selected_sources) {
        if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < 0.1) {
            sources.clear();
            auto& all_sources = w.item_sources.at(item);
            int num_sources = std::uniform_int_distribution<>(1, all_sources.size())(rng);
            std::vector<Source> shuffled = all_sources;
            std::shuffle(shuffled.begin(), shuffled.end(), rng);
            for (int i = 0; i < num_sources && i < shuffled.size(); ++i) {
                sources.push_back(shuffled[i]);
            }
        }
    }
}

Individual genetic_algorithm(int population_size, int generations, const World& w, int item, int quest_giver_knot, const State& initial, std::mt19937& rng) {
    auto population = initialize_population(population_size, w, w.item_sources.size());
    std::unordered_map<std::pair<State, Goal>, std::vector<Action>> plan_cache;
    for (int gen = 0; gen < generations; ++gen) {
        for (auto& ind : population) {
            evaluate_fitness(ind, w, item, quest_giver_knot, initial, plan_cache);
        }
        std::sort(population.begin(), population.end(), [](const auto& a, const auto& b) { return a.fitness > b.fitness; });
        std::vector<Individual> new_population;
        new_population.push_back(population[0]); // Элитизм
        while (new_population.size() < population_size) {
            int i1 = std::uniform_int_distribution<>(0, population_size / 2)(rng);
            int i2 = std::uniform_int_distribution<>(0, population_size / 2)(rng);
            Individual child;
            crossover(population[i1], population[i2], child, rng);
            mutate(child, w, rng);
            new_population.push_back(child);
        }
        population = new_population;
    }
    evaluate_fitness(population[0], w, item, quest_giver_knot, initial, plan_cache);
    return population[0];
}