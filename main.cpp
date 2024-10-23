#include <utility>

#include "bits_stdc++.h"
#include "utils.cpp"

using namespace std;


class SlidingTilePuzzle {

public:
    int zero_i, zero_j, parent_i, parent_j;

    string rank;

    void set_ranking(int **board) {

        string r;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                if (board[i][j] == 0) {
                    zero_i = i;
                    zero_j = j;
                }
                r += to_string(board[i][j]) + " ";
            }

        rank = r;
    }

    int **get_board() const {

        int **board = new int *[4];
        for (int i = 0; i < 4; i++)
            board[i] = new int[4];

        stringstream ss(rank);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                string num;
                ss >> num;
                board[i][j] = stoi(num);
            }
        }

        return board;
    }

    static double get_heuristic(int **board, string &heavy_puzzle) {
        double h = 0;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (board[i][j]) {
                    if (heavy_puzzle == "OFF")
                        h += (abs(board[i][j] / 4 - i) + abs(board[i][j] % 4 - j % 4));
                    else if (heavy_puzzle == "ON")
                        h += (abs(board[i][j] / 4 - i) + abs(board[i][j] % 4 - j % 4)) * board[i][j];
                }
            }
        }
        return h;
    }

    //TODO
    // CHECK THE IDEA OF PASSING THE PARENT H TO SAVE TIME.
    vector<pair<SlidingTilePuzzle, pair<double, double>>> get_successors(int **board, double g_parent, string &heavy_puzzle) const {
        // This functions generates the children, and returns a vector including pairs of the child and pair of its heuristic and g_cost
        vector<pair<SlidingTilePuzzle, pair<double, double>>> successors;

        if (zero_i and not((zero_i - 1 == parent_i) and (zero_j == parent_j))) {
            //UP

            int cell = board[zero_i - 1][zero_j];

            board[zero_i][zero_j] = cell;
            board[zero_i - 1][zero_j] = 0;
            SlidingTilePuzzle child{};
            child.zero_i = zero_i - 1;
            child.zero_j = zero_j;
            child.parent_i = zero_i;
            child.parent_j = zero_j;
            child.set_ranking(board);

            if (heavy_puzzle == "OFF")
                successors.emplace_back(child,
                                        make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle), g_parent + 1));
            else if (heavy_puzzle == "ON")
                successors.emplace_back(child, make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle),
                                                         g_parent + cell));

            board[zero_i][zero_j] = 0;
            board[zero_i - 1][zero_j] = cell;
        }

        if (zero_i != 3 and not((zero_i + 1 == parent_i) and (zero_j == parent_j))) {
            //DOWN

            int cell = board[zero_i + 1][zero_j];

            board[zero_i][zero_j] = cell;
            board[zero_i + 1][zero_j] = 0;
            SlidingTilePuzzle child{};
            child.zero_i = zero_i + 1;
            child.zero_j = zero_j;
            child.parent_i = zero_i;
            child.parent_j = zero_j;
            child.set_ranking(board);

            if (heavy_puzzle == "OFF")
                successors.emplace_back(child,
                                        make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle), g_parent + 1));
            else if (heavy_puzzle == "ON")
                successors.emplace_back(child, make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle),
                                                         g_parent + cell));

            board[zero_i][zero_j] = 0;
            board[zero_i + 1][zero_j] = cell;
        }

        if (zero_j and not((zero_i == parent_i) and (zero_j - 1 == parent_j))) {
            //LEFT

            int cell = board[zero_i][zero_j - 1];
            board[zero_i][zero_j] = cell;
            board[zero_i][zero_j - 1] = 0;
            SlidingTilePuzzle child{};
            child.zero_i = zero_i;
            child.zero_j = zero_j - 1;
            child.parent_i = zero_i;
            child.parent_j = zero_j;
            child.set_ranking(board);

            if (heavy_puzzle == "OFF")
                successors.emplace_back(child,
                                        make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle), g_parent + 1));
            else if (heavy_puzzle == "ON")
                successors.emplace_back(child, make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle),
                                                         g_parent + cell));

            board[zero_i][zero_j] = 0;
            board[zero_i][zero_j - 1] = cell;
        }

        if (zero_j != 3 and not((zero_i == parent_i) and (zero_j + 1 == parent_j))) {
            //RIGHT

            int cell = board[zero_i][zero_j + 1];
            board[zero_i][zero_j] = cell;
            board[zero_i][zero_j + 1] = 0;
            SlidingTilePuzzle child{};
            child.zero_i = zero_i;
            child.zero_j = zero_j + 1;
            child.parent_i = zero_i;
            child.parent_j = zero_j;
            child.set_ranking(board);

            if (heavy_puzzle == "OFF")
                successors.emplace_back(child,
                                        make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle), g_parent + 1));
            else if (heavy_puzzle == "ON")
                successors.emplace_back(child, make_pair(SlidingTilePuzzle::get_heuristic(board, heavy_puzzle),
                                                         g_parent + cell));

            board[zero_i][zero_j] = 0;
            board[zero_i][zero_j + 1] = cell;
        }

        return successors;
    }

    static void print_node(int **stp) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++)
                printf("%*d", 3, stp[i][j]);
            cout << endl;
        }
    }
};

bool operator<(const SlidingTilePuzzle &A, const SlidingTilePuzzle &B) {
    return A.rank < B.rank;
}

class PancakePuzzle {
public:
    string rank;
    int parent_slicePosition;

    void set_ranking(const int *board) {
        string r;

        for (int pancake = 0; pancake < 16; pancake++)
            r += to_string(board[pancake]) + " ";

        rank = r;
    }

    int *get_board() const {
        int *board = new int[16];
        stringstream ss(rank);

        for (int i = 0; i < 16; i++) {
            string num;
            ss >> num;
            board[i] = stoi(num);
        }
        return board;
    }

    static long long int get_heuristic(int *board, string &heavy_puzzle) {
        // In case of heavy problem, the g_cost is min(bottom, top) adjacent pancakes that (bottom-top>1). otherwise it is one.
        long long int h = 0;
        for (int i = 0; i < 15; i++) {
            if (heavy_puzzle == "OFF") {
                if (abs(board[i] - board[i + 1]) > 1)
                    h += 1;
                if (board[0] != 0)
                    h += 1;
            } else if (heavy_puzzle == "ON") {
                if (abs(board[i] - board[i + 1]) > 1)
                    h += min(board[i], board[i + 1]);
                if (board[0] != 0)
                    h += board[0];
            }
        }

        return h;
    }

    vector<pair<PancakePuzzle, pair<double, double>>>
    get_successors(int *board, double g_parent, string &heavy_puzzle) const {
        // This functions generates the children, and returns a vector including pairs of the child and pair of its heuristic and g_cost
        // In case of heavy problem, the g_cost is max(most bottom, most top) flipped pancakes otherwise it is one.
        vector<pair<PancakePuzzle, pair<double, double>>> successors;

        for (int slicePosition = 1; slicePosition < 16; slicePosition++) {
            //slicePosition is starts from zero in which the handle is below the first pancake. At 0, it only flips the first pancake.

            if (slicePosition != parent_slicePosition) {
                int tmp;
                for (int i = 0; i <= floor(slicePosition / 2); i++) {
                    tmp = board[i];
                    board[i] = board[slicePosition - i];
                    board[slicePosition - i] = tmp;
                }

                PancakePuzzle child{};
                child.set_ranking(board);
                child.parent_slicePosition = slicePosition;

                if (heavy_puzzle == "OFF")
                    successors.emplace_back(child,
                                            make_pair(PancakePuzzle::get_heuristic(board, heavy_puzzle), g_parent + 1));
                else if (heavy_puzzle == "ON")
                    successors.emplace_back(child, make_pair(PancakePuzzle::get_heuristic(board, heavy_puzzle),
                                                             g_parent + max(board[0], board[slicePosition])));


                for (int i = 0; i <= floor(slicePosition / 2); i++) {
                    tmp = board[i];
                    board[i] = board[slicePosition - i];
                    board[slicePosition - i] = tmp;
                }
            }
        }

        return successors;
    }

    static void print_node(int *board) {
        for (int i = 0; i < 16; i++)
            cout << board[i] << " ";
        cout << endl;
    }
};

bool operator<(const PancakePuzzle &A, const PancakePuzzle &B) {
    return A.rank < B.rank;
}

double phi(double h, double g, double w, const string &priority_function) {
    if (priority_function == "XDP")
        return double(1 / (2 * w) * (g + (2 * w - 1) * h + sqrt(pow(g - h, 2) + 4 * w * g * h)));
    else if (priority_function == "XUP")
        return double(1 / (2 * w) * (g + h + sqrt(pow(g + h, 2) + 4 * w * (w - 1) * pow(h, 2))));
    else if (priority_function == "WA*")
        return double(g + w * h);
    else
        return 0;
}

bool improved_termination_rule(double f_prime_max, double path_cost, double W, const string &improved_termination) {
    if (improved_termination == "OFF")
        return true;
    else if (improved_termination == "ON") {
        if (path_cost > W * f_prime_max)
            return true;
        else
            return false;
    }
    return true;
}

template<typename T>
long long int IOS(T startNode, double W, const string &priority_function, string heavy_puzzle, const string &improved_termination) {

    // INIT THE NUMBER OF EXPANDED AND GENERATED NODES
    long long int generated = 0, expanded = 0;

    // INIT THE W_f -> W of the Focal List
    double W_f = 2 * W - 1;

    // INIT THE PATH COST FOR THE TERMINATION
    double path_cost = INT32_MAX;

    // INIT THE QUEUES -> pair of f, h, g (the lower h = higher g is preferred) and the object.
    set<pair<double, pair<double, pair<double, T>>>> OpenList, FocalList;

    // INIT THE VISITED (CLOSED LIST)
    map<string, pair<double, string> > visited_Focal, visited_Open;

    // INIT THE PATH I -> pair of ranking of the node to its cost and ranking of its parent.
    map<string, pair<double, string> > path_I;

    // INIT THE INIT VALUES.
    double start_heuristic = T::get_heuristic(startNode.get_board(), heavy_puzzle);
    double f_Open = 0 + start_heuristic;
    double f_Focal = phi(start_heuristic, 0, W_f, priority_function);
    double f_prime_max = f_Focal;

    // PUSH THE START NODE TO THE QUEUES.
    FocalList.insert(make_pair(f_Focal, make_pair(start_heuristic, make_pair(0, startNode))));
    OpenList.insert(make_pair(f_Open, make_pair(start_heuristic, make_pair(0, startNode))));

    visited_Focal[startNode.rank] = make_pair(0, "-1");
    visited_Open[startNode.rank] = make_pair(0, "-1");

    string goal_ranking;

    clock_t start_time, end_time;
    start_time = clock();
    // THE MAIN WHILE
    while ((path_cost > W * (OpenList.begin()->first)) and
           improved_termination_rule(f_prime_max, path_cost, W_f, improved_termination)) {

        if (path_cost > FocalList.begin()->first and path_cost == INT32_MAX) {
            // EXPAND THE BEST NODE ON FOCAL AND APPEND THE CHILDREN TO IT.
            T topNode = FocalList.begin()->second.second.second;
            expanded += 1;


            double g_parent = FocalList.begin()->second.second.first;
            double h_parent = FocalList.begin()->second.first;

            // UPDATE THE F_PRIME_MAX
            f_prime_max = max(f_prime_max, g_parent / W_f + h_parent);


            // GENERATE THE SUCCESSORS OF THE TOP NODE -> pair of Node, its h and its g
            vector<pair<T, pair<double, double>>> successors = topNode.get_successors(topNode.get_board(), g_parent,
                                                                                      heavy_puzzle);


            // POP THE TOP NODE
            auto it1 = FocalList.begin();
            auto it2 = FocalList.begin();
            it2++;
            FocalList.erase(it1, it2);

            for (const auto &child: successors) {
                generated += 1;

                // CHECK IF THIS NODE IS VISITED BEFORE
                if (visited_Focal.find(child.first.rank) == visited_Focal.end()) {
                    // ADD THE CHILD TO FOCAL ONCE-SEEN
                    double g_child = child.second.second;
                    double h_child = child.second.first;
                    visited_Focal[child.first.rank] = make_pair(g_child, topNode.rank);


                    // GOAL CHECK | IF THE CHILD NODE IS THE GOAL, RE-DEFINE THE PATH-COST.
                    if (h_child == 0) // heuristic of the goal is zero
                    {
                        path_cost = g_child;
                        goal_ranking = child.first.rank;

                        // STORING THE PATH
                        string tmp_rank = goal_ranking;
                        while (tmp_rank != "-1") {
                            path_I[tmp_rank] = make_pair(visited_Focal[tmp_rank].first, visited_Focal[tmp_rank].second);
                            tmp_rank = visited_Focal[tmp_rank].second;
                        }
                        break;
                    }

                    // ADD THE CHILD TO THE FOCAL LIST
                    double f_child_Focal = phi(h_child, g_child, W_f, priority_function);
                    FocalList.insert(make_pair(f_child_Focal, make_pair(h_child, make_pair(g_child, child.first))));

                    // UPDATE THE F_PRIME_MAX
//                    f_prime_max = max(f_prime_max, g_child / W_f + h_child);
                }
            }

        } else {
            // NOW WE HAVE FOUND THE SOLUTION, WE WANT TO PROVE IT IS W-OPTIMAL.
            // WE HAVE TO EXPAND FROM THE OPEN LIST AND APPEND THE CHILDREN TO IT.
            T topNode = OpenList.begin()->second.second.second;
            expanded += 1;

            double g_parent = OpenList.begin()->second.second.first;

            // GENERATE THE SUCCESSORS OF THE TOP NODE
            vector<pair<T, pair<double, double>>> successors = topNode.get_successors(topNode.get_board(), g_parent, heavy_puzzle);

            // POP THE TOP NODE
            auto it1 = OpenList.begin();
            auto it2 = OpenList.begin();
            it2++;
            OpenList.erase(it1, it2);

            for (const auto &child: successors) {
                generated += 1;

                if (visited_Open.find(child.first.rank) == visited_Open.end()) {
                    // IF THE G_COST OF ONE OF CHILDREN WAS LESS THAN THE PREVIOUSLY SEEN COST FOR THAT CHILD, UPDATE IT.
                    double g_child = child.second.second;
                    double h_child = child.second.first;

                    if ((path_I.find(child.first.rank) != path_I.end()) and path_I[child.first.rank].first > g_child) {
                        // UPDATE THE PATH COST
                        double cost_reduction = path_I[child.first.rank].first - g_child;
                        path_cost -= cost_reduction;

                        // RECREATE THE PATH -> start from goal, go back to root and find the nodes on the path.
                        map<string, pair<double, string> > copy_path_I;
                        string tmp_rank = goal_ranking;
                        while (tmp_rank != child.first.rank) {
                            copy_path_I[tmp_rank] = make_pair(path_I[tmp_rank].first - cost_reduction,
                                                              path_I[tmp_rank].second);
                            tmp_rank = path_I[tmp_rank].second;
                        }

                        path_I.clear();

                        copy_path_I[child.first.rank] = make_pair(g_child, topNode.rank);
                        tmp_rank = topNode.rank;

                        while (tmp_rank != "-1") {
                            copy_path_I[tmp_rank] = make_pair(visited_Open[tmp_rank].first,
                                                              visited_Open[tmp_rank].second);
                            tmp_rank = visited_Open[tmp_rank].second;
                        }

                        tmp_rank = goal_ranking;
                        while (tmp_rank != "-1") {
                            path_I[tmp_rank] = make_pair(copy_path_I[tmp_rank].first, copy_path_I[tmp_rank].second);
                            tmp_rank = copy_path_I[tmp_rank].second;
                        }
                        copy_path_I.clear();
                    }

                    // ADD THE CHILD TO OPEN ONCE-SEEN
                    visited_Open[child.first.rank] = make_pair(g_child, topNode.rank);

                    // ADD THE CHILD TO THE OPEN LIST
                    double f_child_Open = g_child + h_child;
                    OpenList.insert(make_pair(f_child_Open, make_pair(h_child, make_pair(g_child, child.first))));

                }
            }
        }
    }

    // PRINT THE PATH
//    int state_number = 0;
//    cout << "The Solution is :" << endl;
//    while (goal_ranking != "-1") {
//        T node{};
//        node.rank = goal_ranking;
//        cout << "state " << state_number << ":" << endl;
//        T::print_node(node.get_board());
//
//        cout << endl << "====================" << endl << endl;
//        goal_ranking = path_I[goal_ranking].second;
//        state_number += 1;
//    }


    cout << "Found the solution!..." << endl;
////    cout << state_number - 1 << " is the solution length." << endl;
//    cout << path_cost << " is the path cost." << endl;
//    cout << expanded << " number of nodes expanded." << endl;
//    cout << generated << " number of nodes generated." << endl;
//    end_time = clock();
//    cout << "Total Time was : " << (float) (end_time - start_time) / (float) CLOCKS_PER_SEC << " seconds" << endl;
//    cout << "=========" << endl;

    return expanded;

}

void study_of_Improved_Termination() {
//    cout << endl << endl;
//    cout << "===========================================" << endl;
//    cout << "====== Study of Improved Termination ======" << endl;
//    cout << "===========================================" << endl;
//    cout << endl;
//    cout << "Domain: Heavy Sliding Tile Puzzle." << endl;
//
//
//    int cnt = 3;
//    double bound = 3;
//    double average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    double average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  With condition: "<<average_node_expansion2<<endl;
//    cout << "  No condition:  "<<average_node_expansion<<endl;
//    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
//    cout<<"=========================="<<endl;
//
//    cnt = 3;
//    bound = 2;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  With condition: "<<average_node_expansion2<<endl;
//    cout << "  No condition:  "<<average_node_expansion<<endl;
//    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
//    cout<<"=========================="<<endl;
//
//    cnt = 3;
//    bound = 1.5;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  With condition: "<<average_node_expansion2<<endl;
//    cout << "  No condition:  "<<average_node_expansion<<endl;
//    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
//    cout<<"=========================="<<endl;
//
//    cnt = 3;
//    bound = 1.25;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  With condition: "<<average_node_expansion2<<endl;
//    cout << "  No condition:  "<<average_node_expansion<<endl;
//    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//    cout<<"=========================="<<endl;
//    cout<<"=========================="<<endl;
    cout << endl;
    cout << "Domain: HEAVY Pancake Puzzle." << endl;


    int cnt = 0;
    double bound = 3;
    double average_node_expansion = 0;
//    while (cnt < 51) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
//    }
//    average_node_expansion /= 50;
//
//    cnt = 0;
    double average_node_expansion2 = 0;
//    while (cnt < 51) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 50;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  With condition: "<<average_node_expansion2<<endl;
//    cout << "  No condition:  "<<average_node_expansion<<endl;
//    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
//    cout<<"=========================="<<endl;

    cnt = 0;
    bound = 2;
    average_node_expansion = 0;
    while (cnt < 51) {
        int *startBoard = pick_next_pancake(cnt);
        PancakePuzzle startNode{};
        startNode.set_ranking(startBoard);
        startNode.parent_slicePosition = -1;

        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
    }
    average_node_expansion /= 50;

    cnt = 0;
    average_node_expansion2 = 0;
    while (cnt < 51) {
        int *startBoard = pick_next_pancake(cnt);
        PancakePuzzle startNode{};
        startNode.set_ranking(startBoard);
        startNode.parent_slicePosition = -1;

        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
    }
    average_node_expansion2 /= 50;

    cout << "  Bound: "<<bound<<endl;
    cout << "  With condition: "<<average_node_expansion2<<endl;
    cout << "  No condition:  "<<average_node_expansion<<endl;
    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
    cout<<"=========================="<<endl;

//    cnt = 0;
//    bound = 1.5;
//    average_node_expansion = 0;
//    while (cnt < 51) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
//    }
//    average_node_expansion /= 50;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 51) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 50;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  With condition: "<<average_node_expansion2<<endl;
//    cout << "  No condition:  "<<average_node_expansion<<endl;
//    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
//    cout<<"=========================="<<endl;
//
//    cnt = 0;
//    bound = 1.25;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "WA*", "ON", "OFF");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  With condition: "<<average_node_expansion2<<endl;
//    cout << "  No condition:  "<<average_node_expansion<<endl;
//    cout << "  Gain:  "<< (average_node_expansion - average_node_expansion2)/average_node_expansion*100<<" %"<<endl;
//    cout<<"=========================="<<endl;


}

void study_of_Alternate_Priority_Functions() {

    cout << endl << endl;
    cout << "===========================================" << endl;
    cout << "== Study of Alternate Priority Functions ==" << endl;
    cout << "===========================================" << endl;
//    cout << endl;
//    cout << "Domain: Sliding Tile Puzzle." << endl;
//    cout <<endl<<endl;
//
//    int cnt = 3;
//    double bound = 3;
//    double average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    double average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
//    double average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 3;
//    bound = 2;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
//    average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 3;
//    bound = 1.5;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
//    average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 3;
//    bound = 1.25;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
//    average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout <<"==========================================="<<endl;
//    cout <<"===========================================" << endl;
//    cout << endl;
//    cout << "Domain: Heavy Sliding Tile Puzzle." << endl;
//    cout <<endl<<endl;

    int cnt = 3;
    double bound = 3;
    double average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
    double average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
    double average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;

//    cnt = 3;
//    bound = 2;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
//    average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 3;
//    bound = 1.5;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
//    average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 3;
//    bound = 1.25;
//    average_node_expansion = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 3;
//    average_node_expansion2 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 3;
//    average_node_expansion3 = 0;
//    while (cnt < 103) {
//        int **startBoard = pick_next_stp(cnt);
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;

//    cout <<"==========================================="<<endl;
//    cout <<"===========================================" << endl;
//    cout << endl;
    cout << "Domain: Heavy Pancake Puzzle." << endl;
    cout <<endl<<endl;

    cnt = 0;
    bound = 3;
    average_node_expansion = 0;
    while (cnt < 50) {
        int *startBoard = pick_next_pancake(cnt);
        PancakePuzzle startNode{};
        startNode.set_ranking(startBoard);
        startNode.parent_slicePosition = -1;

        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
        cout<<cnt<<endl;
    }
    average_node_expansion /= 50;
    cout<<average_node_expansion<<endl;

//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 50) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 50;
//    cout<<average_node_expansion2<<endl;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 50) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 50;
//    cout<<average_node_expansion3<<endl;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;

    cout << endl;
//
//    cnt = 0;
//    bound = 2;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 0;
//    bound = 1.5;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 0;
//    bound = 1.25;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "ON", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "ON", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "ON", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout <<"==========================================="<<endl;
//    cout <<"===========================================" << endl;
//    cout << endl;
//    cout << "Domain: Pancake Puzzle." << endl;
//    cout <<endl<<endl;
//
//    cnt = 0;
//    bound = 3;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 0;
//    bound = 2;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 0;
//    bound = 1.5;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;
//
//    cout << endl;
//
//    cnt = 0;
//    bound = 1.25;
//    average_node_expansion = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion += IOS(startNode, bound, "XDP", "OFF", "ON");
//    }
//    average_node_expansion /= 100;
//
//    cnt = 0;
//    average_node_expansion2 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion2 += IOS(startNode, bound, "WA*", "OFF", "ON");
//    }
//    average_node_expansion2 /= 100;
//
//    cnt = 0;
//    average_node_expansion3 = 0;
//    while (cnt < 1001) {
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//
//        average_node_expansion3 += IOS(startNode, bound, "XUP", "OFF", "ON");
//    }
//    average_node_expansion3 /= 100;
//
//    cout << "  Bound: "<<bound<<endl;
//    cout << "  phi=XDP:  "<<average_node_expansion<<endl;
//    cout << "  phi=WA*:  "<<average_node_expansion2<<endl;
//    cout << "  phi=XUP:  "<<average_node_expansion3<<endl;
//    cout<<"=========================="<<endl;

}

int main() {

//    study_of_Improved_Termination();
    study_of_Alternate_Priority_Functions();

//    cout << "Testing STP with WA*..." << endl;
//    int cnt = 3;
//    while (cnt<103)
//    {
//        cout<<"Problem "<<cnt-2<<endl;
//        int **startBoard = pick_next_stp(cnt);
//        // DEBUG
//        {
//            SlidingTilePuzzle::print_node(startBoard);
////            int a;
////            cin >> a;
//        }
//        SlidingTilePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_i = startNode.zero_i;
//        startNode.parent_j = startNode.zero_j;
//
//        IOS(startNode, 1.25, "XDP", "OFF", "OFF");
//
//    }

//    cout << "Testing Pancakes with WA*..." << endl;
//    int cnt = 0;
//    double average_node_expansion =0;
//    while (cnt < 50) {
//        cout << "Problem " << cnt + 1 << endl;
//        int *startBoard = pick_next_pancake(cnt);
//        PancakePuzzle startNode{};
//        startNode.set_ranking(startBoard);
//        startNode.parent_slicePosition = -1;
//        // DEBUG
//        {
//            PancakePuzzle::print_node(startBoard);
////            int a;
////            cin >> a;
//        }
//
//        average_node_expansion += IOS(startNode, 3, "WA*", "ON", "ON");
//
//    }
//    average_node_expansion /=50;
//    cout<<average_node_expansion<<endl;


    return 0;
}
