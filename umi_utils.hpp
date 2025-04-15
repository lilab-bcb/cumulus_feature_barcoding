#ifndef UMI_UTILS
#define UMI_UTILS

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "barcode_utils.hpp"

typedef std::unordered_map<uint64_t, int> UMI2Count;
typedef std::unordered_map<uint64_t, std::vector<int>> UMITable;
typedef UMITable::iterator UMITableIter;

class DisjointSet {
    private:
        std::vector<uint64_t> bid;
        std::vector<int> count;
        std::vector<int> parent;
        std::vector<int> rank;
        std::vector<int> tree_max;
        std::vector<int> argsorted;
        std::vector<int> order;

        std::vector<int> path;
    public:
        DisjointSet() {}

        void init(UMI2Count& umi2count) {
            int n = umi2count.size();

            bid = std::vector<uint64_t>(n, 0);
            count = std::vector<int>(n, 0);
            parent = std::vector<int>(n);
            rank = std::vector<int>(n, 0);
            tree_max = std::vector<int>(n);
            argsorted = std::vector<int>(n);
            order = std::vector<int>(n);

            int i = 0;
            for (auto& kv: umi2count) {
                bid[i] = kv.first;
                count[i] = kv.second;
                parent[i] = i;
                tree_max[i] = i;
                argsorted[i] = i;
                ++i;
            }

            std::sort(argsorted.begin(), argsorted.end(),
                [this](size_t i1, size_t i2) {
                    if (count[i1] != count[i2])
                        return count[i1] > count[i2];
                    else
                        return bid[i1] < bid[i2];
            });

            for (int i = 0; i < argsorted.size(); ++i)
                order[argsorted[i]] = i;
        }

        std::vector<uint64_t>& get_bid() {
            return bid;
        }

        std::vector<int>& get_count() {
            return count;
        }

        std::vector<int>& get_parent() {
            return parent;
        }

        std::vector<int>& get_tree_max() {
            return tree_max;
        }

        std::vector<int>& get_argsorted() {
            return argsorted;
        }

        std::vector<int>& get_order() {
            return order;
        }

        int find_set(int i) { // with path compression optimization
            path.clear();
            while (i != parent[i]) {
                path.push_back(i);
                i = parent[i];
            }
            for (int& v: path)
                parent[v] = i;
            return i;
        }

        void union_sets(int a, int b) {
            a = find_set(a);
            b = find_set(b);

            if (a != b) {
                if (rank[a] <= rank[b]) {
                    parent[a] = b;
                    if (order[tree_max[a]] < order[tree_max[b]])
                        tree_max[b] = tree_max[a];
                } else if (rank[a] > rank[b]) {
                    parent[b] = a;
                    if (order[tree_max[a]] > order[tree_max[b]])
                        tree_max[a] = tree_max[b];
                }
                if (rank[a] == rank[b])
                    ++rank[b];
            }
        }

        int find_union_max(int a) {
            return tree_max[find_set(a)];
        }

        bool connected(int a, int b) {
            return find_set(a) == find_set(b);
        }

        void set_parent(int a, int p) {
            parent[a] = p;
        }

        void set_tree_max(int a, int idx_max) {
            tree_max[a] = idx_max;
        }

        UMI2Count get_roots() {
            UMI2Count res;
            UMI2Count::iterator it;

            for (int i = 0; i < bid.size(); ++i) {
                int root = find_set(i);
                it = res.find(bid[tree_max[root]]);
                if (it != res.end())
                    it->second += count[i];
                else
                    res.insert(std::make_pair(bid[tree_max[root]], count[i]));
            }
            return res;
        }

        void clear() {
            bid.clear();
            count.clear();
            parent.clear();
            rank.clear();
            tree_max.clear();
            argsorted.clear();
            order.clear();
        }
};

class UMICorrectSet {
    private:
        std::vector<UMITable> match_tables;
        DisjointSet ds;
        int umi_length;

        void correct_by_cluster() {
            std::vector<int>& argsorted = this->ds.get_argsorted();
            std::vector<uint64_t>& bids = this->ds.get_bid();

            for (const int& i: argsorted)
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = bids[i] & (~aux_arr[j][NNUC]);
                    UMITableIter it = match_tables[j].find(bid_sub);
                    if (it == match_tables[j].end())
                        match_tables[j].insert(std::make_pair(bid_sub, std::vector<int>(1, i)));
                    else {
                        this->ds.union_sets(it->second[0], i);
                        it->second.push_back(i);
                    }
                }
        }

        void correct_by_adjacency() {
            std::vector<int>& argsorted = this->ds.get_argsorted();
            std::vector<uint64_t>& bids = this->ds.get_bid();
            std::vector<int>& order = this->ds.get_order();
            std::vector<int>& parent = this->ds.get_parent();

            std::vector<int> edge(bids.size(), -1);
            for (const int& i: argsorted) {
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = bids[i] & (~aux_arr[j][NNUC]);
                    UMITableIter it = match_tables[j].find(bid_sub);
                    if (it == match_tables[j].end())
                        match_tables[j].insert(std::make_pair(bid_sub, std::vector<int>(1, i)));
                    else {
                        this->ds.union_sets(it->second[0], i);
                        if (edge[i] == -1 || order[edge[i]] > order[it->second[0]])
                            edge[i] = it->second[0];
                        it->second.push_back(i);
                    }
                }
            }

            //for (int i = 0; i < edge.size(); ++i) {
            //    printf("Node %d: %s, edge = %d, root = %d\n",
            //        i, binary_to_barcode(bids[i], this->umi_length).c_str(), edge[i], this->ds.find_set(i)
            //    );
            //}

            std::vector<int> thresh(bids.size(), -1);
            for (const int& i: argsorted) {
                int root = this->ds.find_set(i);
                if (edge[i] == -1) continue;
                if (thresh[root] == -1 || thresh[root] < order[edge[i]])
                    thresh[root] = order[edge[i]];
            }

            //for (int i = 0; i < thresh.size(); ++i) {
            //    printf("Node %d: (%s, %d), thresh = %d, parent = %d\n",
            //        i, binary_to_barcode(bids[i], this->umi_length).c_str(), this->ds.get_count()[i], thresh[i], this->ds.get_parent()[i]
            //    );
            //}

            for (const int& i: argsorted) {
                if (thresh[parent[i]] == -1 || order[i] <= thresh[parent[i]]) {
                    //printf("Root node (%s, %d)\n",
                    //    binary_to_barcode(bids[i], this->umi_length).c_str(), counts[i]
                    //);
                    this->ds.set_parent(i, i);
                    this->ds.set_tree_max(i, i);
                } else
                    this->ds.set_parent(i, edge[i]);
            }
        }

        void correct_by_directional() {
            std::vector<uint64_t>& bids = this->ds.get_bid();
            std::vector<int>& counts = this->ds.get_count();
            std::vector<int>& argsorted = this->ds.get_argsorted();
            std::vector<int>& order = this->ds.get_order();

            std::vector<int> edge(bids.size(), -1);

            for (const int& i: argsorted) {
                //printf("Consider (%s, %d):\n",
                //    binary_to_barcode(bids[i], this->umi_length).c_str(), counts[i]
                //);
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = bids[i] & (~aux_arr[j][NNUC]);
                    UMITableIter it = match_tables[j].find(bid_sub);

                    if (it == match_tables[j].end())
                        match_tables[j].insert(std::make_pair(bid_sub, std::vector<int>(1, i)));
                    else {
                        for (int& k: it->second) {
                            if (counts[i] == 1 && counts[k] == 1) {
                                //printf("\tUnion (%s, %d) to (%s, %d)\n",
                                //    binary_to_barcode(bids[i], this->umi_length).c_str(), counts[i],
                                //    binary_to_barcode(bids[k], this->umi_length).c_str(), counts[k]
                                //);
                                this->ds.union_sets(i, k);
                            }
                            if (counts[k] > 1 && counts[k] >= 2 * counts[i] - 1) {
                                if (edge[i] == -1 || order[this->ds.find_union_max(edge[i])] > order[this->ds.find_union_max(k)]) {
                                    //printf("\tAdd edge (%s, %d) to (%s, %d)\n",
                                    //    binary_to_barcode(bids[i], this->umi_length).c_str(), counts[i],
                                    //    binary_to_barcode(bids[k], this->umi_length).c_str(), counts[k]
                                    //);
                                    edge[i] = k;
                                }
                            }
                        }
                        it->second.push_back(i);
                    }
                }

                if (counts[i] > 1 && edge[i] != -1) {
                    //printf("\tUnion (%s, %d) to (%s, %d)\n",
                    //    binary_to_barcode(bids[i], this->umi_length).c_str(), counts[i],
                    //    binary_to_barcode(bids[edge[i]], this->umi_length).c_str(), counts[edge[i]]
                    //);
                    this->ds.union_sets(i, edge[i]);
                }
            }

            for (auto umi_it = argsorted.rbegin(); umi_it != argsorted.rend(); ++umi_it) {
                if (counts[*umi_it] != 1) break;
                if (edge[*umi_it] == -1) continue;

                int root = this->ds.find_union_max(*umi_it);
                //printf("\t%s has root %s\n",
                //    binary_to_barcode(bids[*umi_it], this->umi_length).c_str(),
                //    binary_to_barcode(bids[root], this->umi_length).c_str()
                //);


                if (root == *umi_it) {
                    //printf("\tUnion (%s, %d) to (%s, %d)\n",
                    //    binary_to_barcode(bids[*umi_it], this->umi_length).c_str(), counts[*umi_it],
                    //    binary_to_barcode(bids[edge[*umi_it]], this->umi_length).c_str(), counts[edge[*umi_it]]
                    //);
                    this->ds.union_sets(*umi_it, edge[*umi_it]);
                } else if (edge[root] == -1 || order[this->ds.find_union_max(edge[root])] > order[this->ds.find_union_max(edge[*umi_it])]) {
                    //printf("\tAdd edge (%s, %d) to (%s, %d)\n",
                    //    binary_to_barcode(bids[root], this->umi_length).c_str(), counts[root],
                    //    binary_to_barcode(bids[edge[*umi_it]], this->umi_length).c_str(), counts[edge[*umi_it]]
                    //);
                    edge[root] = edge[*umi_it];
                }
            }
        }

    public:
        UMICorrectSet(int umi_length): umi_length(umi_length) {
            this->match_tables = std::vector<UMITable>(umi_length, UMITable());
        }

        void clear() {
            for (int i = 0; i < this->umi_length; ++i)
                this->match_tables[i].clear();
            this->ds.clear();
        }

        void clear_match_tables() {
            for (int i = 0; i < this->umi_length; ++i)
                this->match_tables[i].clear();
        }

        void build_set(UMI2Count& umi2count, std::string& method) {
            this->ds.init(umi2count);
            if (method == "cluster") correct_by_cluster();
            else if (method == "adjacency") correct_by_adjacency();
            else correct_by_directional();
        }

        UMI2Count get_corrected_umi_counts() {
            return this->ds.get_roots();
        }

};

#endif
