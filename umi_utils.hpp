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
        std::vector<int> max_count;
        std::vector<int> order;

        std::vector<int> path;
    public:
        DisjointSet() {}

        void init(UMI2Count& umi2count, bool do_sort) {
            int n = umi2count.size();

            bid = std::vector<uint64_t>(n, 0);
            count = std::vector<int>(n, 0);
            parent = std::vector<int>(n, -1);
            rank = std::vector<int>(n, 0);
            max_count = std::vector<int>(n, 0);
            order = std::vector<int>(n);

            int i = 0;
            for (auto& kv: umi2count) {
                bid[i] = kv.first;
                count[i] = kv.second;
                parent[i] = i;
                max_count[i] = i;
                order[i] = i;
                ++i;
            }

            if (do_sort)
                std::sort(order.begin(), order.end(),
                    [this](size_t i1, size_t i2) {
                        if (count[i1] != count[i2])
                            return count[i1] > count[i2];
                        else
                            return bid[i1] < bid[i2];
                });
        }

        std::vector<uint64_t>& get_bid() {
            return bid;
        }

        std::vector<int>& get_order() {
            return order;
        }

        std::vector<int>& get_count() {
            return count;
        }

        std::vector<int>& get_max_count() {
            return max_count;
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
                    if ((count[max_count[a]] > count[max_count[b]]) || (count[max_count[a]] == count[max_count[b]] && bid[max_count[a]] < bid[max_count[b]]))
                        max_count[b] = max_count[a];
                } else if (rank[a] > rank[b]) {
                    parent[b] = a;
                    if ((count[max_count[a]] < count[max_count[b]]) || (count[max_count[a]] == count[max_count[b]] && bid[max_count[a]] > bid[max_count[b]]))
                        max_count[a] = max_count[b];
                }
                if (rank[a] == rank[b])
                    ++rank[b];
            }
        }

        int find_count_union_max_count(int a) {
            return count[max_count[find_set(a)]];
        }

        uint64_t find_bid_union_max_count(int a) {
            return bid[max_count[find_set(a)]];
        }

        bool connected(int a, int b) {
            return find_set(a) == find_set(b);
        }

        UMI2Count get_roots() {
            UMI2Count res;
            UMI2Count::iterator it;

            for (int i = 0; i < bid.size(); ++i) {
                int root = find_set(i);
                it = res.find(bid[max_count[root]]);
                if (it != res.end())
                    it->second += count[i];
                else
                    res.insert(std::make_pair(bid[max_count[root]], count[i]));
            }
            return res;
        }

        void clear() {
            bid.clear();
            count.clear();
            parent.clear();
            rank.clear();
            max_count.clear();
            order.clear();
        }
};

class UMICorrectSet {
    private:
        std::vector<UMITable> match_tables;
        DisjointSet ds;
        int umi_length;

        void correct_by_cluster() {
            for (const int& i: this->ds.get_order())
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
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
            bool has_match;
            int root;
            for (const int& i: this->ds.get_order()) {
                has_match = false;
                root = -1;
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
                    UMITableIter it = match_tables[j].find(bid_sub);
                    if (it != match_tables[j].end()) {
                        if (has_match) {
                            // Current UMI connects to multiple roots. Set itself as root.
                            has_match = false;
                            break;
                        }
                        has_match = true;
                        root = it->second[0];
                    }
                }

                if (has_match)
                    this->ds.union_sets(i, root);
                else
                    for (int j = 0; j < this->umi_length; ++j) {
                        uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
                        match_tables[j].insert(std::make_pair(bid_sub, std::vector<int>(1, i)));
                    }
            }
        }

        void correct_by_directional() {
            bool has_match;
            int parent;

            int pos_one = -1;
            for (int i = 0; i < this->ds.get_order().size(); ++i)  // TODO: Make it binary search?
                if (this->ds.get_count()[this->ds.get_order()[i]] == 1) {
                    pos_one = i;
                    break;
                }

            // Process UMIs with count 1 as cluster method
            for (auto umi_it = this->ds.get_order().begin() + pos_one; umi_it != this->ds.get_order().end(); ++umi_it) {
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = this->ds.get_bid()[*umi_it] & (~aux_arr[j][NNUC]);
                    std::pair<UMITableIter, bool> ret;
                    ret = match_tables[j].insert(std::make_pair(bid_sub, std::vector<int>(1, *umi_it)));
                    if (!ret.second) {
                        this->ds.union_sets(ret.first->second[0], *umi_it);
                    }
                }
            }

            std::unordered_map<int, std::vector<int>> groups_one;
            int root;
            for (auto umi_it = this->ds.get_order().begin() + pos_one; umi_it != this->ds.get_order().end(); ++umi_it) {
                root = this->ds.find_set(*umi_it);
                if (root == *umi_it)
                    groups_one.insert(std::make_pair(*umi_it, std::vector<int>(1, *umi_it)));
                else
                    groups_one[root].push_back(*umi_it);
            }
            for (auto& kv: groups_one) {
                printf("Root %s: ", binary_to_barcode(this->ds.get_bid()[kv.first], this->umi_length).c_str());
                for (int& i: kv.second) {
                    printf("%s, ", binary_to_barcode(this->ds.get_bid()[i], this->umi_length).c_str());
                }
                printf("\n");
            }
            printf("\n");

            if (pos_one == 0) return;

            // Process UMIs with count > 1
            clear_match_tables();
            for (const int& i: this->ds.get_order()) {
                if (this->ds.get_count()[i] == 1) break;

                has_match = false;
                parent = -1;
                printf("Consider (%s, %d):\n",
                    binary_to_barcode(this->ds.get_bid()[i], this->umi_length).c_str(), this->ds.get_count()[i]
                );
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
                    UMITableIter it = match_tables[j].find(bid_sub);
                    if (it == match_tables[j].end()) {
                        match_tables[j].insert(std::make_pair(bid_sub, std::vector<int>(1, i)));
                        //printf("\tInsert %s into hash table %d\n", binary_to_barcode(bid_sub, this->umi_length).c_str(), j);
                    } else {
                        for (int& k: it->second) {
                            if (this->ds.get_count()[k] >= 2 * this->ds.get_count()[i] - 1) {
                                if (parent == -1
                                 || (this->ds.find_count_union_max_count(parent) < this->ds.find_count_union_max_count(k))
                                 || (this->ds.find_count_union_max_count(parent) == this->ds.find_count_union_max_count(k) && this->ds.find_bid_union_max_count(parent) > this->ds.find_bid_union_max_count(k))
                                ) {
                                    if (parent == -1)
                                        printf("\tSet parent to (%s, %d) of max count %d\n",
                                            binary_to_barcode(this->ds.get_bid()[k], this->umi_length).c_str(), this->ds.get_count()[k], this->ds.find_count_union_max_count(k)
                                        );
                                    else
                                        printf("\tSwitch parent from (%s, %d) of max count %d to (%s, %d) of max count %d\n",
                                            binary_to_barcode(this->ds.get_bid()[parent], this->umi_length).c_str(), this->ds.get_count()[parent], this->ds.find_count_union_max_count(parent),
                                            binary_to_barcode(this->ds.get_bid()[k], this->umi_length).c_str(), this->ds.get_count()[k], this->ds.find_count_union_max_count(k)
                                        );
                                    has_match = true;
                                    parent = k;
                                }
                            }
                        }

                        //printf("\tInsert %s into hash table %d\n", binary_to_barcode(bid_sub, this->umi_length).c_str(), j);
                        it->second.push_back(i);
                    }
                }
                if (has_match) {
                    printf("Union (%s, %d) with (%s, %d)\n",
                            binary_to_barcode(this->ds.get_bid()[i], this->umi_length).c_str(), this->ds.get_count()[i],
                            binary_to_barcode(this->ds.get_bid()[parent], this->umi_length).c_str(), this->ds.get_count()[parent]
                        );
                    this->ds.union_sets(i, parent);
                }
            }

            // Merge the two graphs together
            for (auto& kv: groups_one) {
                has_match = false;
                parent = -1;
                for (const int& i: kv.second) {
                    for (int j = 0; j < this->umi_length; ++j) {
                        uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
                        UMITableIter it = match_tables[j].find(bid_sub);
                        if (it != match_tables[j].end())
                            for (int&k: it->second) {
                                if (this->ds.get_count()[k] >= 2 * this->ds.get_count()[i] - 1) {
                                    if (parent == -1
                                        || (this->ds.find_count_union_max_count(parent) < this->ds.find_count_union_max_count(k))
                                        || (this->ds.find_count_union_max_count(parent) == this->ds.find_count_union_max_count(k) && this->ds.find_bid_union_max_count(parent) > this->ds.find_bid_union_max_count(k))
                                       ) {
                                           if (parent == -1)
                                               printf("\tSet parent to (%s, %d) of max count %d\n",
                                                   binary_to_barcode(this->ds.get_bid()[k], this->umi_length).c_str(), this->ds.get_count()[k], this->ds.find_count_union_max_count(k)
                                               );
                                           else
                                               printf("\tSwitch parent from (%s, %d) of max count %d to (%s, %d) of max count %d\n",
                                                   binary_to_barcode(this->ds.get_bid()[parent], this->umi_length).c_str(), this->ds.get_count()[parent], this->ds.find_count_union_max_count(parent),
                                                   binary_to_barcode(this->ds.get_bid()[k], this->umi_length).c_str(), this->ds.get_count()[k], this->ds.find_count_union_max_count(k)
                                               );
                                           has_match = true;
                                           parent = k;
                                       }
                                }
                            }
                    }
                }
                if (has_match)
                    this->ds.union_sets(kv.first, parent);
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
            if (method == "cluster") {
                this->ds.init(umi2count, false);
                correct_by_cluster();
            } else if (method == "adjacency") {
                this->ds.init(umi2count, true);
                correct_by_adjacency();
            } else {
                this->ds.init(umi2count, true);
                correct_by_directional();
            }
        }

        UMI2Count get_corrected_umi_counts() {
            return this->ds.get_roots();
        }

};

#endif
