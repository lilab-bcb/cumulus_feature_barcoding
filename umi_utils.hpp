#ifndef UMI_UTILS
#define UMI_UTILS

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "barcode_utils.hpp"

typedef std::unordered_map<uint64_t, int> UMI2Count;
typedef std::unordered_map<uint64_t, int> UMITable;
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
                    std::pair<UMITableIter, bool> ret;
                    ret = match_tables[j].insert(std::make_pair(bid_sub, i));
                    if (!ret.second)
                        this->ds.union_sets(ret.first->second, i);
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
                        has_match = true;
                        if (root == -1 || this->ds.get_count()[root] < this->ds.get_count()[it->second] || (this->ds.get_count()[root] == this->ds.get_count()[it->second] && this->ds.get_bid()[root] > this->ds.get_bid()[it->second]))
                            root = it->second;
                    }
                }

                if (has_match)
                    this->ds.union_sets(i, root);
                else
                    for (int j = 0; j < this->umi_length; ++j) {
                        uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
                        match_tables[j].insert(std::make_pair(bid_sub, i));
                    }
            }
        }

        void correct_by_directional() {
            bool has_match;
            int parent;
            for (const int& i: this->ds.get_order()) {
                has_match = false;
                parent = -1;
                for (int j = 0; j < this->umi_length; ++j) {
                    uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
                    std::pair<UMITableIter, bool> ret;
                    ret = match_tables[j].insert(std::make_pair(bid_sub, i));
                    if (!ret.second && this->ds.get_count()[ret.first->second] >= 2 * this->ds.get_count()[i] - 1) {
                        has_match = true;
                        if (parent == -1 || this->ds.get_count()[parent] < this->ds.get_count()[ret.first->second] || (this->ds.get_count()[parent] == this->ds.get_count()[ret.first->second] && this->ds.get_bid()[parent] > this->ds.get_bid()[ret.first->second]))
                            parent = ret.first->second;
                    }
                }
                if (has_match) this->ds.union_sets(i, parent);
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
