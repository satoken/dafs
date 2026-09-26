/*
 *LinearAlignment.cpp*
 The main code for LinearAlignment: Linear-Time Approximation of 
                                    pairwise RNA sequences alignment
                                    and alignment co-incidence probabilities

 author: Sizhen Li
 created by: 06/2020
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string.h>
#include <map>
#include <stdio.h> 
#include <set> 
#include <cstdint>
#include <stdexcept>

#include "BeamAlign.h"

using namespace std;

int partition(vector<int>& A, int p,int q)
{
    int x= A[p];
    int i=p;
    int j;

    for(j=p+1; j<q; j++)
    {
        if(A[j]<=x)
        {
            i=i+1;
            swap(A[i],A[j]);
        }
    }
    swap(A[i],A[p]);
    return i;
}

void quickSort(vector<int>& A, int p,int q)
{
    int r;
    if(p<q)
    {
        r=partition(A, p,q);
        quickSort(A,p,r);  
        quickSort(A,r+1,q);
    }
}

double BeamAlign::beam_prune(std::unordered_map<int, AlignState> &beamstep){
    scores.clear();
    for (auto &item : beamstep) {
        int ik = item.first;
        AlignState &cand = item.second;
        scores.push_back(make_pair(cand.alpha, ik));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    sort(scores.begin(), scores.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first != rhs.first
             ? lhs.first > rhs.first
             : lhs.second < rhs.second;
    });
    const double threshold = scores[beam-1].first;
    // Keep the beam a hard bound even when several candidates have exactly
    // the cutoff score; retaining all ties can otherwise become quadratic.
    for (size_t rank = beam; rank < scores.size(); ++rank)
        beamstep.erase(scores[rank].second);

    return threshold;
}

void update_if_better(AlignState &state, double newscore, int pre_manner, int manner, unsigned step, unsigned i, unsigned k) {
    if (state.alpha < newscore)
        // ++ nos_set_update;
        state.set(newscore, pre_manner, manner, step, i, k);
};

void update(AlignState &state, double newscore, int pre_manner, int manner, unsigned step, unsigned i, unsigned k) {
    state.set(newscore, pre_manner, manner, step, i, k);
};

double BeamAlign::get_trans_emit_prob(int prev_state, int current_state, int i, int k, double** &trans_probs, double** &emit_probs){
    // get_trans_prob
    if (prev_state < 0 || current_state < 0) return xlog(0.0); // N.B.
    double trans_prob = trans_probs[prev_state][current_state];

    // get_emit_prob
    int i_sym;
	int k_sym;

	// Fix symbols to gaps in case of insertions or boundary conditions.
	if(current_state == 0 || k == 0 || k >= seq2_len)
	{
		// Gap is coded into value 4 in the emission table.
		k_sym = 4; 
	}
	else
	{
		k_sym = nucs2[k];
	}

	if(current_state == 1 || i == 0 || i >= seq1_len)
	{
		// Gap is coded into value 4 in the emission table.
		i_sym = 4;
	}
	else
	{
		i_sym = nucs1[i];
	}

	// Compute the symbol index into emission table using the coded nucleotide values:
	// A->0, C->1, G->2, U->3, T->3, .->4
	// This defines a counting system in base of 5. (25 values.)
	// There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
	int sym_index = i_sym * 5 + k_sym;

	// Check for exceptional cases of start and end symbols.
	// The indices correspond to the start symbol?
	if(i == 0 && k == 0)
	{
		sym_index = 25;
	}

	// The indices correspond to the end symbol?
	if(i == seq1_len && k == seq2_len)
	{
		sym_index = 26;
	}
    double emit_prob = emit_probs[sym_index][current_state];

	return(xlog_mul(emit_prob, trans_prob));
}

double BeamAlign::get_match_prior(int i, int k, bool prior)
{
    if (prior){
        if(i == 0 || k == 0 || i == seq1_len || k == seq2_len)
            return(0.0f); 
        
        // Use custom match score function if provided, otherwise use default
        if (match_score_func_) {
            return xlog(match_score_func_(i, k));
        } else {
            // Default: simple similarity-based score
            // This can be improved with more sophisticated scoring
            return xlog(1.0);  // Neutral score
        }
    }
    
    return 0.0f;
}

void BeamAlign::prepare(string &seq1, string &seq2) {
    seq1_len = static_cast<unsigned>(seq1.length() + 1);
    seq2_len = static_cast<unsigned>(seq2.length() + 1);
    max_len = seq1_len >= seq2_len ? seq1_len : seq2_len;

    // Clean up existing arrays before allocating new ones
    if (nucs1) delete[] nucs1;
    if (nucs2) delete[] nucs2;
    
    nucs1 = new int[seq1_len];
    nucs2 = new int[seq2_len];

    scores.reserve(seq2_len);

    replace(seq1.begin(), seq1.end(), 'T', 'U');
    replace(seq2.begin(), seq2.end(), 'T', 'U');
    for (int i = 0; i < seq1_len; ++i){
        if (i == 0) nucs1[i] = 0;
        else nucs1[i] = GET_ACGU_NUM(seq1[i-1]);
        if (nucs1[i] > 3) {
            int r = rand() % 4;
            nucs1[i] = r;
        }
    }
    
    for (int i = 0; i < seq2_len; ++i){
        if (i == 0) nucs2[i] = 0;
        else nucs2[i] = GET_ACGU_NUM(seq2[i-1]);
        if (nucs2[i] > 3) {
            int r = rand() % 4;
            nucs2[i] = r;
        }
    }
}

void BeamAlign::traceback(vector<char> &aln1, vector<char> &aln2){
    unsigned i_seq1 = seq1_len;
    unsigned i_seq2 = seq2_len;
    unsigned step = i_seq1 + i_seq2;
    unsigned key = i_seq1 * max_len + i_seq2;
    int cur_manner = bestALN[step][key].manner;
    double best_score =  bestALN[step][key].alpha;

    while (true) {
        if ((i_seq1 == 0) && (i_seq2 == 0)) break;
        switch (cur_manner) {
            case 3: // ALIGN_ALN:
                if (i_seq1 < seq1_len) aln1.push_back(GET_NUC(nucs1[i_seq1]));
                if (i_seq2 < seq2_len) aln2.push_back(GET_NUC(nucs2[i_seq2]));
                cur_manner = bestALN[step][key].pre;
                step = step - 2;
                i_seq1 --;
                i_seq2 --;
                key = i_seq1 * max_len + i_seq2;
                break;
            case 1: // ALIGN_INS1:
                aln1.push_back(GET_NUC(nucs1[i_seq1]));
                aln2.push_back('-');
                cur_manner = bestINS1[step][key].pre;
                step = step - 1;
                i_seq1 --;
                key = i_seq1 * max_len + i_seq2;
                break;
            case 2: //ALIGN_INS2:
                aln1.push_back('-');
                aln2.push_back(GET_NUC(nucs2[i_seq2]));
                cur_manner = bestINS2[step][key].pre;
                step = step - 1;
                i_seq2 --;
                key = i_seq1 * max_len + i_seq2;
                break;
            default:  // MANNER_NONE or other cases
                printf("wrong manner at %d, %d: manner %d\n", i_seq1, i_seq2, cur_manner); fflush(stdout);
                assert(false);
        }
    }

    // std::reverse(aln1.begin(), aln1.end());
    // std::reverse(aln2.begin(), aln2.end());
    // for(auto e : aln1) cout << e;
    // cout << endl;
    // for(auto e : aln2) cout << e;
    // cout << endl;

    // return get_aln_similarity(aln1, aln2, '-');
}

void BeamAlign::ml_alignment(string &seq1, string &seq2, vector<char> &aln1, vector<char> &aln2, double** &transprobs, double** &emitprobs, bool prior){
    prepare(seq1, seq2);
        
    double trans_emit_prob;
    for(int s = 0; s < seq1_len + seq2_len; ++s){
        int nuc1 = nucs1[s];
        unordered_map<int, AlignState>& beamALN = bestALN[s];
        unordered_map<int, AlignState>& beamINS1 = bestINS1[s];
        unordered_map<int, AlignState>& beamINS2 = bestINS2[s];

        // initial state
        if (s == 0){
            beamALN[0].alpha = xlog(1.0);
            beamALN[0].manner = 3; // ALIGN_ALN;
            beamALN[0].step = 0;
            beamALN[0].i = 0;
            beamALN[0].k = 0;
        }

        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN};
        for (int i=0; i < beams.size(); i++){
            unordered_map<int, AlignState> &beamstep = *beams[i];

            if (beam > 0 && beamstep.size() > beam) beam_prune(beamstep);
            for (auto &item : beamstep) {
                AlignState &state = item.second;
                unsigned i = state.i;
                unsigned k = state.k;
                unsigned step = state.step;
                unsigned next_key;
                int manner = state.manner;
                int next_manner;
                unsigned next_i, next_k, next_step;
                for (int m = 3; m >= 1; m--){
                    next_step = step;
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = m; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step += 2;
                        next_key = next_i * max_len + next_k;
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            trans_emit_prob = xlog_mul(get_match_prior(next_i, next_k, prior), trans_emit_prob); // match score
                            update_if_better(bestALN[next_step][next_key], xlog_mul(state.alpha, trans_emit_prob), manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;
                        next_key = next_i * max_len + next_k;
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            update_if_better(bestINS1[next_step][next_key], xlog_mul(state.alpha, trans_emit_prob), manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;
                        next_key = next_i * max_len + next_k;
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            update_if_better(bestINS2[next_step][next_key], xlog_mul(state.alpha, trans_emit_prob), manner, next_manner, next_step, next_i, next_k);
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    double forward_score = bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].alpha;
    traceback(aln1, aln2);
}

double BeamAlign::max_alignment(unsigned length1, unsigned length2,
                                std::vector<unsigned>& mapping,
                                MatchScoreFunction match_score)
{
    max_pruned_states_ = 0;
    max_score_ = VALUE_MIN;
    const unsigned stride = length2 + 1;
    const auto key_of = [stride](unsigned i, unsigned k) {
        return static_cast<uint64_t>(i) * stride + k;
    };
    const size_t layer_count = static_cast<size_t>(length1) + length2 + 1;
    if (max_layers_.size() < layer_count)
        max_layers_.resize(layer_count);
    for (size_t layer = 0; layer < layer_count; ++layer)
        max_layers_[layer].entries.clear();
    auto& layers = max_layers_;

    // MaxNode scores cannot be NaN here: nodes start at the finite VALUE_MIN
    // sentinel, and update_max stores a candidate only after node.score <
    // score succeeds.  Thus a NaN callback result (or NaN sum) is rejected,
    // and the original finite score/key order remains a strict weak order.
    const auto ranking_before = [](double lhs_score, uint64_t lhs_key,
                                   double rhs_score, uint64_t rhs_key) {
        if (lhs_score != rhs_score)
            return lhs_score > rhs_score;
        return lhs_key < rhs_key;
    };

    const auto node_at = [](const MaxLayer& layer, uint64_t key) {
        const auto position = std::lower_bound(
            layer.entries.begin(), layer.entries.end(), key,
            [](const MaxEntry& entry, uint64_t candidate) {
                return entry.key < candidate;
            });
        return position == layer.entries.end() || position->key != key
             ? static_cast<const MaxNode*>(nullptr)
             : &position->node;
    };

    layers[0].entries.push_back(MaxEntry{0, MaxNode{}});
    layers[0].entries.front().node.score = 0.0;

    const auto update_max = [](MaxNode& node, double score,
                               uint64_t previous, char operation) {
        if (node.score < score) {
            node.score = score;
            node.previous = previous;
            node.operation = operation;
        }
    };

    // At the start of step s, layer s+1 contains only match offers from
    // step s-1.  The current layer supplies two sorted gap streams: X moves
    // a source one row down, while Y leaves its row unchanged.  Merge those
    // streams with the existing match stream, applying the same strict update
    // in offer order (existing match, X, then Y) for equal destination keys.
    const auto merge_gap_successors = [&](unsigned step) {
        MaxLayer& source = layers[step];
        MaxLayer& target = layers[step + 1];

        // Keep the old target entries in one reusable scratch vector while
        // rebuilding target in key order.  No row-sized scratch array or
        // per-layer allocation is needed.
        max_existing_.swap(target.entries);
        target.entries.clear();
        target.entries.reserve(max_existing_.size() + source.entries.size() * 2);

        size_t existing_index = 0;
        // On a fixed diagonal, source keys are ordered by i.  Therefore X
        // successors are a valid prefix (i < length1), while Y successors
        // are a valid suffix (k < length2).  Find those bounds once so the
        // merge itself does not repeatedly divide/modulo every source key.
        size_t x_index = 0;
        size_t x_end = source.entries.size();
        if (step >= length1) {
            const uint64_t first_invalid_x =
                key_of(length1, step - length1);
            x_end = static_cast<size_t>(std::lower_bound(
                source.entries.begin(), source.entries.end(), first_invalid_x,
                [](const MaxEntry& entry, uint64_t candidate) {
                    return entry.key < candidate;
                }) - source.entries.begin());
        }

        size_t y_index = 0;
        if (step >= length2) {
            const uint64_t last_invalid_y =
                key_of(step - length2, length2);
            y_index = static_cast<size_t>(std::upper_bound(
                source.entries.begin(), source.entries.end(), last_invalid_y,
                [](uint64_t candidate, const MaxEntry& entry) {
                    return candidate < entry.key;
                }) - source.entries.begin());
        }

        while (true) {
            if (existing_index == max_existing_.size() &&
                x_index == x_end &&
                y_index == source.entries.size())
                break;

            const bool has_existing = existing_index < max_existing_.size();
            const bool has_x = x_index < x_end;
            const bool has_y = y_index < source.entries.size();

            uint64_t next_key = std::numeric_limits<uint64_t>::max();
            if (has_existing)
                next_key = std::min(next_key,
                                    max_existing_[existing_index].key);
            if (has_x)
                next_key = std::min(next_key,
                                    source.entries[x_index].key + stride);
            if (has_y)
                next_key = std::min(next_key,
                                    source.entries[y_index].key + 1);

            // Construct a default node first.  This deliberately preserves
            // the old node_for-then-update behavior when the first offer is
            // exactly VALUE_MIN or NaN.
            MaxNode node;
            if (has_existing &&
                max_existing_[existing_index].key == next_key) {
                node = max_existing_[existing_index].node;
                ++existing_index;
            }

            // X from source i-1 is offered before Y from source i for the
            // same destination.  update_max is strict, so this also keeps
            // the original winner on equal scores.
            if (has_x && source.entries[x_index].key + stride == next_key) {
                const MaxEntry& entry = source.entries[x_index];
                update_max(node, entry.node.score, entry.key, 'X');
                ++x_index;
            }
            if (has_y && source.entries[y_index].key + 1 == next_key) {
                const MaxEntry& entry = source.entries[y_index];
                update_max(node, entry.node.score, entry.key, 'Y');
                ++y_index;
            }

            target.entries.push_back(MaxEntry{next_key, node});
        }
        max_existing_.clear();
    };

    for (unsigned step = 0; step < length1 + length2; ++step) {
        auto& layer = layers[step];
        if (beam > 0 && layer.entries.size() > static_cast<size_t>(beam)) {
            max_pruned_states_ += layer.entries.size() - static_cast<size_t>(beam);
            max_ranked_.clear();
            max_ranked_.reserve(layer.entries.size());
            for (const MaxEntry& entry : layer.entries)
                max_ranked_.emplace_back(entry.node.score, entry.key);
            const auto ranking_order = [&](const auto& lhs, const auto& rhs) {
                return ranking_before(lhs.first, lhs.second,
                                      rhs.first, rhs.second);
            };
            std::nth_element(
                max_ranked_.begin(), max_ranked_.begin() + beam,
                max_ranked_.end(), ranking_order);

            // The element at `beam` is the first excluded rank.  Since the
            // key makes every rank unique, filtering the already key-sorted
            // layer by this cutoff keeps exactly the same top-beam entries
            // without sorting or separately storing their keys.
            const auto cutoff = max_ranked_[beam];
            size_t write = 0;
            for (size_t read = 0; read < layer.entries.size(); ++read) {
                const MaxEntry& entry = layer.entries[read];
                if (!ranking_before(entry.node.score, entry.key,
                                    cutoff.first, cutoff.second))
                    continue;
                if (write != read)
                    layer.entries[write] = entry;
                ++write;
            }
            layer.entries.resize(write);
        }

        // layer[step+2] is empty here: only the current layer can offer a
        // match to it.  Its successor keys are monotone in source key order,
        // so append directly instead of inserting into a sorted vector.  The
        // final layer cannot have a match successor, so avoid forming the
        // otherwise out-of-range step+2 reference there.
        MaxLayer* match_target = nullptr;
        if (step + 2 < layer_count) {
            match_target = &layers[step + 2];
            match_target->entries.reserve(layer.entries.size());
        }
        for (const MaxEntry& entry : layer.entries) {
            const uint64_t key = entry.key;
            const MaxNode& node = entry.node;
            const unsigned i = key / stride;
            const unsigned k = key % stride;
            if (i < length1 && k < length2) {
                const uint64_t next = key_of(i + 1, k + 1);
                MaxNode match_node;
                update_max(match_node, node.score + match_score(i, k),
                           key, 'M');
                match_target->entries.push_back(MaxEntry{next, match_node});
            }
        }
        merge_gap_successors(step);
    }

    mapping.assign(length1, -1u);
    uint64_t key = key_of(length1, length2);
    unsigned step = length1 + length2;
    const MaxNode* terminal = node_at(layers[step], key);
    if (terminal == nullptr) return VALUE_MIN;
    const double score = terminal->score;
    max_score_ = score;
    while (step > 0) {
        const MaxNode* node = node_at(layers[step], key);
        if (node == nullptr)
            throw std::out_of_range("MaxLayer::at");
        const unsigned i = key / stride;
        const unsigned k = key % stride;
        if (node->operation == 'M') {
            mapping[i - 1] = k - 1;
            step -= 2;
        } else {
            --step;
        }
        key = node->previous;
    }
    return score;
}


double BeamAlign::forward(string seq1, string seq2, double** &trans_probs, double** &emit_probs, bool prior){
    // Clean up beam arrays before prepare
    if (bestINS1) delete[] bestINS1;
    if (bestINS2) delete[] bestINS2;
    if (bestALN) delete[] bestALN;
    
    // Prepare the sequences (this will handle nucs1/nucs2)
    prepare(seq1, seq2);
    
    // Allocate new beam arrays
    bestALN = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];
    bestINS1 = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];
    bestINS2 = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];

    scores.reserve(seq2_len);

    // initial state
    bestALN[0][0].alpha = xlog(1.0);
    bestALN[0][0].manner = 3; // ALIGN_ALN;
    bestALN[0][0].step = 0;
    bestALN[0][0].i = 0;
    bestALN[0][0].k = 0;

    double trans_emit_prob;
    for(int s = 0; s < seq1_len + seq2_len; ++s){
        unordered_map<int, AlignState> &beamALN = bestALN[s];
        unordered_map<int, AlignState> &beamINS1 = bestINS1[s];
        unordered_map<int, AlignState> &beamINS2 = bestINS2[s];

        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN};

        for (int i=0; i < beams.size(); i++){
            unordered_map<int, AlignState> &beamstep = *beams[i];
            if (beam > 0 && beamstep.size() > beam) beam_prune(beamstep);
            for (auto &item : beamstep) {
                AlignState &state = item.second;
                unsigned i = state.i;
                unsigned k = state.k;
                unsigned step = state.step;
                int manner = state.manner;

                int next_manner;
                unsigned next_i, next_k, next_step, next_key;
                for (int m = 1; m <= 3; m++){
                    next_step = step;
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = m; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step += 2;
                        
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, trans_probs, emit_probs);
                            trans_emit_prob = xlog_mul(get_match_prior(next_i, next_k, prior), trans_emit_prob);
                            double newscore = xlog_sum(bestALN[next_step][next_key].alpha, xlog_mul(state.alpha, trans_emit_prob));
                            update(bestALN[next_step][next_key], newscore, manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;

                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))){
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, trans_probs, emit_probs);
                            double newscore = xlog_sum(bestINS1[next_step][next_key].alpha, xlog_mul(state.alpha, trans_emit_prob));
                            update(bestINS1[next_step][next_key], newscore, manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;

                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, trans_probs, emit_probs);
                            double newscore = xlog_sum(bestINS2[next_step][next_key].alpha, xlog_mul(state.alpha, trans_emit_prob));
                            update(bestINS2[next_step][next_key], newscore, manner, next_manner, next_step, next_i, next_k);
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    double forward_score = bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].alpha;
    // printf("forward score: %.5f\n", forward_score);
    return forward_score;
}

double BeamAlign::backward(double** &transprobs, double** &emitprobs, bool prior){
    // Initialize the final state beta = 1.0 (log(1.0) = 0.0)
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].beta = xlog(1.0);
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].manner = 3; //ALIGN_ALN;
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].step = seq1_len + seq2_len;
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].i = seq1_len;
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].k = seq2_len;
    
    // Beta values are already initialized to xlog(0) by AlignState constructor
    // No need to explicitly reinitialize them

    for(int s = seq1_len + seq2_len - 2; s >= 0; --s) {
        double trans_emit_prob;
        unordered_map<int, AlignState> &beamALN = bestALN[s];
        unordered_map<int, AlignState> &beamINS1 = bestINS1[s];
        unordered_map<int, AlignState> &beamINS2 = bestINS2[s];
        vector<unordered_map<int, AlignState>*> beams{&beamALN, &beamINS1, &beamINS2}; // N.B.

        for (int i=0; i < beams.size(); i++){
            unordered_map<int, AlignState> &beamstep = *beams[i]; // N.B.
            for (auto &item : beamstep) {
                AlignState &state = item.second;
                unsigned i = state.i;
                unsigned k = state.k;
                unsigned step = state.step;
                int manner = state.manner;

                int next_manner;
                unsigned next_i, next_k, next_step, next_key;
                for (int m = 1; m <= 3; m++){
                    next_step = step;
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = m; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step += 2;

                        if (((next_i == seq1_len) && (next_k == seq2_len)) || ((next_i < seq1_len) && (next_k < seq2_len))) {
                            next_key = next_i * max_len + next_k;
                            // Check if next state exists before accessing its beta value
                            if (bestALN[next_step].find(next_key) != bestALN[next_step].end()) {
                                trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                                trans_emit_prob = xlog_mul(get_match_prior(next_i, next_k, prior), trans_emit_prob);
                                double next_beta = bestALN[next_step][next_key].beta;
                                if (!std::isnan(trans_emit_prob) && !std::isnan(next_beta)) {
                                    state.beta = xlog_sum(state.beta, xlog_mul(next_beta, trans_emit_prob));
                                }
                            }
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            next_key = next_i * max_len + next_k;
                            // Check if next state exists before accessing its beta value
                            if (bestINS1[next_step].find(next_key) != bestINS1[next_step].end()) {
                                trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                                double next_beta = bestINS1[next_step][next_key].beta;
                                if (!std::isnan(trans_emit_prob) && !std::isnan(next_beta)) {
                                    state.beta = xlog_sum(state.beta, xlog_mul(next_beta, trans_emit_prob));
                                }
                            }
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            next_key = next_i * max_len + next_k;
                            // Check if next state exists before accessing its beta value
                            if (bestINS2[next_step].find(next_key) != bestINS2[next_step].end()) {
                                trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                                double next_beta = bestINS2[next_step][next_key].beta;
                                if (!std::isnan(trans_emit_prob) && !std::isnan(next_beta)) {
                                    state.beta = xlog_sum(state.beta, xlog_mul(next_beta, trans_emit_prob));
                                }
                            }
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    double back_score = bestALN[0][0].beta;
    // printf("backward score: %.5f\n", back_score);
    return back_score;
}

std::unordered_map<int, aln_ret>*  BeamAlign::cal_align_prob(double forward_score, double threshold, std::unordered_map<int, aln_ret>* &aln_results){
    // Always create a new array for this sequence pair - ignore input parameter
    aln_results = new std::unordered_map<int, aln_ret>[seq1_len];
    
    double aln_prob, ins1_prob, ins2_prob;
   
    for (int s = 0; s < seq1_len + seq2_len; s++){
        for(auto &item : bestALN[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            // Add bounds checking before accessing aln_results array
            if (i >= 0 && i < (int)seq1_len && k >= 0) {
                aln_prob = xlog_div(xlog_mul(state.alpha, state.beta), forward_score);
                if (aln_prob > float(-9.91152)) {
                    aln_results[i][k].prob = xlog_sum(aln_results[i][k].prob, aln_prob);
                    aln_results[i][k].aln_prob = aln_prob;
                }
            }
        }

        for(auto &item : bestINS1[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            // Add bounds checking before accessing aln_results array
            if (i >= 0 && i < (int)seq1_len && k >= 0) {
                ins1_prob = xlog_div(xlog_mul(state.alpha, state.beta), forward_score);
                if (ins1_prob > float(-9.91152)) aln_results[i][k].prob = xlog_sum(aln_results[i][k].prob, ins1_prob);
            }
        }

        for(auto &item : bestINS2[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            // Add bounds checking before accessing aln_results array
            if (i >= 0 && i < (int)seq1_len && k >= 0) {
                ins2_prob = xlog_div(xlog_mul(state.alpha, state.beta), forward_score);
                if (ins2_prob > float(-9.91152)) aln_results[i][k].prob = xlog_sum(aln_results[i][k].prob, ins2_prob);
            }
        }
    }

    vector<int> cands;
    for (int i = 1; i < seq1_len; i++) {
        cands.clear();
        for(auto &item : aln_results[i]){
            int k = item.first;
            double prob = item.second.prob;
            if (prob < threshold) {
                cands.push_back(k);
            } else {
                item.second.prob = std::exp(prob);
                item.second.aln_prob = std::exp(item.second.aln_prob);
            }
        }
     
        for (auto &k : cands) {
            aln_results[i].erase(k);
        }
    }
    return aln_results;
}

BeamAlign::BeamAlign(int beam_size)
    : beam(beam_size), match_score_func_(nullptr), seq1_len(0), seq2_len(0), max_len(0) {
    init_pointers();
}

BeamAlign::~BeamAlign(){
    cleanup_arrays();
}

void BeamAlign::init_pointers() {
    bestINS1 = nullptr;
    bestINS2 = nullptr;
    bestALN = nullptr;
    nucs1 = nullptr;
    nucs2 = nullptr;
}

void BeamAlign::cleanup_arrays() {
    if (bestINS1) {
        delete[] bestINS1;
        bestINS1 = nullptr;
    }
    if (bestINS2) {
        delete[] bestINS2;
        bestINS2 = nullptr;
    }
    if (bestALN) {
        delete[] bestALN;
        bestALN = nullptr;
    }
    if (nucs1) {
        delete[] nucs1;
        nucs1 = nullptr;
    }
    if (nucs2) {
        delete[] nucs2;
        nucs2 = nullptr;
    }
}
