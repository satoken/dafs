// $Id:$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cstring>
#include <unistd.h>
#include <vector>
#include <queue>
#include <stack>
#include "cmdline.h"
#include "fa.h"
#include "fold.h"
#include "nussinov.h"
#include "align.h"
#include "needleman_wunsch.h"
#include "ip.h"
#include "typedefs.h"

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
};
};

#define FOREACH(itr, i, v) for (itr i=(v).begin(); i!=(v).end(); ++i)
#define CUTOFF 0.01

#define SPARSE_ALIGNMENT
#define SPARSE_FOLDING
#define SPARSE_UPDATE

class Dusaf
{
private:
  // nodes in the guide tree
  typedef std::pair<float, std::pair<uint,uint> > node_t;
  // alignments
  typedef std::vector< std::pair<uint, std::vector<bool> > > ALN;
  // indices for consensus base pairs
  typedef std::pair< std::pair<uint,uint>, std::pair<uint,uint> > CBP;
  
public:
  Dusaf()
    : w_pct_a_(0.0),
      w_pct_s_(0.0),
      w_pct_f_(0.0),
      n_refinement_(0),
      t_max_(100),
      th_a_(1.0/(1+99)),
      th_s_(1.0/(1+4)),
      w_(1.0),
      eta0_(0.5),
      a_model_(NULL),
      a_decoder_(NULL),
      s_model_(NULL),
      s_decoder_(NULL),
      s_decoder1_(NULL),
      use_alifold_(false),
      // use_bpscore_(false),
      verbose_(0)
  {
  }

  ~Dusaf()
  {
    if (a_model_) delete a_model_;
    if (a_decoder_) delete a_decoder_;
    if (s_model_) delete s_model_;
    if (s_decoder_) delete s_decoder_;
    if (s_decoder1_) delete s_decoder1_;
  }

  Dusaf& parse_options(int& argc, char**& argv);
  int run();

private:
  void relax_matching_probability();
  void relax_basepairing_probability();
  void relax_fourway_consistency();
  void build_tree();
  void print_tree(std::ostream& os, int i) const;
  void project_alignment(ALN& aln, const ALN& aln1, const ALN& aln2, const VU& z) const;
  void project_secondary_structure(VU& xx, VU& yy, const VU& x, const VU& y, const VU& z) const;
  void average_matching_probability(VVF& posterior, const ALN& aln1, const ALN& aln2) const;
  void average_basepairing_probability(VVF& posterior, const ALN& aln) const;
  void align_alignments(ALN& aln, const ALN& aln1, const ALN& aln2) const;
  float align_alignments(VU& ss, ALN& aln, const ALN& aln1, const ALN& aln2) const;
  void solve(VU& x, VU& y, VU& z, const VVF& p_x, const VVF& p_y, const VVF& p_z,
             const ALN& aln1, const ALN& aln2) const
  {
#if defined(WITH_GLPK) || defined(WITH_CPLEX) || defined(WITH_GUROBI)
    t_max_!=0 ?
      solve_by_dd(x, y, z, p_x, p_y, p_z, aln1, aln2) : 
      solve_by_ip(x, y, z, p_x, p_y, p_z, aln1, aln2) ;
#else
    solve_by_dd(x, y, z, p_x, p_y, p_z, aln1, aln2);
#endif
  }
  void solve_by_dd(VU& x, VU& y, VU& z, const VVF& p_x, const VVF& p_y, const VVF& p_z,
                   const ALN& aln1, const ALN& aln2) const;
  void solve_by_ip(VU& x, VU& y, VU& z, const VVF& p_x, const VVF& p_y, const VVF& p_z,
                   const ALN& aln1, const ALN& aln2) const;
  void align(ALN& aln, int ch) const;
  float align(VU& ss, ALN& aln, int ch) const;
  float refine(VU& ss, ALN& aln) const;
  void output_verbose(const VU& x, const VU& y, const VU& z, const ALN& aln1, const ALN& aln2) const;
  void output(std::ostream& os, const ALN& aln, const VU& ss) const;
  void output(std::ostream& os, const ALN& aln) const;
  void output(std::ostream& os, ALN::const_iterator b, ALN::const_iterator e) const;

private:
  float w_pct_a_;               // the weight of PCT for alignment matching probabilities
  float w_pct_s_;               // the weight of PCT for base-pairing probabilities
  float w_pct_f_;               // the weight of four-way PCT
  uint n_refinement_;           // the number of the iterative refinement
  uint t_max_;                  // the maximum number of the iteration of the subgradient update
  float th_a_;                  // the threshold for base-pairing probabilities
  float th_s_;                  // the threshold for alignment matching probabilities
  float w_;                     // the weight for base pairs in the objective function
  float eta0_;                  // the initial step width of the subgradient update
  Align::Model* a_model_;       // alignment model
  Align::Decoder* a_decoder_;   // alignment decoder
  Fold::Model* s_model_;        // folding model
  Fold::Decoder* s_decoder_;    // folding decoder
  Fold::Decoder* s_decoder1_;    // folding decoder for the conclusive folding
  std::vector<Fasta> fa_;       // input sequences
  std::vector<std::vector<MP> > mp_; // alignment matching probability matrices
  std::vector<BP> bp_;          // base-pairing probability matrices
  VVF sim_;                     // simalarity matrix between input sequences
  std::vector<node_t> tree_;    // guide tree
  bool use_alifold_;
  // bool use_bpscore_;
  uint verbose_;
};

static
void
transpose_mp(const MP& mp, MP& mp_trans, uint x, uint y)
{
  assert(mp.size()==x);
  mp_trans.resize(y);
  for (uint i=0; i!=mp.size(); ++i)
    FOREACH(SV::const_iterator, jj, mp[i])
    {
      const uint j=jj->first;
      const float p=jj->second;
      mp_trans[j].push_back(std::make_pair(i, p));
    }
  for (uint j=0; j!=mp_trans.size(); ++j)
    sort(mp_trans[j].begin(), mp_trans[j].end());
}

#if 0
static
void
print_mp(std::ostream& os, const MP& mp)
{
  for (uint i=0; i!=mp.size(); ++i)
  {
    os << i << ":";
    FOREACH (SV::const_iterator, v, mp[i])
      os << " " << v->first << ":" << v->second;
    os << std::endl;
  }
}

static
void
print_bp(std::ostream& os, const BP& bp)
{
  for (uint i=0; i!=bp.size(); ++i)
  {
    os << i << ":";
    FOREACH (SV::const_iterator, v, bp[i])
      os << " " << v->first << ":" << v->second;
    os << std::endl;
  }
}

static
void
print_matching_probability(std::ostream& os, const VVF& p)
{
  for (uint i=0; i!=p.size(); ++i)
  {
    std::cout << i << ":";
    for (uint k=0; k!=p[i].size(); ++k)
      if (p[i][k]>CUTOFF)
        os << " " << k << ":" << p[i][k];
    os << std::endl;
  }
}

static
void
print_basepairing_probability(std::ostream& os, const VVF& p)
{
  for (uint i=0; i!=p.size(); ++i)
  {
    std::cout << i << ":";
    for (uint j=i+1; j!=p.size(); ++j)
      if (p[i][j]>CUTOFF)
        os << " " << j << ":" << p[i][j];
    os << std::endl;
  }
}
#endif

void
Dusaf::
relax_matching_probability()
{
  const uint N=fa_.size();
  std::vector<std::vector<MP> > mp(N, std::vector<MP>(N));
  assert(mp_.size()==N);
  assert(mp_[0].size()==N);
  for (uint x=0; x!=N-1; ++x)
  {
    const uint L1=fa_[x].size();
    for (uint y=x+1; y!=N; ++y)
    {
      const uint L2=fa_[y].size();
      VVF posterior(L1, VF(L2, 0.0));
      assert(L1==mp_[x][y].size());

      float sum_w=0.0;
      for (uint z=0; z!=N; ++z)
      {
        const uint L3=fa_[z].size();
        assert(L3==mp_[z][x].size());
        assert(L3==mp_[z][y].size());
        float w = sim_[z][x] * sim_[z][y];;
        if (w_pct_a_<0.0)
          w *= 1.0/N;
        else if (z==x || z==y)
          w *= (1.0-w_pct_a_)/2;
        else
          w *= w_pct_a_/(N-2);
        sum_w += w;
        for (uint k=0; k!=L3; ++k)
        {
          FOREACH (SV::const_iterator, ii, mp_[z][x][k])
          {
            FOREACH (SV::const_iterator, jj, mp_[z][y][k])
            {
              const uint i=ii->first, j=jj->first;
              const float p_ik=ii->second, p_jk=jj->second;
              assert(i<L1); assert(j<L2);
              posterior[i][j] += p_ik * p_jk * w;
            }
          }
        }
      }
    
      mp[x][y].resize(L1);
      for (uint i=0; i!=L1; ++i)
      {
        for (uint j=0; j!=L2; ++j)
        {
          float v=posterior[i][j]/sum_w;
          if (v>CUTOFF)
            mp[x][y][i].push_back(std::make_pair(j, v));
        }
      }
      transpose_mp(mp[x][y], mp[y][x], L1, L2);
    }
  }

  for (uint x=0; x!=N; ++x)
  {
    mp[x][x].resize(fa_[x].size());
    for (uint i=0; i!=fa_[x].size(); ++i)
      mp[x][x][i].push_back(std::make_pair(i, 1.0f));
  }
  std::swap(mp_, mp);
}

void
Dusaf::
relax_basepairing_probability()
{
  const uint N=bp_.size();
  std::vector<BP> bp(N);
  for (uint x=0; x!=N; ++x)
  {
    const uint L1=bp_[x].size();
    VVF p(L1, VF(L1, 0.0));

    float sum_w = 0.0;
    for (uint y=0; y!=N; ++y)
    {
      const uint L2=bp_[y].size();
      assert(L2==mp_[y][x].size());
      float w = sim_[y][x];
      if (w_pct_s_<0.0)
        w *= 1.0/N;
      else if (y==x)
        w *= 1.0-w_pct_s_;
      else
        w *= w_pct_s_/(N-1);
      sum_w += w;
      for (uint k=0; k!=L2; ++k)
      {
        FOREACH(SV::const_iterator, ll, bp_[y][k])
        {
          const uint l=ll->first;
          const float p_kl=ll->second;
          FOREACH(SV::const_iterator, ii, mp_[y][x][k])
          {
            const uint i=ii->first;
            const float p_ik=ii->second;
            FOREACH(SV::const_iterator, jj, mp_[y][x][l])
            {
              const uint j=jj->first;
              const float p_jl=jj->second;
              if (i<j) p[i][j] += p_kl*p_ik*p_jl*w;
            }
          }
        }
      }
    }
    
    bp[x].resize(L1);
    for (uint i=0; i!=L1-1; ++i)
      for (uint j=i+1; j!=L1; ++j)
      {
        float v = p[i][j]/sum_w;
        if (v>CUTOFF)
          bp[x][i].push_back(std::make_pair(j, v));
      }
  }
  std::swap(bp_, bp);
}

void
Dusaf::
relax_fourway_consistency()
{
  const uint N=fa_.size();
  std::vector<std::vector<MP> > mp(N, std::vector<MP>(N));
  assert(mp_.size()==N);
  assert(mp_[0].size()==N);
  for (uint x=0; x!=N-1; ++x)
  {
    const uint L1=fa_[x].size();
    for (uint y=x+1; y!=N; ++y)
    {
      const uint L2=fa_[y].size();
      VVF posterior(L1, VF(L2, 0.0));
      assert(L1==mp_[x][y].size());

      for (uint i=0; i!=L1; ++i)
      {
        FOREACH (SV::const_iterator, kk, mp_[x][y][i])
        {
          const uint k=kk->first;
          const float p_ik=kk->second;
          posterior[i][k] += p_ik * (1.0-w_pct_f_);
          FOREACH (SV::const_iterator, jj, bp_[x][i])
          {
            const uint j=jj->first;
            const float p_ij=jj->second;
            SV::const_iterator ll1=mp_[x][y][j].begin();
            SV::const_iterator ll2=bp_[y][k].begin();
            while (ll1!=mp_[x][y][j].end() && ll2!=bp_[y][k].end())
            {
              if (ll1->first<ll2->first) ++ll1;
              else if (ll1->first>ll2->first) ++ll2;
              else /* if (ll1->first==ll2->first) */
              {
                const uint l=ll1->first;
                const float p_jl=ll1->second;
                const float p_kl=ll2->second;
                posterior[i][k] += p_ij * p_kl * p_jl * w_pct_f_;
                posterior[j][l] += p_ij * p_kl * p_ik * w_pct_f_;
                ++ll1; ++ll2;
              }
            }
          }
        }
      }
      
      mp[x][y].resize(L1);
      for (uint i=0; i!=L1; ++i)
      {
        for (uint j=0; j!=L2; ++j)
        {
          float v=posterior[i][j];
          if (v>CUTOFF)
            mp[x][y][i].push_back(std::make_pair(j, v));
        }
      }
      transpose_mp(mp[x][y], mp[y][x], L1, L2);
    }
  }

  for (uint x=0; x!=N; ++x)
  {
    mp[x][x].resize(fa_[x].size());
    for (uint i=0; i!=fa_[x].size(); ++i)
      mp[x][x][i].push_back(std::make_pair(i, 1.0f));
  }
  std::swap(mp_, mp);
}

void
Dusaf::
build_tree()
{
  uint n=fa_.size();
  tree_.resize(2*n-1);
  std::fill(tree_.begin(), tree_.end(), std::make_pair(0.0, std::make_pair(-1u, -1u)));
  
  VVF d(n, VF(n, 0.0));
  VU idx(2*n-1, -1u);
  for (uint i=0; i!=n; ++i) idx[i]=i;

  std::priority_queue<node_t> pq;
  for (uint i=0; i!=n-1; ++i)
  {
    for (uint j=i+1; j!=n; ++j)
    {
      d[i][j] = d[j][i] = sim_[i][j];
      pq.push(std::make_pair(sim_[i][j], std::make_pair(i,j)));
    }
  }

  while (!pq.empty())
  {
    node_t t=pq.top(); pq.pop();
    if (idx[t.second.first]!=-1u && idx[t.second.second]!=-1u)
    {
      assert(n<tree_.size());
      const uint l = idx[t.second.first];
      const uint r = idx[t.second.second];
      idx[t.second.first] = idx[t.second.second] = -1u;
      for (uint i=0; i!=n; ++i)
      {
        if (idx[i]!=-1u)
        {
          uint ii=idx[i];
          d[ii][l] = d[l][ii] = (d[ii][l]+d[ii][r])*t.first/2;
          pq.push(std::make_pair(d[ii][l], std::make_pair(i,n)));
        }
      }
      tree_[n] = t;
      idx[n++] = l;
    }
  }
  assert(n==tree_.size());
}

// print the guide tree
void
Dusaf::
print_tree(std::ostream& os, int i) const
{
  if (tree_[i].second.first==-1u)
  {
    assert(tree_[i].second.second==-1u);
    os << fa_[i].name();
  }
  else
  {
    os << "[ " << tree_[i].first << " ";
    print_tree(os, tree_[i].second.first);
    os << " ";
    print_tree(os, tree_[i].second.second);
    os << " ]";
  }
}

void
Dusaf::
average_matching_probability(VVF& posterior, const ALN& aln1, const ALN& aln2) const
{
  const uint L1 = aln1.front().second.size();
  const uint L2 = aln2.front().second.size();
  const uint N1 = aln1.size();
  const uint N2 = aln2.size();
  VVF p(L1, VF(L2, 0.0));
  FOREACH (ALN::const_iterator, it1, aln1)
  {
    assert(L1==it1->second.size());
    FOREACH (ALN::const_iterator, it2, aln2)
    {
      assert(L2==it2->second.size());
      const MP& m = mp_[it1->first][it2->first];
      for (uint i=0, ii=0; i!=L1; ++i)
      {
        if (!it1->second[i]) continue;
        assert(ii<m.size());
        SV::const_iterator x=m[ii].begin();
        for (uint j=0, jj=0; j!=L2 && x!=m[ii].end(); ++j)
        {
          if (!it2->second[j]) continue;
          if (jj==x->first)
          {
            p[i][j] += x->second/(N1*N2);
            ++x;
          }
          ++jj;
        }
        ++ii;
      }
    }
  }
  for (uint i=0; i!=p.size(); ++i)
    for (uint j=0; j!=p[i].size(); ++j)
    {
      if (p[i][j]<=CUTOFF) p[i][j]=0.0;
      assert(p[i][j]<=1.0);
    }
  std::swap(posterior, p);
}

void
Dusaf::
average_basepairing_probability(VVF& posterior, const ALN& aln) const
{
  const uint L=aln.front().second.size();
  const uint N=aln.size();
  VVF p(L, VF(L, 0.0));
  FOREACH (ALN::const_iterator, it, aln)
  {
    assert(L==it->second.size());
    uint s=it->first;
    VU idx(fa_[s].size());
    for (uint i=0, j=0; i!=L; ++i)
      if (it->second[i]) idx[j++]=i;
    const BP& bp=bp_[s];
    for (uint i=0; i!=bp.size(); ++i)
      for (uint j=0; j!=bp[i].size(); ++j)
        p[idx[i]][idx[bp[i][j].first]] += bp[i][j].second/N;
  }

  if (use_alifold_)
  {
    char** seqs = new char*[N+1];
    seqs[N] = NULL;
    char* str = new char[L+1];
    char** s = seqs;
    FOREACH (ALN::const_iterator, it, aln)
    {
      *s = new char[L+1];
      for (uint i=0, j=0; i!=L; ++i)
        (*s)[i] = it->second[i] ? fa_[it->first].seq()[j++] : '-';
      (*s)[L] = 0;
      ++s;
    }

    // scaling parameters to avoid overflow
    double min_en = Vienna::alifold(seqs, str);
    double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
    Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
    Vienna::plist* pi;
#else
    Vienna::pair_info* pi;
#endif
    Vienna::alipf_fold(seqs, str, &pi);
    for (uint k=0; pi[k].i!=0; ++k)
      p[pi[k].i-1][pi[k].j-1] += pi[k].p;
    free(pi);

    Vienna::free_alipf_arrays();
    delete[] str;
    for (s=seqs; *s!=NULL; ++s)
      delete[] *s;
    delete[] seqs;

    for (uint i=0; i!=L-1; ++i)
      for (uint j=i+1; j!=L; ++j)
        p[i][j] /= 2;
  }

  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
    {
      if (p[i][j]<=CUTOFF) p[i][j]=0.0;
      assert(p[i][j]<=1.0);
    }
  std::swap(p, posterior);
}

static
float
calculate_similarity_score(const MP& mp, uint L1, uint L2)
{
  assert(mp.size()==L1);

  VVF dp(L1+1, VF(L2+1, 0.0));
  VVI tr(L1+1, VI(L2+1, 0));
  for (uint i=1; i!=L1+1; ++i)
  {
    uint j=1;
    FOREACH (SV::const_iterator, jj, mp[i-1])
    {
      for (; j-1<jj->first; ++j)
      {
        dp[i][j] = dp[i][j-1];
        tr[i][j] = tr[i][j-1]+1;
        if (dp[i][j]<dp[i-1][j])
        {
          dp[i][j] = dp[i-1][j];
          tr[i][j] = tr[i-1][j]+1;
        }
      }

      dp[i][j] = dp[i-1][j-1]+jj->second;
      tr[i][j] = tr[i-1][j-1]+1;
      if (dp[i][j]<dp[i][j-1])
      {
        dp[i][j] = dp[i][j-1];
        tr[i][j] = tr[i][j-1]+1;
      }
      if (dp[i][j]<dp[i-1][j])
      {
        dp[i][j] = dp[i-1][j];
        tr[i][j] = tr[i-1][j]+1;
      }
      ++j;
    }

    for (; j<L2+1; ++j)
    {
      dp[i][j] = dp[i][j-1];
      tr[i][j] = tr[i][j-1]+1;
      if (dp[i][j]<dp[i-1][j])
      {
        dp[i][j] = dp[i-1][j];
        tr[i][j] = tr[i-1][j]+1;
      }
    }
  }
  //return dp[L1][L2]/(std::min(L1, L2)); // simple version
  return dp[L1][L2]/tr[L1][L2];
}

void
Dusaf::
project_alignment(ALN& aln, const ALN& aln1, const ALN& aln2, const VU& z) const
{
  const uint L1=aln1[0].second.size();
  const uint L2=aln2[0].second.size();
  aln.resize(aln1.size()+aln2.size());
  uint c=0;
  for (uint i=0; i!=z.size(); ++i) if (z[i]!=-1u) c++;
  const uint L=L1+L2-c;
  ALN::iterator p=aln.begin();
  FOREACH (ALN::const_iterator, q, aln1)
  {
    p->first=q->first;
    p->second.resize(L, false);
    uint r=0, k=0;
    for (uint i=0; i!=z.size(); ++i)
    {
      if (z[i]!=-1u)
      {
        while (k<z[i]) { p->second[r++]=false; k++; }
        p->second[r++] = q->second[i]; ++k;
      }
      else
        p->second[r++] = q->second[i];
    }
    while (k<L2) { p->second[r++] = false; k++; }
    ++p;
  }
  FOREACH (ALN::const_iterator, q, aln2)
  {
    p->first=q->first;
    p->second.resize(L, false);
    uint k=0, r=0;
    for (uint i=0; i!=z.size(); ++i)
    {
      if (z[i]!=-1u)
      {
        while (k<z[i]) p->second[r++] = q->second[k++];
        p->second[r++] = q->second[k++];
      }
      else
        p->second[r++] = false;
    }
    while (k<L2) p->second[r++] = q->second[k++];
    ++p;
  }
}

void
Dusaf::
project_secondary_structure(VU& xx, VU& yy, const VU& x, const VU& y, const VU& z) const
{
  const uint L1=x.size();
  const uint L2=y.size();
  VU idx1(L1, -1u), idx2(L2, -1u);
  uint r=0, k=0;
  for (uint i=0; i!=z.size(); ++i)
  {
    if (z[i]!=-1u)
    {
      while (k<z[i])
      {
        idx2[k] = r;
        ++r; ++k;
      }
      idx1[i] = r;
      idx2[k] = r;
      ++r; ++k;
    }
    else
    {
      idx1[i] = r;
      ++r;
    }
  }
  while (k<L2)
  {
    idx2[k] = r;
    ++r; ++k;
  }
  const uint L=r;

  xx.resize(L); std::fill(xx.begin(), xx.end(), -1u);
  yy.resize(L); std::fill(yy.begin(), yy.end(), -1u);
  for (uint i=0; i!=L1; ++i)
    if (x[i]!=-1u)
      xx[idx1[i]]=idx1[x[i]];
  for (uint k=0; k!=L2; ++k)
    if (y[k]!=-1u)
      yy[idx2[k]]=idx2[y[k]];
}

void
Dusaf::
output_verbose(const VU& x, const VU& y, const VU& z, const ALN& aln1, const ALN& aln2) const
{
  ALN aln;
  project_alignment(aln, aln1, aln2, z);
  VU xx, yy;
  project_secondary_structure(xx, yy, x, y, z);
  const uint L=aln[0].second.size();

  std::string x_str(L, '.'), y_str(L, '.');
  for (uint i=0; i!=L; ++i)
    if (xx[i]!=-1u) { x_str[i]='('; x_str[xx[i]]=')'; }
  for (uint k=0; k!=L; ++k)
    if (yy[k]!=-1u) { y_str[k]='('; y_str[yy[k]]=')'; }
  
  output(std::cout, aln.begin(), aln.begin()+aln1.size());
  std::cout << x_str << std::endl;

  output(std::cout, aln.begin()+aln1.size(), aln.end());
  std::cout << y_str << std::endl;

  std::cout << std::endl;
}

void
Dusaf::
align_alignments(ALN& aln, const ALN& aln1, const ALN& aln2) const
{
  // calculate posteriors
  VVF p_x, p_y, p_z;
  average_basepairing_probability(p_x, aln1);
  average_basepairing_probability(p_y, aln2);
  average_matching_probability(p_z, aln1, aln2);

  // solve the problem
  VU x, y, z;
  solve(x, y, z, p_x, p_y, p_z, aln1, aln2);

  // build the result
  project_alignment(aln, aln1, aln2, z);
}

float
Dusaf::
align_alignments(VU& ss, ALN& aln, const ALN& aln1, const ALN& aln2) const
{
  // calculate posteriors
  VVF p_x, p_y, p_z;
  average_basepairing_probability(p_x, aln1);
  average_basepairing_probability(p_y, aln2);
  average_matching_probability(p_z, aln1, aln2);

  // solve the problem
  VU x, y, z;
  solve(x, y, z, p_x, p_y, p_z, aln1, aln2);

  // build the result
  project_alignment(aln, aln1, aln2, z);
  VU xx, yy;
  project_secondary_structure(xx, yy, x, y, z);
  assert(xx.size()==yy.size());
  ss.resize(xx.size());
  std::fill(ss.begin(), ss.end(), -1u);
  for (uint i=0; i!=ss.size(); ++i)
    if (xx[i]==yy[i]) ss[i]=xx[i];

  // calculate alignment score
  float a_score=0.0;
  for (uint i=0; i!=z.size(); ++i)
    if (z[i]!=-1u) a_score += p_z[i][z[i]]-th_a_;

  // calculate folding score
  assert(x.size()==z.size());
  float s_score=0.0;
  for (uint i=0; i!=x.size(); ++i)
  {
    if (x[i]!=-1u && z[i]!=-1u)
    {
      const uint j=x[i];
      const uint k=z[i];
      if (z[j]!=-1u && z[j]==y[k])
      {
        const uint l=y[k];
        s_score += p_x[i][j]-th_s_;
        s_score += p_y[k][l]-th_s_;
      }
    }
  }

  return w_*s_score + a_score;
}

void
Dusaf::
solve_by_dd(VU& x, VU& y, VU& z,
            const VVF& p_x, const VVF& p_y, const VVF& p_z,
            const ALN& aln1, const ALN& aln2) const
{
  const uint L1=p_x.size();
  const uint L2=p_y.size();

  // enumerate the candidates of consensus base-pairs
  std::vector<CBP> cbp;         // consensus base-pairs
#ifdef SPARSE_UPDATE
  VVU c_x(L1), c_y(L2), c_z(L1); // project consensus base-pairs into each structure and alignment
#endif
  for (uint i=0; i!=L1-1; ++i)
    for (uint j=i+1; j!=L1; ++j)
      if (p_x[i][j]>CUTOFF)
        for (uint k=0; k!=L2-1; ++k)
          if (p_z[i][k]>CUTOFF)
            for (uint l=k+1; l!=L2; ++l)
              if (p_y[k][l]>CUTOFF && p_z[j][l]>CUTOFF)
              {
                assert(p_x[i][j]<=1.0);
                assert(p_y[k][l]<=1.0);
                float p=(p_x[i][j]+p_y[k][l])/2;
                float q=(p_z[i][k]+p_z[j][l])/2;
                if (p-th_s_>0.0 && w_*(p-th_s_)+(q-th_a_)>0.0)
                {
                  cbp.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
#ifdef SPARSE_UPDATE
                  c_x[i].push_back(j);
                  c_y[k].push_back(l);
                  c_z[i].push_back(k);
                  c_z[j].push_back(l);
#endif
                }
              }
#ifdef SPARSE_UPDATE
  for (uint i=0; i!=c_x.size(); ++i)
  {
    std::sort(c_x[i].begin(), c_x[i].end());
    c_x[i].erase(std::unique(c_x[i].begin(), c_x[i].end()), c_x[i].end());
  }
  for (uint k=0; k!=c_y.size(); ++k)
  {
    std::sort(c_y[k].begin(), c_y[k].end());
    c_y[k].erase(std::unique(c_y[k].begin(), c_y[k].end()), c_y[k].end());
  }
  for (uint i=0; i!=c_z.size(); ++i)
  {
    std::sort(c_z[i].begin(), c_z[i].end());
    c_z[i].erase(std::unique(c_z[i].begin(), c_z[i].end()), c_z[i].end());
  }
#endif

  // precalculate the range for alignment, i.e. alignment envelope
  a_decoder_->initialize(p_z);

  // multipliers
  VVF q_x(L1, VF(L1, 0.0));
  VVF q_y(L2, VF(L2, 0.0));
  VVF q_z(L1, VF(L2, 0.0));

  //uint c=0;
  float c=0.0;
  float eta=eta0_, s_prev=0.0;
  uint violated=0;
  uint t;
  for (t=0; t!=t_max_; ++t)
  {
    // solve the subproblems
    float s = 0.0;
    s += s_decoder_->decode(p_x, q_x, x);
    s += s_decoder_->decode(p_y, q_y, y);
    s += a_decoder_->decode(p_z, q_z, z);

    if (verbose_>=2) output_verbose(x, y, z, aln1, aln2);

    // update the multipliers
    violated=0;
    VVI t_x(L1, VI(L1, 0));
    VVI t_y(L2, VI(L2, 0));
    VVI t_z(L1, VI(L2, 0));
    for (uint u=0; u!=cbp.size(); ++u)
    {
      const uint i=cbp[u].first.first, j=cbp[u].first.second;
      const uint k=cbp[u].second.first, l=cbp[u].second.second;
      const float s_w = q_x[i][j]+q_y[k][l]-q_z[i][k]-q_z[j][l];
      const int w_ijkl = s_w>0.0f ? 1 : 0;
      if (w_ijkl)
      {
        s += s_w * w_ijkl;
        t_x[i][j] += w_ijkl;
        t_y[k][l] += w_ijkl;
        t_z[i][k] += w_ijkl;
        t_z[j][l] += w_ijkl;
      }
    }

#ifdef SPARSE_UPDATE // efficient implementation using sparsity
    for (uint i=0; i!=L1; ++i)
    {
      const uint j = x[i];
      if (j!=-1u && t_x[i][j]!=1)
      {
        violated++;
        q_x[i][j] -= eta*(t_x[i][j]-1);
      }
      for (uint jj=0; jj!=c_x[i].size(); ++jj)
      {
        const uint j=c_x[i][jj];
        if (x[i]!=j && t_x[i][j]!=0)
        {
          violated++;
          q_x[i][j] -= eta*t_x[i][j];
        }
      }
    }

    for (uint k=0; k!=L2; ++k)
    {
      const uint l = y[k];
      if (l!=-1u && t_y[k][l]!=1)
      {
        violated++;
        q_y[k][l] -= eta*(t_y[k][l]-1);
      }
      for (uint ll=0; ll!=c_y[k].size(); ++ll)
      {
        const uint l=c_y[k][ll];
        if (y[k]!=l && t_y[k][l]!=0)
        {
          violated++;
          q_y[k][l] -= eta*t_y[k][l];
        }
      }
    }

    for (uint i=0; i!=L1; ++i)
    {
      const uint k = z[i];
      if (k!=-1u && t_z[i][k]>1)
      {
        violated++;
        q_z[i][k] = std::max(0.0f, q_z[i][k]-eta*(1-t_z[i][k]));
      }
      for (uint kk=0; kk!=c_z[i].size(); ++kk)
      {
        const uint k=c_z[i][kk];
        if (z[i]!=k && t_z[i][k]>0)
        {
          violated++;
          q_z[i][k] = std::max(0.0f, q_z[i][k]+eta*t_z[i][k]);
        }
      }
    }
#else // naive implementation
    for (uint i=0; i!=L1-1; ++i)
      for (uint j=i+1; j!=L1; ++j)
      {
        const int x_ij = x[i]==j ? 1 : 0;
        if (t_x[i][j]-x_ij!=0)
        {
          violated++;
          q_x[i][j] -= eta*(t_x[i][j]-x_ij);
        }
      }
    for (uint k=0; k!=L2-1; ++k)
      for (uint l=k+1; l!=L2; ++l)
      {
        const int y_kl = y[k]==l ? 1 : 0;
        if (t_y[k][l]-y_kl!=0)
        {
          violated++;
          q_y[k][l] -= eta*(t_y[k][l]-y_kl);
        }
      }
    for (uint i=0; i!=L1; ++i)
      for (uint k=0; k!=L2; ++k)
      {
        const int z_ik = z[i]==k ? 1 : 0;
        if (z_ik-t_z[i][k]<0)
        {
          violated++;
          q_z[i][k] = std::max(0.0f, q_z[i][k]-eta*(z_ik-t_z[i][k]));
        }
      }
#endif

    if (verbose_>=2)
      std::cout << "Step: " << t << std::endl
                << "eta: " << eta << std::endl
                << "L: " << s << std::endl
                << "Violated: " << violated << std::endl << std::endl;

    if (violated==0)  break;     // all constraints were satisfied.

    // update the step width
    if (s<s_prev || t==0)
    {
      c += std::max(0.0f, 4.0f*cbp.size()-violated)/(4.0*cbp.size());
      eta = eta0_/(1.0+sqrt(c));
    }
    s_prev = s;
  }
  if (verbose_==1) std::cout << "Step: " << t << ", Violated: " << violated << std::endl;
}

void
Dusaf::
solve_by_ip(VU& x, VU& y, VU& z,
            const VVF& p_x, const VVF& p_y, const VVF& p_z,
            const ALN& aln1, const ALN& aln2) const
{
  const uint L1=p_x.size();
  const uint L2=p_y.size();

  // integer programming
  IP ip(IP::MAX, 1);

  // variables
  VVI v_x(L1, VI(L1, -1));
  VVI v_y(L2, VI(L2, -1));
  VVI v_z(L1, VI(L2, -1));
  VI v_w;

  // enumerate the candidates of aligned bases
  for (uint i=0; i!=L1; ++i)
    for (uint k=0; k!=L2; ++k)
      if (p_z[i][k]>CUTOFF)
        v_z[i][k] = ip.make_variable(p_z[i][k]-th_a_);
  
  // enumerate the candidates of consensus base-pairs
  std::vector<CBP> cbp;
  for (uint i=0; i!=L1-1; ++i)
    for (uint j=i+1; j!=L1; ++j)
      if (p_x[i][j]>CUTOFF)
        for (uint k=0; k!=L2-1; ++k)
          if (p_z[i][k]>CUTOFF)
            for (uint l=k+1; l!=L2; ++l)
              if (p_y[k][l]>CUTOFF && p_z[j][l]>CUTOFF)
              {
                assert(p_x[i][j]<=1.0);
                assert(p_y[k][l]<=1.0);
                float p=(p_x[i][j]+p_y[k][l])/2;
                float q=(p_z[i][k]+p_z[j][l])/2;
                if (p-th_s_>0.0 && w_*(p-th_s_)+(q-th_a_)>0.0)
                {
                  cbp.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
                  v_w.push_back(ip.make_variable(0.0));
                  if (v_x[i][j]<0)
                    v_x[i][j] = ip.make_variable(w_*(p_x[i][j]-th_s_));
                  if (v_y[k][l]<0)
                    v_y[k][l] = ip.make_variable(w_*(p_y[k][l]-th_s_));
                }
              }
  ip.update();

  // constraints: each base is paired with at most one base (in a)
  for (uint i=0; i<L1; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint j=0; j<i; ++j)
      if (v_x[j][i]>=0)
        ip.add_constraint(row, v_x[j][i], 1);
    for (uint j=i+1; j<L1; ++j)
      if (v_x[i][j]>=0)
        ip.add_constraint(row, v_x[i][j], 1);
  }
  
  // constraints: no pseudoknots are allowed (in a)
  for (uint i=0; i<L1-1; ++i)
    for (uint j=i+1; j<L1; ++j)
      if (v_x[i][j]>=0)
        for (uint k=i+1; k<j; ++k)
          for (uint l=j+1; l<L1; ++l)
            if (v_x[k][l]>=0)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, v_x[i][j], 1);
              ip.add_constraint(row, v_x[k][l], 1);
            }

  // constraints: each base is paired with at most one base (in b)
  for (uint i=0; i<L2; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint j=0; j<i; ++j)
      if (v_y[j][i]>=0)
        ip.add_constraint(row, v_y[j][i], 1);
    for (uint j=i+1; j<L2; ++j)
      if (v_y[i][j]>=0)
        ip.add_constraint(row, v_y[i][j], 1);
  }

  // constraints: no pseudoknots are allowed (in b)
  for (uint i=0; i<L2-1; ++i)
    for (uint j=i+1; j<L2; ++j)
      if (v_y[i][j]>=0)
        for (uint k=i+1; k<j; ++k)
          for (uint l=j+1; l<L2; ++l)
            if (v_y[k][l]>=0)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, v_y[i][j], 1);
              ip.add_constraint(row, v_y[k][l], 1);
            }

  // constraints: each base is aligned with at most one base
  for (uint i=0; i<L1; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint k=0; k<L2; ++k)
      if (v_z[i][k]>=0)
        ip.add_constraint(row, v_z[i][k], 1);
  }
  for (uint k=0; k<L2; ++k)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint i=0; i<L1; ++i)
      if (v_z[i][k]>=0)
        ip.add_constraint(row, v_z[i][k], 1);
  }

  // constraints: no crossing matches are allowed
  for (uint i=0; i<L1; ++i)
    for (uint k=0; k<L2; ++k)
      if (v_z[i][k]>=0)
        for (uint j=i+1; j<L1; ++j)
          for (uint l=0; l<k; ++l)
            if (v_z[j][l]>=0)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, v_z[i][k], 1);
              ip.add_constraint(row, v_z[j][l], 1);
            }

  // constraints for consensus base pairs
  VVI r_x(L1, VI(L1, -1));
  for (uint i=0; i<L1-1; ++i)
    for (uint j=i+1; j<L1; ++j)
      if (v_x[i][j]>=0)
      {
        r_x[i][j] = ip.make_constraint(IP::FX, 0, 0);
        ip.add_constraint(r_x[i][j], v_x[i][j], 1);
      }

  VVI r_y(L2, VI(L2, -1));
  for (uint i=0; i<L2-1; ++i)
    for (uint j=i+1; j<L2; ++j)
      if (v_y[i][j]>=0)
      {        
        r_y[i][j] = ip.make_constraint(IP::FX, 0, 0);
        ip.add_constraint(r_y[i][j], v_y[i][j], 1);
      }

  VVI r_z(L1, VI(L2, -1));
  for (uint i=0; i<L1; ++i)
    for (uint k=0; k<L2; ++k)
      if (v_z[i][k]>=0)
      {
        r_z[i][k] = ip.make_constraint(IP::LO, 0, 0);
        ip.add_constraint(r_z[i][k], v_z[i][k], 1);
      }

  for (uint u=0; u!=cbp.size(); ++u)
  {
    const uint i=cbp[u].first.first, j=cbp[u].first.second;
    const uint k=cbp[u].second.first, l=cbp[u].second.second;
    assert(r_x[i][j]>=0 && v_x[i][j]>=0);
    ip.add_constraint(r_x[i][j], v_w[u], -1);
    assert(r_y[k][l]>=0 && v_y[k][l]>=0);
    ip.add_constraint(r_y[k][l], v_w[u], -1);
    assert(r_z[i][k]>=0 && v_z[i][k]>=0);
    ip.add_constraint(r_z[i][k], v_w[u], -1);
    assert(r_z[j][l]>=0 && v_z[j][l]>=0);
    ip.add_constraint(r_z[j][l], v_w[u], -1);
  }

  // execute optimization
  ip.solve();
  
  // build the result
  x.resize(L1);
  std::fill(x.begin(), x.end(), -1u);
  for (uint i=0; i<L1-1; ++i)
    for (uint j=i+1; j<L1; ++j)
      if (v_x[i][j]>=0 && ip.get_value(v_x[i][j])>0.5)
        x[i]=j;

  y.resize(L2);
  std::fill(y.begin(), y.end(), -1u);
  for (uint i=0; i<L2-1; ++i)
    for (uint j=i+1; j<L2; ++j)
      if (v_y[i][j]>=0 && ip.get_value(v_y[i][j])>0.5)
        y[i]=j;

  z.resize(L1);
  std::fill(z.begin(), z.end(), -1u);
  for (uint i=0; i<L1; ++i)
    for (uint k=0; k<L2; ++k)
      if (v_z[i][k]>=0 && ip.get_value(v_z[i][k])>0.5)
        z[i]=k;
}

void
Dusaf::
align(ALN& aln, int ch) const
{
  if (tree_[ch].second.first==-1u)
  {
    assert(tree_[ch].second.second==-1u);
    aln.resize(1);
    aln[0].first = ch;
    aln[0].second = std::vector<bool>(fa_[ch].size(), true);
  }
  else
  {
    ALN aln1, aln2;
    align(aln1, tree_[ch].second.first);
    align(aln2, tree_[ch].second.second);
    align_alignments(aln, aln1, aln2);
  }
}

float
Dusaf::
align(VU& ss, ALN& aln, int ch) const
{
  float s=0.0;
  if (tree_[ch].second.first==-1u)
  {
    assert(tree_[ch].second.second==-1u);
    aln.resize(1);
    aln[0].first = ch;
    aln[0].second = std::vector<bool>(fa_[ch].size(), true);
  }
  else
  {
    ALN aln1, aln2;
    align(aln1, tree_[ch].second.first);
    align(aln2, tree_[ch].second.second);
    s=align_alignments(ss, aln, aln1, aln2);
  }
  return s;
}

float
Dusaf::
refine(VU& ss, ALN& aln) const
{ 
 VU group[2];
  do {
    group[0].clear();
    group[1].clear();
    for (uint i=0; i!=aln.size(); ++i)
      group[rand()%2].push_back(i);
  } while (group[0].empty() || group[1].empty());

  ALN a[2];
  for (uint i=0; i!=2; ++i)
  {
    uint N=group[i].size();
    a[i].resize(N);
    uint L=aln[group[i][0]].second.size();
    for (uint j=0; j!=N; ++j)
      a[i][j].first = aln[group[i][j]].first;
    for (uint k=0; k!=L; ++k)
    {
      bool gap=true;
      for (uint j=0; j!=N; ++j)
        gap &= !aln[group[i][j]].second[k];
      if (!gap)
      {
        for (uint j=0; j!=N; ++j)
          a[i][j].second.push_back(aln[group[i][j]].second[k]);
      }
    }
  }

  ALN r;
  float s=align_alignments(ss, r, a[0], a[1]);
  std::swap(r, aln);
  return s;
}

void
Dusaf::
output(std::ostream& os, const ALN& aln, const VU& ss) const
{
  os << ">SS_cons" << std::endl;
  std::string s(ss.size(), '.');
  for (uint i=0; i!=ss.size(); ++i)
    if (ss[i]!=-1u)
    {
      s[i]='(';
      s[ss[i]]=')';
    }
  std::cout << s << std::endl;
  output(os, aln);
}

void
Dusaf::
output(std::ostream& os, const ALN& aln) const
{
  output(os, aln.begin(), aln.end());
}

void
Dusaf::
output(std::ostream& os, ALN::const_iterator b, ALN::const_iterator e) const
{
  for (ALN::const_iterator a=b; a!=e; ++a)
  {
    uint s=a->first;
    os << ">" << " " << fa_[s].name() << std::endl;
    for (uint j=0, k=0; j!=a->second.size(); ++j)
    {
      if (a->second[j])
        os << fa_[s].seq()[k++];
      else
        os << '-';
    }
    os << std::endl;
  }
}

Dusaf&
Dusaf::
parse_options(int& argc, char**& argv)
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info)!=0) exit(1);

  // general options
  n_refinement_ = args_info.refinement_arg;
  w_ = args_info.weight_arg;
  eta0_ = args_info.eta_arg;
  t_max_ = args_info.max_iter_arg;
  w_pct_f_ = args_info.fourway_pct_arg;
  verbose_ = args_info.verbose_arg;

  // options for alignments
  th_a_ = args_info.align_th_arg;
  w_pct_a_ = args_info.align_pct_arg;
  // use_bpscore_ = args_info.use_bpscore_flag!=0;
  // std::string arg_x;
  // if (args_info.extra_given) arg_x = std::string(args_info.extra_arg);
  if (strcasecmp(args_info.align_model_arg, "CONTRAlign")==0)
    a_model_ = new CONTRAlign(th_a_);
  else if (strcasecmp(args_info.align_model_arg, "ProbCons")==0)
    a_model_ = new ProbCons(th_a_);
#if 0
  else if (strcasecmp(args_info.align_model_arg, "PartAlign")==0)
    a_model_ = new PartAlign(th_a_, arg_x);
#endif
  assert(a_model_!=NULL);
  a_decoder_ = new SparseNeedlemanWunsch(th_a_);

  // options for folding
  w_pct_s_ = args_info.fold_pct_arg;
  th_s_ = args_info.fold_th_arg[0];
  if (args_info.gamma_given) th_s_ = 1.0/(1.0+args_info.gamma_arg[0]);
  use_alifold_ = args_info.use_alifold_flag!=0;
  if (strcasecmp(args_info.fold_model_arg, "Boltzmann")==0)
    s_model_ = new RNAfold(true, NULL, CUTOFF);
  else if (strcasecmp(args_info.fold_model_arg, "Vienna")==0)
    s_model_ = new RNAfold(false, NULL, CUTOFF);
  else if (strcasecmp(args_info.fold_model_arg, "CONTRAfold")==0)
    s_model_ = new CONTRAfold(CUTOFF);
  assert(s_model_!=NULL);

  s_decoder_ = new SparseNussinov(w_, th_s_);
  if (args_info.fold_th1_given)
    s_decoder1_ = new SparseNussinov(w_, args_info.fold_th1_arg[0]);
  else if (args_info.gamma1_given)
    s_decoder1_ = new SparseNussinov(w_, 1.0/(1.0+args_info.gamma1_arg[0]));
  else
    s_decoder1_ = new SparseNussinov(w_, th_s_);

  if (args_info.inputs_num==0)
  {
    cmdline_parser_print_help();
    cmdline_parser_free(&args_info);
    exit(1);
  }

  // read sequences
  Fasta::load(fa_, args_info.inputs[0]);

  cmdline_parser_free(&args_info);
  return *this;
}

int
Dusaf::
run()
{
  const uint N=fa_.size();
  
  // calculate base-pairing probabilities
  bp_.resize(N);
  for (uint i=0; i!=N; ++i)
    s_model_->calculate(fa_[i].seq(), bp_[i]);

  // calculate matching probabilities
  mp_.resize(N, std::vector<MP>(N));
  for (uint i=0; i!=N; ++i)
  {
    mp_[i][i].resize(fa_[i].size());
    for (uint x=0; x!=fa_[i].size(); ++x)
      mp_[i][i][x].push_back(std::make_pair(x, 1.0f));
    for (uint j=i+1; j!=N; ++j)
    {
      // if (use_bpscore_)
      //   a_model_->calculate(fa_[i].seq(), fa_[j].seq(), bp_[i], bp_[j], mp_[i][j]);
      // else
      //   a_model_->calculate(fa_[i].seq(), fa_[j].seq(), mp_[i][j]);
      a_model_->calculate(fa_[i].seq(), fa_[j].seq(), mp_[i][j]);
      transpose_mp(mp_[i][j], mp_[j][i], fa_[i].size(), fa_[j].size());
    }
  }

  // four-way probabilistic consistency tranformation
  if (w_pct_f_!=0.0)
    relax_fourway_consistency();

  // calculate probabilistic similarity scores
  // which are used for building guide trees and PCTs
  sim_.resize(N, VF(N));
  for (uint i=0; i!=N; ++i)
  {
    sim_[i][i] = 1.0;
    for (uint j=i+1; j!=N; ++j)
      sim_[i][j] = sim_[j][i] = calculate_similarity_score(mp_[i][j], fa_[i].size(), fa_[j].size());
  }

  // probabilistic consistency tranformation for base-pairing probabilitiy matrix
  if (w_pct_s_!=0.0)
    relax_basepairing_probability();

  // probabilistic consistency tranformation for matching probability matrix
  if (w_pct_a_!=0.0)
    relax_matching_probability();
  
  // compute the guide tree
  build_tree();
  print_tree(std::cout, tree_.size()-1);
  std::cout << std::endl;

  // compute progressive alignments along with the guide tree
  VU ss;
  ALN aln;
  float s;
  s=align(ss, aln, tree_.size()-1);
  
  // iterative refinement
  for (uint i=0; i!=n_refinement_; ++i)
  {
    VU ss_temp=ss;
    ALN aln_temp=aln;
    float s_temp;

    s_temp=refine(ss_temp, aln_temp);
    //std::cout << s << " " << s_temp << std::endl;
    if (s_temp>s)
    {
      s=s_temp;
      std::swap(ss, ss_temp);
      std::swap(aln, aln_temp);
    }
  }
  if (s_decoder1_!=NULL)
  {
    // compute the common secondary structures from the averaged base-pairing matrix
    VVF p;
    average_basepairing_probability(p, aln);
    s_decoder1_->decode(p, ss);
  }
  
  // output the alignment
  std::sort(aln.begin(), aln.end());
  output(std::cout, aln, ss);

  return 0;
}

int
main(int argc, char* argv[])
{
  Dusaf dusaf;
  return dusaf.parse_options(argc, argv).run();
}
