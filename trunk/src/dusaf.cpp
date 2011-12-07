// $Id:$

#include "config.h"
#include <cstring>
#include <unistd.h>
#include <vector>
#include <queue>
#include <stack>
#include "fa.h"
#include "fold.h"
#include "align.h"
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
    : n_pct_a_(0),
      n_pct_s_(0),
      n_refinement_(0),
      t_max_(100),
      th_a_(1.0/(1+99)),
      th_s_(1.0/(1+4)),
      w_(1.0),
      eta0_(0.5),
      en_a_(NULL),
      en_s_(NULL),
      use_alifold_(false),
      use_bpscore_(false),
      verbose_(0)
  {
  }

  ~Dusaf()
  {
    if (en_a_) delete en_a_;
    if (en_s_) delete en_s_;
  }

  Dusaf& parse_options(int& argc, char**& argv);
  void usage(const char* progname) const { };
  int run(const char* filename);

private:
  void relax_matching_probability();
  void relax_basepairing_probability();
  void build_tree();
  void print_tree(std::ostream& os, int i) const;
  float decode_alignment(const VVF& p, const VVF& q,
                         const std::vector< std::pair<uint,uint> >& r, VU& al) const;
  float decode_secondary_structure(const VVF& p, const VVF& q, VU& ss) const;
  void project_alignment(ALN& aln, const ALN& aln1, const ALN& aln2, const VU& z) const;
  void project_secondary_structure(VU& xx, VU& yy, const VU& x, const VU& y, const VU& z) const;
  void average_matching_probability(VVF& posterior, const ALN& aln1, const ALN& aln2) const;
  void average_basepairing_probability(VVF& posterior, const ALN& aln) const;
  void align_alignments(ALN& aln, const ALN& aln1, const ALN& aln2) const;
  void align_alignments(VU& ss, ALN& aln, const ALN& aln1, const ALN& aln2) const;
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
  void align(VU& ss, ALN& aln, int ch) const;
  void refine(VU& ss, ALN& aln) const;
  void output_verbose(const VU& x, const VU& y, const VU& z, const ALN& aln1, const ALN& aln2) const;
  void output(std::ostream& os, const ALN& aln, const VU& ss) const;
  void output(std::ostream& os, const ALN& aln) const;
  void output(std::ostream& os, ALN::const_iterator b, ALN::const_iterator e) const;
  void print_basepairing_probability(std::ostream& os, const VVF& p) const;
  void print_matching_probability(std::ostream& os, const VVF& p) const;

private:
  uint n_pct_a_;                // the number of PCT for alignment matching probabilities
  uint n_pct_s_;                // the number of PCT for base-pairing probabilities
  uint n_refinement_;           // the number of the iterative refinement
  uint t_max_;                  // the maximum number of the iteration of the subgradient update
  float th_a_;                  // the threshold for base-pairing probabilities
  float th_s_;                  // the threshold for alignment matching probabilities
  float w_;                     // the weight for base pairs in the objective function
  float eta0_;                  // the initial step width of the subgradient update
  Align* en_a_;                 // alignment engine
  Fold* en_s_;                  // folding engine
  std::vector<Fasta> fa_;       // input sequences
  std::vector<std::vector<MP> > mp_; // alignment matching probability matrices
  std::vector<BP> bp_;          // base-pairing probability matrices
  std::vector<float> dist_;     // distance matrix of input sequences
  std::vector<node_t> tree_;    // guide tree
  bool use_alifold_;
  bool use_bpscore_;
  uint verbose_;
};

void
Dusaf::
relax_matching_probability()
{
  const uint N=fa_.size();
  std::vector<std::vector<MP> > mp(N, std::vector<MP>(N));
  assert(mp_.size()==N);
  assert(mp_[0].size()==N);
  for (uint x=0; x!=N; ++x)
  {
    const uint L1=fa_[x].size();
    for (uint y=0; y!=N; ++y)
    {
      const uint L2=fa_[y].size();
      VVF posterior(L1, VF(L2, 0.0));
      assert(L1==mp_[x][y].size());

      for (uint z=0; z!=N; ++z)
      {
        const uint L3=fa_[z].size();
        assert(L3==mp_[z][x].size());
        assert(L3==mp_[z][y].size());
        for (uint k=0; k!=L3; ++k)
        {
          FOREACH (SV::const_iterator, ii, mp_[z][x][k])
          {
            FOREACH (SV::const_iterator, jj, mp_[z][y][k])
            {
              const uint i=ii->first, j=jj->first;
              const float p_ik=ii->second, p_jk=jj->second;
              assert(i<L1); assert(j<L2);
              posterior[i][j] += p_ik * p_jk;
            }
          }
        }
      }
    
      mp[x][y].resize(L1);
      for (uint i=0; i!=L1; ++i)
      {
        for (uint j=0; j!=L2; ++j)
        {
          float v=posterior[i][j]/N;
          if (v>CUTOFF)
            mp[x][y][i].push_back(std::make_pair(j, v));
        }
      }
    }
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

    for (uint y=0; y!=N; ++y)
    {
      const uint L2=bp_[y].size();
      assert(L2==mp_[y][x].size());
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
              if (i<j) p[i][j] += p_kl*p_ik*p_jl/N;
            }
          }
        }
      }
    }
    
    bp[x].resize(L1);
    for (uint i=0; i!=L1-1; ++i)
      for (uint j=i+1; j!=L1; ++j)
        if (p[i][j]>CUTOFF)
          bp[x][i].push_back(std::make_pair(j, p[i][j]));
  }
  std::swap(bp_, bp);
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
  for (uint i=0, k=0; i!=n-1; ++i)
  {
    for (uint j=i+1; j!=n; ++j, ++k)
    {
      d[i][j] = d[j][i] = dist_[k];
      pq.push(std::make_pair(dist_[k], std::make_pair(i,j)));
    }
  }

  while (!pq.empty())
  {
    node_t t=pq.top(); pq.pop();
    if (idx[t.second.first]!=-1u && idx[t.second.second]!=-1u)
    {
      assert(n<tree_.size());
      uint l = idx[t.second.first];
      uint r = idx[t.second.second];
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

float
Dusaf::
decode_alignment(const VVF& p, const VVF& q,
                 const std::vector< std::pair<uint,uint> >& r, VU& al) const
{
  const uint L1 = p.size();
  const uint L2 = p[0].size();

  VVF dp(L1+1, VF(L2+1, 0.0));
  VVC tr(L1+1, VC(L2+1, ' '));
  for (uint i=1; i!=L1+1; ++i) tr[i][0] = 'X';
  for (uint k=1; k!=L2+1; ++k) tr[0][k] = 'Y';
#ifdef SPARSE_ALIGNMENT // efficient implementation using sparsity
  for (uint i=1; i!=L1+1; ++i)
  {
    for (uint k=r[i].first; k<=r[i].second; ++k)
    {
      if (k==0) continue;
      float v = dp[i-1][k-1]+p[i-1][k-1]-th_a_+q[i-1][k-1];
      char t = 'M';
      if (v<dp[i-1][k])
      {
        v = dp[i-1][k];
        t = 'X';
      }
      if (v<dp[i][k-1])
      {
        v = dp[i][k-1];
        t = 'Y';
      }
      dp[i][k] = v;
      tr[i][k] = t;
    }
  }
#else // naive implementation
  for (uint i=1; i!=L1+1; ++i)
  {
    for (uint k=1; k!=L2+1; ++k)
    {
      float v = dp[i-1][k-1]+p[i-1][k-1]-th_a_+q[i-1][k-1];
      char t = 'M';
      if (v<dp[i-1][k])
      {
        v = dp[i-1][k];
        t = 'X';
      }
      if (v<dp[i][k-1])
      {
        v = dp[i][k-1];
        t = 'Y';
      }
      dp[i][k] = v;
      tr[i][k] = t;
    }
  }
#endif

  // traceback
  std::string rpath;
  int i = L1, k = L2;
  while (i>0 || k>0)
  {
    rpath.push_back(tr[i][k]);
    switch (tr[i][k])
    {
      case 'M': --i; --k; break;
      case 'X': --i; break;
      case 'Y': --k; break;
      default: assert(!"unreachable"); break;
    }
  }
  std::string vpath(rpath.size(), ' ');
  std::reverse_copy(rpath.begin(), rpath.end(), vpath.begin());

  // decode the resultant alignment
  al.resize(L1, -1u);
  for (uint i=0, k=0, p=0; p!=vpath.size(); ++p)
  {
    switch (vpath[p])
    {
      case 'M':
        assert(i<L1); assert(k<L2);
        al[i++]=k++;
        break;
      case 'X':
        assert(i<L1); assert(k<=L2);
        al[i++]=-1u;
        break;
      case 'Y':
        assert(i<=L1); assert(k<L2);
        k++;
        break;
      default: break;
    }
  }

  return dp[L1][L2];
}

void
alignment_envelope(const VVF& p, float th, std::vector< std::pair<uint,uint> >& r)
{
#ifdef SPARSE_ALIGNMENT
  const uint L1 = p.size();
  const uint L2 = p[0].size();
  r.resize(L1+1);
  std::fill(r.begin(), r.end(), std::make_pair(0, 0));

  for (uint i=1; i!=L1+1; ++i)
  {
    // find the first alignable point
    for (uint k=1; k!=L2+1; ++k)
    {
      if (p[i-1][k-1]-th>=0.0)
      {
        r[i-1].first = std::min(r[i-1].first, k-1);
        r[i].first = k;
        break;
      }
    }
    // no alignable point
    if (r[i].first==0)          
    {
      r[i].first = r[i-1].first;
      r[i].second = r[i-1].second;
      continue;
    }

    // find the last alignable point
    for (uint k=L2; k!=0; --k)
    {
      if (p[i-1][k-1]-th>=0.0)
      {
        r[i-1].second = std::max(r[i-1].second, k-1);
        r[i].second = k;
        break;
      }
    }
  }
  assert(r[0].first==0);
  r[L1].second=L2;

  // force the envelope to be monotonic
  for (uint i=L1, v=L2; i!=0; --i)
    r[i].first = v = std::min(v, r[i].first);
  for (uint i=0, v=0; i!=L1+1; ++i)
    r[i].second = v = std::max(v, r[i].second);

  // ensure the conectivity
  for (uint i=1; i!=L1+1; ++i)
  {
    if (r[i-1].second<r[i].first)
      r[i].first = r[i-1].second;
  }
#endif
}

float
Dusaf::
decode_secondary_structure(const VVF& p, const VVF& q, VU& ss) const
{
  uint L=p.size();
  assert(p[0].size()==L);

#ifdef SPARSE_FOLDING // efficient implementation using sparsity
  VVF dp(L, VF(L, 0.0));
  BP bp(L);
  VVU tr(L, VU(L, 0));
  for (uint l=1; l<L; ++l)
  {
    for (uint i=0; i+l<L; ++i)
    {
      uint j=i+l;
      float v=0.0;
      int t=0;
      if (i+1<j)
      {
        v=dp[i+1][j];
        t=1;
      }
      if (i<j-1 && v<dp[i][j-1])
      {
        v=dp[i][j-1];
        t=2;
      }
      if (i+1<j-1)
      {
        float s=w_*(p[i][j]-th_s_)-q[i][j];
        if (s>0.0)
        {
          bp[j].push_back(std::make_pair(i,dp[i+1][j-1]+s));
          if (v<dp[i+1][j-1]+s)
          {
            v=dp[i+1][j-1]+s;
            t=3;
          }
        }
      }
      for (SV::const_iterator x=bp[j].begin(); x!=bp[j].end(); ++x)
      {
        const uint k=x->first;
        const float s=x->second;
        if (i<k)
        {
          if (v<dp[i][k-1]+s)
          {
            v=dp[i][k-1]+s;
            t=k-i+3;
          }
        }
      }
      dp[i][j]=v;
      tr[i][j]=t;
    }
  }

  // trace back
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  std::stack<std::pair<uint,uint> > st;
  st.push(std::make_pair(0, L-1));
  while (!st.empty())
  {
    const std::pair<uint,uint> p=st.top(); st.pop();
    const int i=p.first, j=p.second;
    switch (tr[i][j])
    {
      case 0:
        break;
      case 1:
        st.push(std::make_pair(i+1, j));
        break;
      case 2:
        st.push(std::make_pair(i, j-1));
        break;
      case 3:
        ss[i]=j;
        st.push(std::make_pair(i+1, j-1));
        break;
      default:
        const int k=i+tr[i][j]-3;
        st.push(std::make_pair(i, k-1));
        ss[k]=j;
        st.push(std::make_pair(k+1, j-1));
        break;
    }
  }

  return dp[0][L-1];

#else // naive implementation
  // calculate scoring matrices for the current step
  VVF sm(L, VF(L, 0.0));
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      sm[i][j] = w_*(p[i][j]-th_s_)-q[i][j];

  VVF dp(L, VF(L, 0.0));
  VVU tr(L, VU(L, 0));
  for (uint l=1; l<L; ++l)
  {
    for (uint i=0; i+l<L; ++i)
    {
      uint j=i+l;
      float v=0.0;
      int t=0;
      if (i+1<j)
      {
        v=dp[i+1][j];
        t=1;
      }
      if (i<j-1 && v<dp[i][j-1])
      {
        v=dp[i][j-1];
        t=2;
      }
      if (i+1<j-1 && v<dp[i+1][j-1]+sm[i][j])
      {
        v=dp[i+1][j-1]+sm[i][j];
        t=3;
      }
      for (uint k=i+1; k<j; ++k)
      {
        if (v<dp[i][k]+dp[k+1][j])
        {
          v=dp[i][k]+dp[k+1][j];
          t=k-i+3;
        }        
      }
      dp[i][j]=v;
      tr[i][j]=t;
    }
  }

  // trace back
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  std::stack<std::pair<uint,uint> > st;
  st.push(std::make_pair(0, L-1));
  while (!st.empty())
  {
    const std::pair<uint,uint> p=st.top(); st.pop();
    const int i=p.first, j=p.second;
    switch (tr[i][j])
    {
      case 0:
        break;
      case 1:
        st.push(std::make_pair(i+1, j));
        break;
      case 2:
        st.push(std::make_pair(i, j-1));
        break;
      case 3:
        ss[i]=j;
        st.push(std::make_pair(i+1, j-1));
        break;
      default:
        const int k=i+tr[i][j]-3;
        st.push(std::make_pair(i, k));
        st.push(std::make_pair(k+1, j));
        break;
    }
  }
  return dp[0][L-1];
#endif
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

void
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

  // precalculate the range for alignment
  std::vector< std::pair<uint,uint> > r;
  alignment_envelope(p_z, th_a_, r);

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
    s += decode_secondary_structure(p_x, q_x, x);
    s += decode_secondary_structure(p_y, q_y, y);
    s += decode_alignment(p_z, q_z, r, z);

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

void
Dusaf::
align(VU& ss, ALN& aln, int ch) const
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
    align_alignments(ss, aln, aln1, aln2);
  }
}

void
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
  align_alignments(ss, r, a[0], a[1]);
  std::swap(r, aln);
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
  int ch;
  const char* a_en = "CONTRAlign";
  const char* s_en = "Boltzmann";
  std::string arg_x;
  while ((ch=getopt(argc, argv, "hv:a:s:p:q:r:t:g:u:w:e:f:m:lbx:"))!=-1)
  {
    switch (ch)
    {
      case 'a':
        a_en = optarg;
        break;
      case 's':
        s_en = optarg;
        break;
      case 'p':
        n_pct_a_ = atoi(optarg);
        break;
      case 'q':
        n_pct_s_ = atoi(optarg);
        break;
      case 'r':
        n_refinement_ = static_cast<uint>(atoi(optarg));
        break;
      case 't':
        th_s_ = atof(optarg);
        break;
      case 'g':
        th_s_ = 1.0/(atof(optarg)+1.0);
        break;
      case 'u':
        th_a_ = atof(optarg);
        break;
      case 'w':
        w_ = atof(optarg);
        break;
      case 'e':
        eta0_ = atof(optarg);
        break;
      case 'm':
        t_max_ = atoi(optarg);
        break;
      case 'l':
        use_alifold_ = true;
        break;
      case 'b':
        use_bpscore_ = true;
        break;
      case 'x':
        arg_x = optarg;
        break;
      case 'v':
        verbose_ = atoi(optarg);
        break;
      case 'h': case '?': default:
        usage(argv[0]);
        exit(0);
        break;
    }
  }
  argc -= optind;
  argv += optind;

  if (strcasecmp(a_en, "CONTRAlign")==0)
    en_a_ = new CONTRAlign(th_a_);
  else if (strcasecmp(a_en, "ProbCons")==0)
    en_a_ = new ProbCons(th_a_);
  else if (strcasecmp(a_en, "PartAlign")==0)
    en_a_ = new PartAlign(th_a_, arg_x);

  if (strcasecmp(s_en, "Boltzmann")==0)
    en_s_ = new RNAfold(true, NULL, CUTOFF);
  else if (strcasecmp(s_en, "Vienna")==0)
    en_s_ = new RNAfold(false, NULL, CUTOFF);
  else if (strcasecmp(s_en, "CONTRAfold")==0)
    en_s_ = new CONTRAfold(CUTOFF);

  if (!en_a_ || !en_s_)
  {
    usage(argv[0]);
    exit(0);
  }

  return *this;
}

void
Dusaf::
print_matching_probability(std::ostream& os, const VVF& p) const
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

void
Dusaf::
print_basepairing_probability(std::ostream& os, const VVF& p) const
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

int
Dusaf::
run(const char* filename)
{
  // read sequences
  uint N=Fasta::load(fa_, filename);

  // calculate base-pairing probabilities
  bp_.resize(N);
  for (uint i=0; i!=N; ++i)
    en_s_->fold(fa_[i].seq(), bp_[i]);

  // calculate matching probabilities
  mp_.resize(N, std::vector<MP>(N));
  dist_.resize(N*(N-1)/2);
  for (uint i=0, k=0; i!=N-1; ++i)
  {
    for (uint j=i+1; j!=N; ++j, ++k)
    {
      float d = use_bpscore_ ?
        en_a_->align(fa_[i].seq(), fa_[j].seq(), bp_[i], bp_[j], mp_[i][j]) :
        en_a_->align(fa_[i].seq(), fa_[j].seq(), mp_[i][j]);
      transpose_mp(mp_[i][j], mp_[j][i], fa_[i].size(), fa_[j].size());
      dist_[k] = d/(std::min(fa_[i].size(), fa_[j].size()));
    }
  }
  for (uint i=0; i!=N; ++i)
  {
    mp_[i][i].resize(fa_[i].size());
    for (uint x=0; x!=fa_[i].size(); ++x)
      mp_[i][i][x].push_back(std::make_pair(x, 1.0f));
  }

  // probabilistic consistency tranformation for matching probability matrix
  for (uint i=0; i!=n_pct_a_; ++i)
    relax_matching_probability();

  // probabilistic consistency tranformation for base-pairing probabilitiy matrix
  for (uint i=0; i!=n_pct_s_; ++i)
    relax_basepairing_probability();
  
  // compute the guide tree
  build_tree();
  print_tree(std::cout, tree_.size()-1);
  std::cout << std::endl;

  // compute progressive alignments along with the guide tree
  VU ss;
  ALN aln;
  align(ss, aln, tree_.size()-1);
  
  // iterative refinement
  for (uint i=0; i!=n_refinement_; ++i)
    refine(ss, aln);

#if 0
  if (do_refolding_)
  {
    // compute the common secondary structures from the averaged base-pairing matrix
    const uint L=aln[0].second.size();
    VVF p;
    average_basepairing_probability(p, aln);
    for (uint i=0; i!=L-1; ++i)
      for (uint j=i+1; j!=L; ++j)
        p[i][j] -= th_s_;
    decode_secondary_structure(p, ss);
  }
#endif
  
  // output the alignment
  std::sort(aln.begin(), aln.end());
  output(std::cout, aln, ss);

  return 0;
}

int
main(int argc, char* argv[])
{
  Dusaf dusaf;
  dusaf.parse_options(argc, argv);
  return dusaf.run(argv[0]);
}
