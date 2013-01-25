#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include <lemon/network_simplex.h>

using namespace lemon;
using namespace std;

typedef NetworkSimplex<ListDigraph, int, double> AlgType;

void apply_lambda(double lambda, int r, int c,
    ListDigraph::ArcMap<double>* cost,
    const vector<vector<vector<ListDigraph::Arc> > >& colarcs);
int extract_emd_cost(int r, int c, const AlgType& alg,
    const vector<vector<vector<ListDigraph::Arc> > >& colarcs,
    void (*output_function)(const char*));
double extract_amp_sum(int r, int c, const AlgType& alg,
    const vector<vector<double> >& a,
    const vector<vector<ListDigraph::Arc> >& nodearcs,
    void (*output_function)(const char*));

void emd_flow(
    const vector<vector<double> >& a,
    int k,
    int emd_bound,
    vector<vector<bool> >* result,
    int* emd_cost,
    double* amp_sum,
    double* final_lambda,
    void (*output_function)(const char*),
    bool verbose) {

  clock_t total_time_begin = clock();

  const int kOutputBufferSize = 1000;
  char output_buffer[kOutputBufferSize];

  // number of rows
  int r = a.size();
  // number of columns
  int c = a[0].size();

  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "r = %d,  c = %d,  k = %d,  "
        "emd_bound = %d\n", r, c, k, emd_bound);
    output_function(output_buffer);
  }

  clock_t graph_construction_time_begin = clock();
  
  // nodes corresponding to the matrix entries
  std::vector<std::vector<ListDigraph::Node> > innode;
  std::vector<std::vector<ListDigraph::Node> > outnode;
  // arcs from innodes to outnodes
  std::vector<std::vector<ListDigraph::Arc> > nodearcs;
  // arcs corresponding to the EMD
  std::vector<std::vector<std::vector<ListDigraph::Arc> > > colarcs;
  // graph
  ListDigraph g;
  ListDigraph::ArcMap<int > capacity(g);
  ListDigraph::ArcMap<double> cost(g);

  // source
  ListDigraph::Node s = g.addNode();
  // sink
  ListDigraph::Node t = g.addNode();

  // create two nodes for each entry in the matrix
  innode.resize(r);
  outnode.resize(r);
  for (int ii = 0; ii < r; ++ii) {
    innode[ii].resize(c);
    outnode[ii].resize(c);
    for (int jj = 0; jj < c; ++jj) {
      innode[ii][jj] = g.addNode();
      outnode[ii][jj] = g.addNode();
    }
  }

  // add arcs from innodes to outnodes
  nodearcs.resize(r);
  for (int ii = 0; ii < r; ++ii) {
    nodearcs[ii].resize(c);
    for (int jj = 0; jj < c; ++jj) {
      nodearcs[ii][jj] = g.addArc(innode[ii][jj], outnode[ii][jj]);
      cost[nodearcs[ii][jj]] = - abs(a[ii][jj]);
      capacity[nodearcs[ii][jj]] = 1;
    }
  }

  // add arcs from source to column 1
  for (int ii = 0; ii < r; ++ii) {
    ListDigraph::Arc a = g.addArc(s, innode[ii][0]);
    cost[a] = 0;
    capacity[a] = 1;
  }

  // add arcs from column c to sink
  for (int ii = 0; ii < r; ++ii) {
    ListDigraph::Arc a = g.addArc(outnode[ii][c - 1], t);
    cost[a] = 0;
    capacity[a] = 1;
  }

  // add arcs between columns
  colarcs.resize(r);
  for (int row = 0; row < r; ++row) {
    colarcs[row].resize(c - 1);
    for (int col = 0; col < c - 1; ++col) {
      colarcs[row][col].resize(r);
      for (int dest = 0; dest < r; ++dest) {
        colarcs[row][col][dest] =
            g.addArc(outnode[row][col], innode[dest][col + 1]);
        cost[colarcs[row][col][dest]] = abs(row - dest);
        capacity[colarcs[row][col][dest]] = 1;
      }
    }
  }

  clock_t graph_construction_time = clock() - graph_construction_time_begin;

  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "The graph has %d nodes and %d "
        "arcs.\n", countNodes(g), countArcs(g));
    output_function(output_buffer);
    snprintf(output_buffer, kOutputBufferSize, "Total construction time: %lf "
        "s\n ", static_cast<double>(graph_construction_time) / CLOCKS_PER_SEC);
    output_function(output_buffer);

  }

  // cost scaling
  AlgType alg(g);
  alg.upperMap(capacity);
  alg.stSupply(s, t, k);

  // make lambda larger until we find a solution that fits into the EMD budget
  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize,
        "Finding large enough value of lambda ...\n");
    output_function(output_buffer);
  }

  double lambda_high = 0.01;
  while (true) {
    apply_lambda(lambda_high, r, c, &cost, colarcs);
    alg.costMap(cost);

    alg.run();
    int cur_emd_cost = extract_emd_cost(r, c, alg, colarcs, output_function);

    if (verbose) {
      snprintf(output_buffer, kOutputBufferSize, "l: %lf  emd: %d\n",
          lambda_high, cur_emd_cost);
      output_function(output_buffer);
    }

    if (cur_emd_cost <= emd_bound) {
      break;
    } else {
      lambda_high = lambda_high * 2;
    }
  }

  // binary search on lambda
  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Binary search on lambda ...\n");
    output_function(output_buffer);
  }

  double lambda_low = 0;
  double lambda_eps = 0.00001;
  while(lambda_high - lambda_low > lambda_eps) {
    double cur_lambda = (lambda_high + lambda_low) / 2;
    apply_lambda(cur_lambda, r, c, &cost, colarcs);
    alg.costMap(cost);

    alg.run();
    int cur_emd_cost = extract_emd_cost(r, c, alg, colarcs, output_function);

    if (verbose) {
      snprintf(output_buffer, kOutputBufferSize, "l_cur: %lf  (l_low: %lf, "
          "l_high: %lf)  emd: %d\n", cur_lambda, lambda_low, lambda_high,
          cur_emd_cost);
      output_function(output_buffer);
    }

    if (cur_emd_cost <= emd_bound) {
      lambda_high = cur_lambda;
    } else {
      lambda_low = cur_lambda;
    }
  }

  apply_lambda(lambda_high, r, c, &cost, colarcs);
  alg.costMap(cost);

  alg.run();
  *emd_cost = extract_emd_cost(r, c, alg, colarcs, output_function);
  *amp_sum = extract_amp_sum(r, c, alg, a, nodearcs, output_function);
  *final_lambda = lambda_high;

  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Final l: %lf, total cost: %lf, "
        "amp sum: %lf, EMD cost: %d\n", lambda_high, alg.totalCost(), *amp_sum,
        *emd_cost);
    output_function(output_buffer);
  }

  result->clear();
  result->resize(r);
  for (int ii = 0; ii < r; ++ii) {
    (*result)[ii].resize(c);
    for (int jj = 0; jj < c; ++jj) {
      (*result)[ii][jj] = (alg.flow(nodearcs[ii][jj]) > 0);
    }
  }

  clock_t total_time = clock() - total_time_begin;
  if (verbose) {
    snprintf(output_buffer, kOutputBufferSize, "Total time %lf s\n",
        static_cast<double>(total_time) / CLOCKS_PER_SEC);
    output_function(output_buffer);
  }

  return;
}

void apply_lambda(double lambda, int r, int c,
    ListDigraph::ArcMap<double>* cost,
    const vector<vector<vector<ListDigraph::Arc> > >& colarcs) {
  for (int row = 0; row < r; ++row) {
    for (int col = 0; col < c - 1; ++col) {
      for (int dest = 0; dest < r; ++dest) {
        (*cost)[colarcs[row][col][dest]] = lambda * abs(row - dest);
      }
    }
  }
}

int extract_emd_cost(int r, int c, const AlgType& alg, 
    const vector<vector<vector<ListDigraph::Arc> > >& colarcs,
    void (*output_function)(const char*)) {
  int emd_cost = 0;
  for (int row = 0; row < r; ++row) {
    for (int col = 0; col < c - 1; ++col) {
      for (int dest = 0; dest < r; ++dest) {
        if (alg.flow(colarcs[row][col][dest]) > 0) {
          emd_cost += abs(row - dest);
          if (alg.flow(colarcs[row][col][dest]) != 1) {
            output_function("ERROR: nonzero flow on a column edge is not 1.\n");
          }
        }
      }
    }
  }
  return emd_cost;
}

double extract_amp_sum(int r, int c, const AlgType& alg,
    const vector<vector<double> >& a,
    const vector<vector<ListDigraph::Arc> >& nodearcs,
    void (*output_function)(const char*)) {
  double amp_sum = 0;
  for (int row = 0; row < r; ++row) {
    for (int col = 0; col < c; ++col) {
      if (alg.flow(nodearcs[row][col]) > 0) {
        amp_sum += abs(a[row][col]);
        if (alg.flow(nodearcs[row][col]) != 1) {
          output_function("ERROR: nonzero flow on a node edge is not 1.\n");
        }
      }
    }
  }
  return amp_sum;
}
