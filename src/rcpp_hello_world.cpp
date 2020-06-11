// [[Rcpp::plugins("cpp17")]]

#include <Rcpp.h>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "node.hpp"

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
string ReadParentVector(string file_path, size_t mutation_count) {
  Rcout << "Reading...\n";

  ifstream dat_file (file_path);

  NumericVector vec(mutation_count);
  string line;
  // Optimal parent vector is on line 2 * mutation_count + 4.
  size_t line_no = 2 * mutation_count + 4;
  for (size_t i = 0; i < line_no; i++) {
    getline(dat_file, line);
  }
  Rcout << line_no << ": " << line << "\n";
  dat_file.close();
  return line;
}

//' @export
// [[Rcpp::export]]
NumericVector GetChains(NumericVector parent_vector) {
  size_t node_count = parent_vector.size();
  vector<Node *> nodes;
  NumericVector chain_labels(node_count);
  for (size_t i = 0; i < node_count; i++) {
    nodes.push_back(new Node(i));
  }
  Node *root = new Node(node_count);
  nodes.push_back(root);
  for (size_t i = 0; i < node_count; i++) {
    size_t parent_idx = parent_vector[i];
    nodes[parent_idx]->add_child(nodes[i]);
    Rcout << "node: " << nodes[i]->get_id() << ", parent: " << parent_idx << "\n";
  }

  for (size_t i = 0; i <= node_count; i++) {
      Rcout << "node: " << nodes[i]->get_id() << ", children: [ ";
      for (auto child : nodes[i]->get_children()) {
          Rcout << child->get_id() << " ";
      }
      Rcout << "]\n";
  }

  size_t chain_id = 0;
  deque<Node *> queue;
  for (auto child : root->get_children()) {
      queue.push_back(child);
  }
  while (!queue.empty()) {
    auto node = queue.back();
    Rcout << node->get_id() << "\n";
    queue.pop_back();
    while (true) {
        chain_labels(node->get_id()) = chain_id;
        Rcout << "num children: " << node->get_children().size() << "\n";
        if (node->get_children().size() == 1) {
            node = node->get_children()[0];
        } else {
            break;
        }
    }
    chain_labels(node->get_id()) = chain_id;
    chain_id++;

    for (auto child_node : node->get_children()) {
        queue.push_back(child_node);
    }
  }

  return chain_labels;
}
