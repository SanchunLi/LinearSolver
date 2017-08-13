#ifndef BAC_H
#define BAC_H

#include<iostream>
#include<fstream>
#include<sstream>

#include<vector>
#include<list>
#include<map>
#include<string>

#include<memory>
#include<algorithm>
#include<functional>
#include<iterator>
#include<ctime>

//store the input data
typedef std::vector<int> Vectori;
typedef std::vector<std::vector<int>> Mati;
typedef std::vector<double> Vectorf;
typedef std::vector<std::vector<double>> Matf;

typedef std::pair<double, Vectorf> Solve;

using namespace std;

class BranchAndCut
{
public:

	BranchAndCut(string input_file, string output_file);

	//read_file to matrix stored in row ,write file to txt in style
	bool read_file(const string file, Matf& mat_a, Matf& mat_b);

	/*************************************************************************************/


	//
	void constraint_A(Matf& A,
		Vectori& row_index,
		const int& index,
		const Vectorf& array1,
		const Vectorf& array2,
		const Vectorf& node_max_flow,
		const double net_node_num,
		const double edge_num
		);

	//Ax=b,construct the A matrix, and b, and A'(nonzero colume)
	void construct_A(Vectorf& cj,
		Matf& A,
		Vectorf& b,
		const Vectorf& node_max_flow,
		const double& min_servers,
		const int constraint_num,
		const Matf& mat_a_colume,
		const Matf& mat_b_colume,
		Mati& A_nonzero_index);

	//check whether the x is integer
	bool is_integer(double x)
	{ 
		return (((x - floor(x)) <= 0.00001 ? 1 : 0) || ((x - floor(x)) >= 0.99999 ? 1 : 0));
	}

	//make A cj b to standard simplex style
	bool simplex_standard(Vectorf& cj, Matf& A, Vectorf &b);

	// transform the matrix stored in rows to colunme
	Matf store_by_colume(const Matf&);

	// the three functions nonuse
	bool delete_one_row(Matf&, int);
	bool delete_one_colume(Matf&, int);
	bool fix(const Vectorf&, Vectorf&, Matf&, Vectorf&);
	
	//fix some value in x, and recovery x from the computed x
	bool recovery_x(const Solve&, Solve&);

	// 
	Mati fix_A_nonzero_colume(const Vectorf &fixed_x);

	// when some varibles were determined, mark the deleted rows and columes, not really erase
	int fixed(Vectorf& cj, const Matf& A, Vectorf& b, const Vectorf& fixed_vector,
		Mati&, Vectori&, Vectori&);

	//caution : the input Matf_a , should be stored in colume, and then compute every node's max flow  
	Vectorf one_node_max_flow(const Matf&);

	//the minimum server's number, the number may not a feasible solution
	int min_server_number(const Vectorf& node, const double sum_demand);

	//check whether the soulution is feasible
	bool is_feasible_solution(const Vectorf& x);
	
	/*************************************************************************************/


	//Because the linear equation's structure (object is min, constrains <= someone),
	//we use the dual simplex method to get the optimal solution.
	//functions below are called by dual simplex method.

	//get the out base varible index
	int get_out_basevarible(const Vectorf& b,
		const Matf& A,
		const int& row,
		const int& col,
		const Vectori&,
		const Vectori&);

	//get the in base varible index 
	int get_in_basevarible(const Vectorf& b,
		const Matf& A,
		const Vectorf& sigma,
		const int & out_varible,
		const int &row, const int &col,
		const Mati& A_nonzero,
		const Vectori& cj_index);

	//print the midium process 
	void printmid(const Vectorf& b,
		const Vectori& B,
		const Matf& A,
		const Vectorf& cj,
		const int& row, const int& col);

	//print the result
	void print_result(Solve& _solve,
		const Vectorf& cj,
		const Vectorf& cb,
		const Vectorf& b,
		const Vectori& B,
		const int& row, const int&col,
		const Vectori& cb_index,
		const Vectori& cj_index,
		const int fixed_number
		);

	//dual simplex method
	int dualsimplex(Vectorf& cj_,
		Matf& A_,
		Vectorf& b_,
		Solve &solve,
		bool ifprint,
		Mati& A_nonzero_index,
		const Vectori&,
		const Vectori&,
		const int fixed_number
		);

	/*************************************************************************************/


	//get valid edges for each net_node
	vector<std::pair<double, vector<double>>> net_valid_edge(const double& index, const vector<double>& array1,
		const vector<double>& array2, const double& edge_num,
		const vector<double>& l);

	// get output path from solves
	void func(vector<std::pair<double, vector<double>>> net_edge, double size, vector<double> vetr,
		double consume_demd_temp, double idx, double cons_net_node_idx_init,
		vector<double> net_node_idx,
		double consume_demd,
		const vector<double>& array1,
		const vector<double>& array2, const double& edge_num, vector<double>&l,
		const string& outfile, vector<double> iterater_vector, double last_edge_index,
		int& path_count, Matf& store_path);

	// write the path to output_file
	void write(const Vectorf& solve, const string& _output_file,
		const Matf& mat_a_column,
		const Matf& mat_b_column);

	/*************************************************************************************/


	//data member
	int edge_number_, net_node_number_, consumer_node_number_;
	double price_server_;
	Matf _A, _mat_a, _mat_b;
	Vectorf _cj, _b;
	string _input_file, _output_file;
	Mati A_nonzero_index;
	Vectorf _cj_index, _cb_index;

	/*************************************************************************************/
};

#endif
