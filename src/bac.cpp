#include"bac.h"

BranchAndCut::BranchAndCut(string input_file, string output_file) :
_input_file(input_file), _output_file(output_file)
{
	if (!read_file(input_file, _mat_a, _mat_b))
	{
		cout << "read case error" << endl;
	}
}


// read_file to matrix stored in row
bool BranchAndCut::read_file(const string file, Matf& mat_a, Matf& mat_b)
{

	std::fstream file_stream;
	std::stringstream iostring;
	string line;

	file_stream.open(file, fstream::in | fstream::app);

	std::getline(file_stream, line);
	iostring << line;
	iostring >> net_node_number_ >> edge_number_ >> consumer_node_number_;

	std::getline(file_stream, line); //跳过空行
	std::getline(file_stream, line);
	iostring.clear();
	iostring << line;
	iostring >> price_server_;

	std::getline(file_stream, line);//跳过空行

	double temp;
	multimap<int, Vectorf> store_in_map;

	for (int i = 0; i < edge_number_; i++)
	{
		Vectorf temp_vector;
		temp_vector.reserve(4);
		std::getline(file_stream, line);
		iostring.clear();
		iostring << line;
		while (iostring >> temp)
		{
			temp_vector.push_back(temp);
		}
		if (temp_vector[0] > temp_vector[1])
		{
			swap(temp_vector[0], temp_vector[1]);
		}
		//store_in_map[temp_vector[0]] = temp_vector;
		store_in_map.insert({ temp_vector[0], temp_vector }); //map是为了使得第一列（网络节点）递增排序
	}

	// 网络节点部分的数据（四列）
	mat_a.reserve(edge_number_);
	for (const auto &i : store_in_map)
	{
		mat_a.push_back(i.second);
	}

	std::getline(file_stream, line);  //跳过空行

	// 消费节点部分的数据（三列）
	mat_b.reserve(consumer_node_number_);
	for (int i = 0; i < consumer_node_number_; i++)
	{

		Vectorf temp_b_vector;
		temp_b_vector.reserve(3);

		std::getline(file_stream, line);
		iostring.clear();
		iostring << line;

		while (iostring >> temp)
			temp_b_vector.push_back(temp);

		mat_b.push_back(temp_b_vector);

	}
	return true;
}


// transform the matrix stored in row to colunme
Matf BranchAndCut::store_by_colume(const Matf& mat_a)
{
	Matf mat_colume_store;
	Vectorf a;

	int colume = mat_a[0].size();
	int row = mat_a.size();
	for (int i = 0; i < colume; i++)
	{
		mat_colume_store.push_back(a);

		for (int j = 0; j < row; j++)
		{
			mat_colume_store[i].push_back(mat_a[j][i]);

		}
	}
	return mat_colume_store;
}




// this function nonuse----------------------------------------------------------see construct_A
/*
Transform the original equation to the simplex standard input type

original question：
min(cx)
subject:
A(m*n)x>=b

m row structure：
1 all nodes : input-output+maxflow*x>=b
2 sum(x)>=min_server_number
3 x<=1
4 l<=lmax
注：4 has two plans.a: At the first,l's equations hasn't been added the model,
use a funtion that check whether l cross the border.
b:directly add l into the model,disadvantages is the number of rows too big.

After standard:
max(-c(n+m)x)
subject:
-A(m*(n+m))x<=-b
*/
bool BranchAndCut::simplex_standard(Vectorf& cj, Matf&A, Vectorf&b)
{

	if (!A.size() || !A[0].size())
	{
		cout << "A is empty matrix" << endl;
		return false;
	}
	int row = A.size(), colume = A[0].size();

	for (int i = 0; i < row; i++)
	{
		cj.push_back(0);
		b[i] = -b[i];

		for (int j = 0; j < colume; j++)
		{
			A[i][j] = -A[i][j];
		}
	}

	for (int i = 0; i < colume; i++)
		cj[i] = -cj[i];

	for (int i = 0; i < row; i++)
		for (int j = 0; j < row; j++)
			A[i].push_back(0);

	for (int i = 0; i < row; i++)
	{
		A[i][colume + i] = 1;
	}

	return true;
}
//--------------------------------------------------------------------------------------------end



// the three functions nonuse----------------------------------------------------------------
//delete A_'s some row
bool BranchAndCut::delete_one_row(Matf& A_, int row_)
{
	if (row_ >= 0)
	{
		A_.erase(A_.begin() + row_);
		return true;
	}
	return false;

}

//delete A_'s some colume
bool BranchAndCut::delete_one_colume(Matf& A_, int colume)
{
	if (colume >= 0)
	{
		for (int i = 0; i < A_.size(); i++)
		{
			A_[i].erase(A_[i].begin() + colume);
		}
		return true;
	}
	return false;

}

// 
bool BranchAndCut::fix(const Vectorf &fixed_x, Vectorf& cj, Matf& A_, Vectorf& b)
{
	int count_row = 0, count_colume = 0;

	for (int i = 0; i < fixed_x.size(); i++)
	{
		if (fixed_x[i] != -2)
		{
			if (fixed_x[i] == 0)
			{
				A_[i - count_row][i - count_row] = 0;
			}
			else if (fixed_x[i] == 1)
			{
				b[net_node_number_ - count_row]--;
				delete_one_row(A_, i - count_row);
				b.erase(b.begin() + i - count_row);
				count_row++;
			}

			delete_one_colume(A_, i - count_colume);
			cj.erase(cj.begin() + i - count_colume);

			++count_colume;
		}

	}
	return true;

}
//--------------------------------------------------------------------------------------------end




/*

computed_solve: Along with the iteration, some varibles were determined, so the size of "computed_solve" is 
less than or equal to the number of solves(n).

final_result: some determined varibles inside already, then add the "computed_solve" to obtain the final 
solve(x).

final_result.second: a vector, including -2(undetermined position), and 0(determined position)) or 
1(determined position)

*/
bool BranchAndCut::recovery_x(const Solve& computed_solve, Solve& final_result)
{

	int count = 0;
	Vectorf::iterator iter = final_result.second.begin();
	while ((iter = find(iter, final_result.second.end(), 1)) != final_result.second.end())
	{
		count++;
		iter++;
	}
	final_result.first = -computed_solve.first + count*price_server_;
	Vectorf& x = final_result.second;

	// final_result=determinde solve+computed_solve
	// final_result: the final solve
	int j = 0;
	for (auto &i : x)
	{
		if (i == -2)
		{
			i = computed_solve.second[j];
			j++;
		}
	}
	return true;
}



// this function nonsense---------------------------------------------------------------------
Mati BranchAndCut::fix_A_nonzero_colume(const Vectorf &fixed_x)
{
	int count_row = 0, count_colume = 0;
	Mati fixed_A_non = A_nonzero_index;
	int fix_number = 0;
	Vectori fixed_vector, _row_vector;

	for (int i = 0; i < fixed_x.size(); i++)
	{
		if (fixed_x[i] != -2)
		{
			if (fixed_x[i] == 1)
			{
				_row_vector.push_back(i);
			}
			fixed_vector.push_back(i);
			fix_number++;
		}
	}

	for (int i = 0; i < fixed_A_non.size(); ++i)
	{

		int count = 0;
		bool is_delete = 0;

		for (auto &j : fixed_vector)
		{
			if (j < fixed_A_non[i][0])
				++count;
			else if (j == fixed_A_non[i][0])
			{
				fixed_A_non[i].erase(fixed_A_non[i].begin());
				fixed_A_non[i][0] -= fix_number;
				is_delete = 1;
				break;
			}
			else if (j>fixed_A_non[i][0])
			{
				break;
			}

		}
		if (!is_delete)
		{
			fixed_A_non[i][0] -= count;
		}

		for (int index = 1; index < fixed_A_non[i].size(); index++)
		{
			fixed_A_non[i][index] -= fix_number;
		}
	}

	for (auto i : _row_vector)
	{
		fixed_A_non.erase(fixed_A_non.begin() + i - count_row);
		++count_row;
	}

	if (!_row_vector.empty())
	{
		int base_varible_coeff = net_node_number_ + 2 * edge_number_ - fix_number;
		for (auto &i : fixed_A_non)
		{
			*(i.end() - 1) = base_varible_coeff;
			++base_varible_coeff;
		}
	}

	return fixed_A_non;

}
//------------------------------------------------------------------------------------------end



// The second optimization point
/*

when some solves(varibles more exact) have been determined, then some rows and some columes should be deleted,
but we mark them as 0 instead of eraseing them for better efficiency.

cj_index: initialized 1. if 0, indicate this colume was deleted.
cb_index: initialized 1. if 0, indicate this row was deleted.

erase's disadvantage:
erase(): delete and copy_forward.  this operation is time-cosuming, and make last indexs in A invalid.

*/
int BranchAndCut::fixed(Vectorf& cj, const Matf& A, Vectorf& b, const Vectorf& fixed_vector,
	Mati& A_nonzero_mat, Vectori& cb_index, Vectori &cj_index)
{

	Vectorf index_0, index_1;
	int fixed_number = 0;

	int fixed_vector_size = fixed_vector.size();
	int varible_number = net_node_number_ + 2 * edge_number_;
	int col = 1 + 3 * net_node_number_ + 4 * edge_number_;
	int row = 1 + 2 * net_node_number_ + 2 * edge_number_;

	for (int i = 0; i < fixed_vector_size; i++)
	{
		if (fixed_vector[i] != -2)
		{
			if (fixed_vector[i] == 0)
			{
				cj_index[i] = 0;
			}
			if (fixed_vector[i] == 1)// 放置服务器时，该行约束恒成立，所以删掉该行
			{
				//if cj[i]colume is marked 0, the meaning is the colume is deleted.
				cj_index[i] = 0;

				//if i row is marked 0, the meaning is the row is deleted.
				cb_index[i] = 0;
				b[net_node_number_]++;  // 若该顶点放置服务器则第二类约束（一行）不等式右边发生变化（见construct_A函数）

				//because the row is deleted ,the row of base varible should also be deleted
				cj_index[varible_number + i] = 0;
			}
			fixed_number++;
		}
	}
	return fixed_number;

}


// select a base varible to be replaced out.(out_varible) 
// convert to nonbase varible
int BranchAndCut::get_out_basevarible(const Vectorf& b,
	const Matf& A,
	const int& row,
	const int& col,
	const Vectori& cb_index,
	const Vectori& cj_index)
{

	int col_index, row_index, row_basevarible;

	int count = 0;

	for (row_index = 0; row_index < row; row_index++)
	{
		if (count == 0 && cb_index[row_index] == 1)
		{
			row_basevarible = row_index;
			++count;
		}
		if (cb_index[row_index] == 1)
		{
			if (b[row_index] < b[row_basevarible])
				row_basevarible = row_index;
		}
	}
	//if every number in the b colume >0,return  optimal solution
	if (b[row_basevarible] >= 0.0)
		return -2;

	for (col_index = 0; col_index < col; col_index++)
	{
		if (cj_index[col_index])
		{
			if (A[row_basevarible][col_index] < 0.0)
				break;
		}
	}

	//if every number in the output base varible's row >0,no optimal solution
	if (col_index >= col)
		return -1;
	return row_basevarible;

}


// select a nonbase varible to be replaced in.(in_varible) 
// convert to base varible
int BranchAndCut::get_in_basevarible(const Vectorf& b,
	const Matf& A,
	const Vectorf& sigma,
	const int & out_varible,
	const int &row, const int &col,
	const Mati& A_nonzero,
	const Vectori& cj_index)
{
	int col_index, find(0);
	double min_sigma = 99999;
	int col_basevarible;

	int col_count = A_nonzero[out_varible].size();
	for (col_index = 0; col_index < col_count; col_index++)
	{
		if (cj_index[A_nonzero[out_varible][col_index]])
		{
			if (A[out_varible][A_nonzero[out_varible][col_index]] < 0)
			{
				find++;
				if ((find == 1) || (min_sigma > sigma[A_nonzero[out_varible][col_index]] / A[out_varible][A_nonzero[out_varible][col_index]]))
				{
					double sita = sigma[A_nonzero[out_varible][col_index]] / A[out_varible][A_nonzero[out_varible][col_index]];
					min_sigma = sita;
					col_basevarible = A_nonzero[out_varible][col_index];
				}
			}
		}

	}
	return col_basevarible;

}


// for debug
void BranchAndCut::printmid(const Vectorf& b,
	const Vectori& B,
	const Matf& A,
	const Vectorf& cj,
	const int& row, const int& col)
{

	cout << "\t\t";
	for (auto i : cj)
	{
		cout << i << "\t";
	}
	cout << endl;
	for (int row_index = 0; row_index < row; row_index++)
	{
		cout << B[row_index] << "\t" << b[row_index] << "\t";
		for (int col_index = 0; col_index < col; col_index++)
		{
			cout << A[row_index][col_index] << "\t";
		}
		cout << endl;
	}

}


//solve in the main funcion , please initialize it, and allocate the memory
void BranchAndCut::print_result(Solve& _solve,
	const Vectorf& cj,
	const Vectorf& cb,
	const Vectorf& b,
	const Vectori& B,
	const int& row, const int&col,
	const Vectori& cb_index,
	const Vectori& cj_index,
	const int fixed_number
	)
{
	int row_index, col_index;
	double opti = 0;

	double& optimum = _solve.first;

	Vectorf result(col, 0);

	Vectorf x(col - row - fixed_number, 0);

	for (row_index = 0; row_index < row; row_index++)
	{
		if (cb_index[row_index])
		{
			if (cj_index[B[row_index]])
			{
				result[B[row_index]] = b[row_index];
			}
		}
	}

	optimum = 0;
	for (int i = 0; i < col - row; i++)
	{
		if (cj_index[i])
		{
			optimum += result[i] * cj[i];
		}
	}

	int xcount = 0;
	for (int col_index = 0; col_index < col - row; col_index++)
	{
		if (cj_index[col_index])
		{
			x[xcount] = result[col_index];
			++xcount;
		}
	}
	_solve.second = x;

}


// dualsimplex (important)
int BranchAndCut::dualsimplex(Vectorf& cj,
	Matf& A,
	Vectorf& b,
	Solve &solve,
	bool ifprint,
	Mati& A_nonzero_index,
	const Vectori&cb_index,
	const Vectori&cj_index,
	const int fixed_number)
{
	int row = A.size();
	int col = A[0].size();
	Vectorf sigma(col), cb(row, 0);

	Vectori  B(row, 0);
	for (int i = 0; i < row; i++)
	{
		B[i] = col - row + i;
	}

	int res, out_basevarible, in_basevarible;
	double temp_b;

	sigma = cj;

	int count = 0;
	while (true)
	{
		++count;
		res = get_out_basevarible(b, A, row, col, cb_index, cj_index);

		if (res == -2)
		{
			print_result(solve, cj, cb, b, B, row, col, cb_index, cj_index, fixed_number);
			cout << "update " << count << endl;
			return true;
		}
		else if (res == -1)
		{

			cout << "update " << count << endl;
			solve.first = -1;
			cout << "无最优解！" << endl;
			return false;

		}
		else
		{
			out_basevarible = res;
			in_basevarible = get_in_basevarible(b, A, sigma, out_basevarible, row, col, A_nonzero_index, cj_index);
			//if (ifprint)
			//printmid();
			B[out_basevarible] = in_basevarible;
			cb[out_basevarible] = cj[in_basevarible];
			temp_b = b[out_basevarible] / A[out_basevarible][in_basevarible];

			//将出基变量所在行先进行变换，也就是所有元素除以A中出基变量所在行和入基变量所在列的元素
			b[out_basevarible] = temp_b;

			//保存出基行和入基列的主元素
			double principle_value = A[out_basevarible][in_basevarible];

			//for (int col_index = 0; col_index < col; col_index++)
			//	A[out_basevarible][col_index] /= principle_value;
			int col_size = A_nonzero_index[out_basevarible].size();
			for (int col_index = 0; col_index <col_size; col_index++)
			{
				if (cj_index[A_nonzero_index[out_basevarible][col_index]])
				{
					A[out_basevarible][A_nonzero_index[out_basevarible][col_index]] /= principle_value;
				}
			}


			//计算出基变量外的行元素
			for (int row_index = 0; row_index < row; row_index++)
			{

				if (cb_index[row_index])
				{
					//计算出基变量所在行
					if (row_index != out_basevarible)
					{
						//	b[row_index] -= temp_b*A[row_index][in_basevarible];
						//	double out_basevarible_scale = A[row_index][in_basevarible];
						//	for (int col_index = 0; col_index < col; col_index++)
						//	A[row_index][col_index] -= A[out_basevarible][col_index] * out_basevarible_scale;

						double out_basevarible_scale = A[row_index][in_basevarible];
						//b[row_index] -= temp_b*A[row_index][in_basevarible];
						if (out_basevarible_scale != 0)
						{
							//更新b
							b[row_index] -= temp_b*out_basevarible_scale;

							int out_basevarible_number = A_nonzero_index[out_basevarible].size();
							for (int col_index = 0; col_index <out_basevarible_number; col_index++)
							{
								if (cj_index[A_nonzero_index[out_basevarible][col_index]])
								{
									double temp = A[row_index][A_nonzero_index[out_basevarible][col_index]] - A[out_basevarible][A_nonzero_index[out_basevarible][col_index]] * out_basevarible_scale;

									if (temp != 0)
									{
										if (A[row_index][A_nonzero_index[out_basevarible][col_index]] == 0)
										{
											A_nonzero_index[row_index].push_back(A_nonzero_index[out_basevarible][col_index]);
										}
									}

									A[row_index][A_nonzero_index[out_basevarible][col_index]] = temp;

									if (temp == 0)
									{
										Vectori::iterator index_iter = find(A_nonzero_index[row_index].begin(), A_nonzero_index[row_index].end(), A_nonzero_index[out_basevarible][col_index]);
										A_nonzero_index[row_index].erase(index_iter);
									}
								}
							}
						}
					}
				}
			}

			//计算判别数
			double coeff_sigma = sigma[in_basevarible];
			// int sigma_count = A_nonzero_index[out_basevarible].size();
			for (int col_index = 0; col_index < col_size; col_index++)
			{
				if (cj_index[A_nonzero_index[out_basevarible][col_index]])
				{
					sigma[A_nonzero_index[out_basevarible][col_index]] -= coeff_sigma*A[out_basevarible][A_nonzero_index[out_basevarible][col_index]];
				}
			}

		}

	}

	return true;
}


// note: this Matf_a should be stored in colume, compute the max out-flow of each net_node.(use the capacity 
// of each adjacent edge)
Vectorf BranchAndCut::one_node_max_flow(const Matf& mat_a)
{
	Vectorf node(net_node_number_, 0);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j<edge_number_; j++)
		{
			node[mat_a[i][j]] += mat_a[2][j];
		}
	}
	return node;
}


// find the minimum number of servers. the out-flow of net_nodes with server must be gerater than or equal to 
// the sum of demand(all comsumer_node).
int BranchAndCut::min_server_number(const Vectorf& node, const double sum_demand)
{
	double sum = 0;
	Vectorf sort_node = node;
	int server_number = 0;

	sort(sort_node.begin(), sort_node.end());

	for (int i = sort_node.size() - 1; i >= 0; i--)
	{
		if (sum < sum_demand)
		{
			sum += sort_node[i];
			server_number++;
		}
		else
			break;
	}
	return server_number;

}


//
bool BranchAndCut::is_feasible_solution(const Vectorf& x)
{
	for (auto i : x)
	{
		if (!is_integer(i))
			return false;
	}
	return true;
}


// compute A
void BranchAndCut::constraint_A(Matf& A,
	Vectori& row_index,
	const int& index,
	const Vectorf& array1,
	const Vectorf& array2,
	const Vectorf& node_max_flow,
	const double net_node_num,
	const double edge_num
	)
{

	int cnt = count(array1.cbegin(), array1.cend(), index - 1);
	if (cnt != 0)
	{
		Vectorf::const_iterator beg = find(array1.cbegin(), array1.cend(), index - 1);
		Vectorf::difference_type idx = beg - array1.cbegin();
		for (int i = 0; i < cnt; ++i)
		{
			A[index - 1][idx + i + net_node_num] = -1;
			A[index - 1][idx + i + edge_num + net_node_num] = 1;

			row_index.push_back(idx + i + net_node_num);
			row_index.push_back(idx + i + edge_num + net_node_num);
		}
		int cnt2 = count(array2.cbegin(), array2.cend(), index - 1);
		if (cnt2 != 0)
		{
			Vectorf::const_iterator start_beg = array2.cbegin();
			while (cnt2)
			{
				Vectorf::const_iterator beg2 = find(start_beg, array2.cend(), index - 1);
				Vectorf::difference_type idx2 = beg2 - array2.cbegin();
				A[index - 1][idx2 + edge_num + net_node_num] = -1;
				A[index - 1][idx2 + net_node_num] = 1;

				row_index.push_back(idx2 + edge_num + net_node_num);
				row_index.push_back(idx2 + net_node_num);
				start_beg = beg2 + 1;
				--cnt2;
			}

		}

	}
	else
	{
		int cnt3 = count(array2.cbegin(), array2.cend(), index - 1);
		Vectorf::const_iterator start_beg3 = array2.cbegin();
		while (cnt3)
		{
			Vectorf::const_iterator beg3 = find(start_beg3, array2.cend(), index - 1);
			Vectorf::difference_type idx3 = beg3 - array2.cbegin();
			A[index - 1][idx3 + edge_num + net_node_num] = -1;
			A[index - 1][idx3 + net_node_num] = 1;

			row_index.push_back(idx3 + edge_num + net_node_num);
			row_index.push_back(idx3 + net_node_num);

			start_beg3 = beg3 + 1;
			--cnt3;
		}
	}
	A[index - 1][index - 1] = node_max_flow[index - 1];
}


// really compute A, and compute b, and compute A'(record the nonzero colume)
void BranchAndCut::construct_A(Vectorf& cj,
	Matf& A,
	Vectorf& b,
	const Vectorf& node_max_flow,
	const double& min_servers,
	const int constraint_num,
	const Matf& mat_a_colume,
	const Matf& mat_b_colume,
	Mati& A_nonzero_index)
{
	int base_varible_row = 0;

	// 标准A的初始化  0矩阵
	//放外面Matf A(constraint_num);  A已经初始化
	for (int k = 0; k < constraint_num; ++k)
	{
		Vectorf A_i(1 + 3 * net_node_number_ + 4 * edge_number_, 0); //A(m*(n+m))
		A[k] = A_i;
	}

	//1:in-out+w*i>=b
	//standard:-in+out-w*i<=-b
	for (int i = 1; i <= net_node_number_; ++i)
	{
		Vectori row_index;
		row_index.reserve(50);  //lsc
		//  A[1] 开始求的  in-out+max*xi>=0(bi)
		constraint_A(A, row_index, i, mat_a_colume[0], mat_a_colume[1], node_max_flow,
			net_node_number_, edge_number_);

		//构造A非0值的标号，in-out+max*x>bm行
		row_index.push_back(i - 1);
		base_varible_row = net_node_number_ + 2 * edge_number_ + i - 1;
		row_index.push_back(base_varible_row);
		sort(row_index.begin(), row_index.end());
		//A_nonzero_index外面初始化。
		//A_nonzero_index.reserve(constraint_num);
		//A_nonzero_index.push_back(row_index); //A_nonzero_index知道大小吗？知道的话还可以改
		A_nonzero_index[i - 1] = row_index;

	}

	//2:sum(x)>=min_servers
	//standard:-sum(x)<=-min_servers
	Vectorf temp_A_line(net_node_number_, -1);
	copy(temp_A_line.cbegin(), temp_A_line.cend(), A[net_node_number_].begin());

	//for (int i = 0; i < net_node_number_; ++i)
	//{
	//	A[net_node_number_][i] = 1;   // sum(xi)>=1;
	//}
	//构造A非0值的标号，sum(x)<min_servers行
	Vectori sigma_row_index(net_node_number_ + 1, 0);
	for (int i = 0; i < net_node_number_; i++)
		sigma_row_index[i] = i;
	sigma_row_index[net_node_number_] = ++base_varible_row;
	A_nonzero_index[net_node_number_] = sigma_row_index;

	//3:x<=1
	int node_edge_index = 0;
	for (int i = 0; i < net_node_number_; i++)
	{

		A[net_node_number_ + 1 + i][i] = 1;
		cj[i] = (-price_server_);

		//构造A非0值的标号，x<1 行
		Vectori x_row_index;
		x_row_index.push_back(node_edge_index++);
		x_row_index.push_back(++base_varible_row);
		A_nonzero_index[net_node_number_ + i + 1] = (x_row_index);
	}

	//4:l<=lmax (edge_number_ *2 rows)
	int row_l = 2 * net_node_number_ + 1;

	for (int i = 0; i < edge_number_; i++)
	{
		A[2 * net_node_number_ + 1 + i][net_node_number_ + i] = 1;
		cj[net_node_number_ + i] = -mat_a_colume[3][i];

		//构造A非0值的标号，l<lmax行
		Vectori lout_row_index;
		lout_row_index.push_back(node_edge_index++);
		lout_row_index.push_back(++base_varible_row);

		A_nonzero_index[row_l + i] = lout_row_index;
	}

	row_l = 2 * net_node_number_ + 1 + edge_number_;
	for (int i = 0; i < edge_number_; i++)
	{
		A[2 * net_node_number_ + 1 + edge_number_ + i][edge_number_ + net_node_number_ + i] = 1;
		cj[net_node_number_ + edge_number_ + i] = -mat_a_colume[3][i];

		//构造A非0值的标号，l<lmax行
		Vectori lin_row_index;
		lin_row_index.push_back(node_edge_index++);
		lin_row_index.push_back(++base_varible_row);
		A_nonzero_index[row_l + i] = lin_row_index;
	}

	//5:compute b
	for (int i = 0; i < consumer_node_number_; i++)
	{
		b[mat_b_colume[1][i]] = -mat_b_colume[2][i];
	}
	b[net_node_number_] = -min_servers;
	for (int i = 0; i < net_node_number_; i++)
	{
		b[net_node_number_ + 1 + i] = 1;
	}
	for (int i = 0; i < edge_number_; i++)
	{
		b[2 * net_node_number_ + 1 + i] = mat_a_colume[2][i];
	}

	for (int i = 0; i < edge_number_; i++)
	{
		b[2 * net_node_number_ + 1 + i + edge_number_] = mat_a_colume[2][i];
	}

	//....
	for (int i = 0; i < net_node_number_; i++)
	{
		int size = A_nonzero_index[i].size();
		for (int j = 0; j < size; j++)
		{
			A[i][A_nonzero_index[i][j]] = -A[i][A_nonzero_index[i][j]];
		}
	}

	// add slack variable (A unit array)
	int row = 2 * net_node_number_ + 1 + 2 * edge_number_; //m
	int varible_number = net_node_number_ + 2 * edge_number_; //n
	for (int i = 0; i < row; i++)
	{
		A[i][varible_number++] = 1;
	}
}


// output path------------------------------------------------------------------------
// for a net_node, compute it's out_edges and in_edges
vector<std::pair<double, vector<double>>> BranchAndCut::net_valid_edge(const double& index, 
	const vector<double>& array1,
	const vector<double>& array2, 
	const double& edge_num,
	const vector<double>& l)
{

	vector<double> edge_out;
	vector<double> edge_in;
	std::pair<double, vector<double>> edge_out_pair;
	std::pair<double, vector<double>> edge_in_pair;
	vector<std::pair<double, vector<double>>> result;


	double cnt = count(array1.cbegin(), array1.cend(), index - 1);
	if (cnt != 0)
	{
		vector<double>::const_iterator beg = find(array1.cbegin(), array1.cend(), index - 1);
		vector<double>::difference_type idx = beg - array1.cbegin();
		for (double i = 0; i < cnt; ++i)
		{
			//output += l[idx + i];

			//input += l[idx + i + edge_num];
			if (l[idx + i] != 0)
				edge_out.push_back(idx + i);
			if (l[idx + i + edge_num] != 0)
				edge_in.push_back(idx + i + edge_num);
		}
		double cnt2 = count(array2.cbegin(), array2.cend(), index - 1);
		if (cnt2 != 0)
		{
			vector<double>::const_iterator start_beg = array2.cbegin();
			while (cnt2)
			{
				vector<double>::const_iterator beg2 = find(start_beg, array2.cend(), index - 1);
				vector<double>::difference_type idx2 = beg2 - array2.cbegin();
				//output += l[idx2 + edge_num];
				//input += l[idx2];
				if (l[idx2 + edge_num] != 0)
					edge_out.push_back(idx2 + edge_num);
				if (l[idx2] != 0)
					edge_in.push_back(idx2);
				start_beg = beg2 + 1;
				--cnt2;
			}
		}
	}
	else
	{
		double cnt3 = count(array2.cbegin(), array2.cend(), index - 1);
		vector<double>::const_iterator start_beg3 = array2.cbegin();
		while (cnt3)
		{
			vector<double>::const_iterator beg3 = find(start_beg3, array2.cend(), index - 1);
			vector<double>::difference_type idx3 = beg3 - array2.cbegin();
			//output += l[idx3 + edge_num];
			//input += l[idx3];
			if (l[idx3 + edge_num] != 0)
				edge_out.push_back(idx3 + edge_num);
			if (l[idx3] != 0)
				edge_in.push_back(idx3);
			start_beg3 = beg3 + 1;
			--cnt3;
		}
	}

	edge_out_pair = std::make_pair(index - 1, edge_out);
	edge_in_pair = std::make_pair(index - 1, edge_in);
	result.push_back(edge_out_pair);
	result.push_back(edge_in_pair);
	return result;
}


// note: from the consumer_node to server
/*
for a consumer_node, find the paths to servers and store them.(path may more than one, multiple 
servers server this consumer_node)
*/
void BranchAndCut::func(vector<std::pair<double, vector<double>>> net_edge, double size, vector<double> vetr,
	double consume_demd_temp, double idx, double cons_net_node_idx_init,
	vector<double> net_node_idx,
	double consume_demd,
	const vector<double>& array1,
	const vector<double>& array2, const double& edge_num, vector<double>&l,
	const string& outfile, vector<double> iterater_vector, double last_edge_index,
	int& path_count, Matf& store_path)
{
	//double cons_net_node_idx = 0;

	net_edge =
		net_valid_edge(cons_net_node_idx_init + 1, array1,
		array2, edge_num, l);           // 统计该节点不为零的出边和入边

	size = net_edge[1].second.size();
	vetr = net_edge[1].second;

	if (vetr.empty())                // vetr 存储入边
	{
		// todo...
		//consume_demd = std::min(consume_demd, consume_demd_temp);
		net_node_idx.push_back(cons_net_node_idx_init);
		iterater_vector.push_back(last_edge_index);

		double min = 999;
		for (auto i : iterater_vector)
		{
			if (i != 999)
			{
				if (min >= l[i])
					min = l[i];
			}
		}
		if (min != 0)
		{

			if (min == 999)
				net_node_idx.push_back(consume_demd);
			//outfile_stream << consume_demd << endl;
			else
			{
				net_node_idx.push_back(min);
				//outfile_stream << min << endl;

			}
			path_count++;

			//outfile_stream.close();

			store_path.push_back(net_node_idx);
			for (auto i : iterater_vector)
			{
				if (i != 999)
					l[i] -= min;
			}
		}
	}
	else
	{
		//consume_demd = std::min(consume_demd, consume_demd_temp);
		net_node_idx.push_back(cons_net_node_idx_init);

		iterater_vector.push_back(last_edge_index);
		for (idx = 0; idx < size; ++idx)   // size=每个节点的入边数目 net_edge[1].second.size()
		{
			if (vetr[idx] >= edge_num)
				cons_net_node_idx_init = array2[vetr[idx] - edge_num];   // 消费节点
			else
				cons_net_node_idx_init = array1[vetr[idx]];   // 消费节点

			consume_demd_temp = l[vetr[idx]];
			last_edge_index = vetr[idx];

			func(net_edge, size, vetr, consume_demd_temp, idx, cons_net_node_idx_init,
				net_node_idx, consume_demd, array1, array2,
				edge_num, l, outfile, iterater_vector, last_edge_index, path_count, store_path);
		}
	}
}


// for all consumer_nodes, print the paths(including all net_nodes in this path)
void BranchAndCut::write(const Vectorf& solve, 
	const string& _output_file,
	const Matf& mat_a_column,
	const Matf& mat_b_column)
{
	vector<double> x, l;
	x.reserve(net_node_number_);
	l.reserve(edge_number_ * 2);

	if (solve.empty())
		exit(0);

	copy(solve.cbegin(), solve.cbegin() + net_node_number_, std::back_inserter(x));
	copy(solve.cbegin() + net_node_number_, solve.cend(), std::back_inserter(l));


	int path_count = 0; //输出总路径数
	Matf store_path;

	vector<double> iterator_vector;

	for (double i = 0; i < consumer_node_number_; ++i)
	{
		double cons_net_node_idx_init = mat_b_column[1][i];   // 消费节点

		vector<std::pair<double, vector<double>>> net_edge =
			net_valid_edge(cons_net_node_idx_init + 1, mat_a_column[0],
			mat_a_column[1], edge_number_, l);           // 统计该节点不为零的出边和入边

		double size = net_edge[1].second.size();
		vector<double> vetr = net_edge[1].second;
		double consume_demd = mat_b_column[2][i];
		vector<double> net_node_idx;
		vector<double> array1 = mat_a_column[0];
		vector<double> array2 = mat_a_column[1];

		net_node_idx.push_back(i);

		func(net_edge, size, vetr, consume_demd, 0, cons_net_node_idx_init, net_node_idx,
			consume_demd,
			array1, array2, edge_number_, l, _output_file, iterator_vector, 999, path_count, store_path);


	}

	ofstream file_stream(_output_file);
	file_stream << path_count << endl;
	file_stream << endl;

	for (const auto&i : store_path)
	{
		std::for_each(i.crbegin() + 1, i.crend(),
			[&file_stream](const double& elem){file_stream << elem << " "; });

		file_stream << *(i.end() - 1) << endl;
	}

}

